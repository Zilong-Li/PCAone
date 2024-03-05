#include "LD.hpp"

using namespace std;

MyVector calc_sds(const MyMatrix & X)
{
    // compute degree of freedom
    const int df = X.rows() - 1; // N-1
    return (X.array().square().colwise().sum() / df).sqrt();
}

// given a list of snps, find its index per chr in the original pos
// assume chromosomes are continuous
// TODO check duplicated POS
std::tuple<Int2D, Int2D> get_target_snp_idx(const std::string & filebim,
                                            const Int1D & pos,
                                            const Int1D & chr_pos_end,
                                            const std::vector<std::string> & chrs,
                                            bool header,
                                            Int1D colidx)
{
    Int1D t_pos, t_chr_pos_end, idx, bp;
    std::vector<std::string> t_chrs;
    get_snp_pos_bim(filebim, t_pos, t_chr_pos_end, t_chrs, header, colidx);
    UMapIntInt mpos;
    int c, s, e, p, i;
    Int2D ret(t_chrs.size());
    Int2D ret2(t_chrs.size());
    for(int tc = 0; tc < (int)t_chrs.size(); tc++)
    {
        for(c = 0; c < (int)chrs.size(); c++)
            if(chrs[c] == t_chrs[tc]) break;
        e = chr_pos_end[c];
        s = c > 0 ? chr_pos_end[c - 1] : 0;
        for(i = s; i <= e; i++) mpos[pos[i]] = i;
        e = t_chr_pos_end[tc];
        s = tc > 0 ? t_chr_pos_end[tc - 1] : 0;
        for(i = s; i <= e; i++)
        {
            p = t_pos[i];
            if(mpos.count(p))
            {
                bp.push_back(p);
                idx.push_back(mpos[p]);
            }
        }
        ret[tc] = idx;
        ret2[tc] = bp;
        bp.clear();
        idx.clear();
        mpos.clear();
    }
    return std::make_tuple(ret, ret2);
}

void calc_ld_metrics(const std::string & fileout,
                     const std::string & filebim,
                     MyMatrix & G,
                     const MyVector & F,
                     const Int1D & snp_pos,
                     const Int1D & chr_pos_end,
                     int ld_window_bp,
                     double r2_tol,
                     bool verbose = false)
{
    const bool pick_random_one = F.size() > 0 ? false : true;
    cao << tick.date() << "start ld pruning and pick_random_one =" << pick_random_one << std::endl;
    G.rowwise() -= G.colwise().mean(); // Centering first
#if defined(DEBUG)
    std::ofstream ofs_win(fileout + ".ld.window");
    ofs_win << "#window\tchr\tpos_start\tpos_end\tnsites" << std::endl;
#endif
    Int1D ws, we;
    int nsnp = snp_pos.size();
    int j{0}, c{0}, w{0}, nsites;
    for(int i = 0; i < nsnp; i++)
    {
        if(snp_pos[i] == snp_pos[chr_pos_end[c]])
        {
            c++;
            continue;
        }
        for(j = i; j <= chr_pos_end[c]; j++)
            if(snp_pos[j] - snp_pos[i] > ld_window_bp) break;
        nsites = j - i;
        ws.push_back(i); // start pos in the window
        we.push_back(nsites); // the number of sites
#if defined(DEBUG)
        ofs_win << w++ << "\t" << c + 1 << "\t" << snp_pos[i] << "\t" << snp_pos[j - 1] << "\t" << nsites
                << std::endl;
#endif
    }
    MyVector sds = 1.0 / calc_sds(G).array();
    ArrayXb keep = ArrayXb::Constant(G.cols(), true);
    const double df = 1.0 / (G.rows() - 1); // N-1
    for(w = 0; w < (int)ws.size(); w++)
    {
        int i = ws[w];
        if(!keep(i)) continue;
#pragma omp parallel for
        for(int j = 1; j < we[w]; j++)
        {
            int k = i + j;
            if(!keep(k)) continue;
            double r = (G.col(i).array() * G.col(k).array() * (sds(i) * sds(k))).sum() * df;
            if(r * r > r2_tol)
            {
                int o = k; // or i
                if(!pick_random_one) o = MAF(F(k)) > MAF(F(i)) ? i : k;
                keep(o) = false;
            }
        }
    }
    std::ifstream fin(filebim);
    if(!fin.is_open()) throw invalid_argument("can not open " + fileout + ".kept.bim");
    std::ofstream ofs_out(fileout + ".ld.prune.out");
    std::ofstream ofs_in(fileout + ".ld.prune.in");
    std::string line;
    int i = 0;
    while(getline(fin, line))
    {
        if(keep(i))
            ofs_in << line << std::endl;
        else
            ofs_out << line << std::endl;
        i++;
    }
}

Int1D valid_assoc_file(const std::string & fileassoc, const std::string & colnames)
{
    std::ifstream fin(fileassoc);
    if(!fin.is_open()) throw invalid_argument("can not open " + fileassoc);
    std::string line, sep{"\t"}, sep2{","};
    getline(fin, line);
    const auto fields = split_string(line, sep);
    std::vector<std::string> fields_users{"CHR", "BP", "P"};
    if(!colnames.empty()) fields_users = split_string(colnames, sep2);
    Int1D idx(3, -1);
    int j = 0;
    for(auto col : fields)
    {
        if(col == fields_users[0])
        {
            idx[0] = j;
        }
        else if(col == fields_users[1])
        {
            idx[1] = j;
        }
        else if(col == fields_users[2])
        {
            idx[2] = j;
        }
        j++;
    }
    j = 0;
    for(auto i : idx)
    {
        if(i < 0) cao.error("the assoc-like file has no " + fields_users[j] + " column");
        j++;
    }
    return idx;
}

std::vector<UMapIntDouble> map_index_snps(const std::string & fileassoc,
                                          const Int1D & colidx,
                                          double clump_p2)
{
    std::ifstream fin(fileassoc);
    if(!fin.is_open()) throw invalid_argument("can not open " + fileassoc);
    std::string line, chr_cur, chr_prev, sep{"\t"};
    getline(fin, line);
    vector<UMapIntDouble> vm;
    UMapIntDouble m;
    int c = 0, bp;
    double pval;
    while(getline(fin, line))
    {
        auto tokens = split_string(line, sep);
        chr_cur = tokens[colidx[0]];
        bp = std::stoi(tokens[colidx[1]]);
        pval = std::stod(tokens[colidx[2]]);
        if(pval <= clump_p2) m.insert({bp, pval});
        if(!chr_prev.empty() && chr_prev != chr_cur)
        {
            c++;
            vm.push_back(m);
            m.clear();
        }
        chr_prev = chr_cur;
    }
    vm.push_back(m); // add the last chr
    return vm;
}

std::vector<UMapIntString> map_assoc_file(const std::string & fileassoc, const Int1D & colidx)
{
    std::ifstream fin(fileassoc);
    if(!fin.is_open()) throw invalid_argument("can not open " + fileassoc);
    vector<UMapIntString> vm;
    UMapIntString m;
    int c = 0, bp;
    std::string line, chr_cur, chr_prev, sep{"\t"};
    getline(fin, line);
    m.insert({-1, line});
    while(getline(fin, line))
    {
        auto tokens = split_string(line, sep);
        chr_cur = tokens[colidx[0]];
        bp = std::stoi(tokens[colidx[1]]);
        m.insert({bp, line});
        if(!chr_prev.empty() && chr_prev != chr_cur)
        {
            c++;
            vm.push_back(m);
            m.clear();
        }
        chr_prev = chr_cur;
    }
    vm.push_back(m); // add the last chr
    return vm;
}

void calc_ld_clump(std::string fileout,
                   std::string fileassoc,
                   std::string colnames,
                   int clump_bp,
                   double clump_r2,
                   double clump_p1,
                   double clump_p2,
                   const MyMatrix & G,
                   const Int1D & snp_pos,
                   const Int1D & chr_pos_end,
                   const std::vector<std::string> & chrs)
{
    cao << tick.date() << "start do LD-based clumping -> p1=" << clump_p1 << ", p2=" << clump_p2
        << ", r2=" << clump_r2 << ", bp=" << clump_bp << ", assoc file:" + fileassoc << std::endl;
    auto colidx = valid_assoc_file(fileassoc, colnames); // 0: chr, 1: pos, 2: pvalue
    Int2D idx_per_chr, bp_per_chr;
    std::tie(idx_per_chr, bp_per_chr) =
        get_target_snp_idx(fileassoc, snp_pos, chr_pos_end, chrs, true, colidx);
    const auto pvals_per_chr = map_index_snps(fileassoc, colidx, clump_p2);
    const auto line_per_chr = map_assoc_file(fileassoc, colidx);
    // sort by pvalues and get new idx
    const MyVector sds = 1.0 / calc_sds(G).array();
    const double df = 1.0 / (G.rows() - 1); // N-1
#pragma omp parallel for
    for(int c = 0; c < (int)idx_per_chr.size(); c++)
    {
        const auto idx = idx_per_chr[c];
        const auto bp = bp_per_chr[c];
        auto mbp = vector2map(bp);
        auto lines = line_per_chr[c];
        std::ofstream ofs(fileout + ".clump.chr" + std::to_string(c + 1));
        ofs << line_per_chr[0].at(-1) << "\tSP2\n";
        // greedy clumping algorithm
        auto mpp = pvals_per_chr[c]; // key: pos, val: pval
        Double1D pp;
        Int1D ps;
        for(auto it = mpp.begin(); it != mpp.end(); it++)
        {
            if(it->second <= clump_p1)
            {
                ps.push_back(it->first);
                pp.push_back(it->second);
            }
        }
        for(auto i : sortidx(pp))
        { // snps sorted by p value
            int p = ps[i], p2;
            if(mpp.count(p) == 0)
                continue; // if snps with pval < clump_p1 are already clumped with previous snps
            Int1D clumped;
            bool backward = true;
            const size_t j = mbp[p]; // j:current
            size_t k = j; //  k:forward or backward
            while(true)
            {
                if(backward)
                {
                    --k;
                    if(k < 0 || (bp[k] < p - clump_bp))
                    {
                        backward = false;
                        k = j;
                        continue;
                    }
                }
                else
                {
                    ++k;
                    if(k >= bp.size() || (bp[k] > p + clump_bp)) break;
                }
                p2 = bp[k];
                if(mpp.count(p2) == 0) continue;
                double r =
                    (G.col(idx[j]).array() * G.col(idx[k]).array() * (sds(idx[j]) * sds(idx[k]))).sum() * df;
                if(r * r >= clump_r2)
                {
                    clumped.push_back(p2);
                    mpp.erase(p2);
                }
            }
            // what we do with clumped SNPs. sort them by pval?
            ofs << lines[p] << "\t";
            if(clumped.empty())
            {
                ofs << "NONE";
            }
            else
            {
                Double1D opp;
                for(auto op : clumped) opp.push_back(pvals_per_chr[c].at(op));
                for(auto oi : sortidx(opp)) ofs << clumped[oi] << ",";
            }
            ofs << std::endl;
        }
        // end current chr
    }
    std::ofstream ofs(fileout + ".clump.txt");
    string fn;
    for(int c = 0; c < (int)idx_per_chr.size(); c++)
    {
        fn = fileout + ".clump.chr" + std::to_string(c + 1);
        std::ifstream ifs(fn);
        ofs << ifs.rdbuf();
        ifs.close();
        std::remove(fn.c_str());
    }
}

void run_ld_stuff(const Param & params)
{
    std::vector<std::string> chromosomes; // for ld stuff
    Int1D snp_pos; // for ld stuff
    Int1D chr_pos_end; // store the last SNP in snp_pos in each chromosom
    get_snp_pos_bim(params.filebim, snp_pos, chr_pos_end, chromosomes);
    MyMatrix G; // get G from .residuals file
    MyVector F; // allele frequency
    if(params.clump.empty())
        calc_ld_metrics(params.fileout, params.filebim, G, F, snp_pos, chr_pos_end, params.ld_bp,
                        params.ld_r2, params.verbose);
    else
        calc_ld_clump(params.fileout, params.clump, params.assoc_colnames, params.clump_bp, params.clump_r2,
                      params.clump_p1, params.clump_p2, G, snp_pos, chr_pos_end, chromosomes);
}

void run_ld_stuff(const Param & params, Data * data)
{
    std::vector<std::string> chromosomes; // for ld stuff
    Int1D snp_pos; // for ld stuff
    Int1D chr_pos_end; // store the last SNP in snp_pos in each chromosom
    get_snp_pos_bim(params.filebim, snp_pos, chr_pos_end, chromosomes);
    MyVector F; // allele frequency
    if(params.clump.empty())
        calc_ld_metrics(params.fileout, params.filebim, data->G, F, snp_pos, chr_pos_end, params.ld_bp,
                        params.ld_r2, params.verbose);
    else
    {
        std::string sep{","};
        const auto assocfiles = split_string(params.clump, sep);
        for(size_t i = 0; i < assocfiles.size(); i++)
            calc_ld_clump(params.fileout + ".pheno" + to_string(i), assocfiles[i], params.assoc_colnames,
                          params.clump_bp, params.clump_r2, params.clump_p1, params.clump_p2, data->G,
                          snp_pos, chr_pos_end, chromosomes);
    }
}
