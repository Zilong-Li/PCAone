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

void calc_ld_metrics(std::string fileout,
                     const MyMatrix & G,
                     const MyVector & F,
                     const Int1D & snp_pos,
                     const Int1D & chr_pos_end,
                     int ld_window_bp,
                     double r2_tol,
                     bool verbose = false)
{
    cao << tick.date() << "start calculating ld  metrics" << std::endl;
    // G.rowwise() -= G.colwise().mean(); // Centering
#if defined(DEBUG)
    std::ofstream ofs_res(fileout + ".residuals");
    ofs_res.write((char *)G.data(), G.size() * sizeof(double));
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
                int o = MAF(F(k)) > MAF(F(i)) ? i : k;
                keep(o) = false;
            }
        }
    }
    std::ifstream fin(fileout + ".kept.bim");
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

void calc_ld_pairs(std::string fileout,
                   std::string filebim,
                   const MyMatrix & G,
                   const MyVector & F,
                   const Int1D & snp_pos,
                   const Int1D & chr_pos_end,
                   const std::vector<std::string> & chrs)
{
    cao << tick.date() << "start calculating pairwise ld r2 given a list of SNPs " << std::endl;
    // G.rowwise() -= G.colwise().mean(); // Centering
    MyVector sds = 1.0 / calc_sds(G).array();
    const double df = 1.0 / (G.rows() - 1); // N-1
    Int2D idx_per_chr, bp_per_chr;
    std::tie(idx_per_chr, bp_per_chr) = get_target_snp_idx(filebim, snp_pos, chr_pos_end, chrs);
#pragma omp parallel for
    for(int i = 0; i < (int)idx_per_chr.size(); i++)
    {
        auto idx = idx_per_chr[i];
        std::ofstream ofs(fileout + ".ld.chr." + std::to_string(i + 1), std::ios::binary);
        int m = idx.size();
        ofs.write((char *)&m, sizeof(m));
        // calc pairwise r2 for G[,idx]
        for(int j = 0; j < m; j++)
        {
            // output diagnal, k = j
            for(int k = j; k < m; k++)
            {
                double r =
                    (G.col(idx[j]).array() * G.col(idx[k]).array() * (sds(idx[j]) * sds(idx[k]))).sum() * df;
                r *= r;
                ofs.write((char *)&r, sizeof(r));
            }
        }
    }
}

Int1D valid_assoc_file(const std::string & fileassoc)
{
    std::ifstream fin(fileassoc);
    if(!fin.is_open()) throw invalid_argument("can not open " + fileassoc);
    std::string line, sep{"\t"};
    getline(fin, line);
    const auto fields = split_string(line, sep);
    Int1D idx(5, -1);
    int j = 0;
    for(auto col : fields)
    {
        if(col == "CHR")
        {
            idx[0] = j;
        }
        else if(col == "BP")
        {
            idx[1] = j;
        }
        else if(col == "P")
        {
            idx[2] = j;
        }
        else if(col == "A1")
        {
            idx[3] = j;
        }
        else if(col == "A2")
        {
            idx[4] = j;
        }
        j++;
    }
    j = 0;
    for(auto i : idx)
    {
        if(i < 0) cao.error("the " + to_string(j) + "-th field did not exists");
        j++;
    }
    return idx;
}

std::vector<UMapIntDouble> map_index_snps(const std::string & fileassoc,
                                          const Int1D & colidx,
                                          double clump_p1)
{
    std::ifstream fin(fileassoc);
    if(!fin.is_open()) throw invalid_argument("can not open " + fileassoc);
    std::string line, chr_cur, chr_prev, sep{" \t"};
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
        if(pval <= clump_p1) m.insert({bp, pval});
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
    std::string line, chr_cur, chr_prev, sep{" \t"};
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
                   int clump_bp,
                   double clump_r2,
                   double clump_p1,
                   double clump_p2,
                   const MyMatrix & G,
                   const MyVector & F,
                   const Int1D & snp_pos,
                   const Int1D & chr_pos_end,
                   const std::vector<std::string> & chrs)
{
    cao << tick.date() << "start do LD-based clumping -> p1=" << clump_p1 << ", p2=" << clump_p2
        << ", r2=" << clump_r2 << ", bp=" << clump_bp << std::endl;
    auto colidx = valid_assoc_file(fileassoc); // 0: chr, 1: pos, 2: pvalue
    Int2D idx_per_chr, bp_per_chr;
    std::tie(idx_per_chr, bp_per_chr) =
        get_target_snp_idx(fileassoc, snp_pos, chr_pos_end, chrs, true, colidx);
    const auto pvals_per_chr = map_index_snps(fileassoc, colidx, clump_p1);
    const auto line_per_chr = map_assoc_file(fileassoc, colidx);
    // sort by pvalues and get new idx
    MyVector sds = 1.0 / calc_sds(G).array();
    const double df = 1.0 / (G.rows() - 1); // N-1
#pragma omp parallel for
    for(int c = 0; c < (int)idx_per_chr.size(); c++)
    {
        const auto idx = idx_per_chr[c];
        const auto bp = bp_per_chr[c];
        auto pvals = pvals_per_chr[c]; // key: pos, val: pval
        auto mbp = vector2map(bp);
        auto lines = line_per_chr[c];
        std::ofstream ofs(fileout + ".clump.chr" + std::to_string(c + 1));
        ofs << line_per_chr[0].at(-1) << "\tSP2\n";
        // greedy clumping algorithm
        for(auto i : sortidx(pvals))
        { // snps sorted by p value
            int p = bp[i];
            if(pvals.count(p) == 0)
                continue; // if snps with pval < clump_p1 are already clumped with previous snps
            int p2, k = mbp[p], j = mbp[p]; // j:cur, k:forward or backward
            Int1D clumped;
            while(--k && k >= 0)
            {
                p2 = bp[k];
                if(p2 < p - clump_bp) break;
                double r =
                    (G.col(idx[j]).array() * G.col(idx[k]).array() * (sds(idx[j]) * sds(idx[k]))).sum() * df;
                if(r * r >= clump_r2 && pvals[p2] <= clump_p2)
                {
                    clumped.push_back(p2);
                    if(pvals[p2]) pvals.erase(p2);
                }
            }
            k = j;
            while(++k && k < bp.size())
            {
                p2 = bp[k];
                if(p2 > p + clump_bp) break;
                double r =
                    (G.col(idx[j]).array() * G.col(idx[k]).array() * (sds(idx[j]) * sds(idx[k]))).sum() * df;
                if(r * r >= clump_r2 && pvals[p2] <= clump_p2)
                {
                    clumped.push_back(p2);
                    if(pvals[p2]) pvals.erase(p2);
                }
            }
            // what we do with clumped SNPs. output them!
            ofs << lines[p] << "\t"; // or copy line from the input file
            for(auto o : clumped) ofs << o << ",";
            ofs << std::endl;
        }
        // end current chr
    }
}
