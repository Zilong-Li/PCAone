#include "LD.hpp"

using namespace std;

MyVector calc_sds(const MyMatrix & X)
{
    // compute degree of freedom
    const int df = X.rows() - 1; // N-1
    return (X.array().square().colwise().sum() / df).sqrt();
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
    Int2D idx_per_chr = get_target_snp_idx(filebim, snp_pos, chr_pos_end, chrs);
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
    auto colidx = valid_assoc_file(fileassoc); // 0: chr, 1: pos, 2: pvalue
    Int2D idx_per_chr = get_target_snp_idx(fileassoc, snp_pos, chr_pos_end, chrs, true, colidx);
    // sort by pvalues and get new idx
    MyVector sds = 1.0 / calc_sds(G).array();
    const double df = 1.0 / (G.rows() - 1); // N-1
#pragma omp parallel for
    for(int i = 0; i < (int)idx_per_chr.size(); i++)
    {
        auto idx = idx_per_chr[i];
        int m = idx.size();
        ArrayXb keep = ArrayXb::Constant(m, true);
        for(int j = 0; j < m; j++)
        {
            if(!keep(j)) continue;
            for(int k = j + 1; k < m; k++)
            {
                double r =
                    (G.col(idx[j]).array() * G.col(idx[k]).array() * (sds(idx[j]) * sds(idx[k]))).sum() * df;
                if(r * r > clump_r2) keep(k) = false;
            }
        }
        // now output current chr
    }
}
