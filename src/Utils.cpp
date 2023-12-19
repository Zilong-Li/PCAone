#include "Utils.hpp"

#include <stdexcept>
#include <unordered_map>
#include <vector>

using namespace std;

std::string get_machine()
{
    struct utsname unameData;
    if(uname(&unameData) != 0)
    {
        perror("uname");
        exit(EXIT_FAILURE);
    }
    std::string machine{unameData.machine};
    std::string node{unameData.nodename};
    std::string release{unameData.release};
    std::string version{unameData.version};
    std::string sysname{unameData.sysname};
    return "Machine name: " + machine + "\nNode name: " + node + "\nOperating system release: " + release
           + "\nOperating system version: " + version + "\nOperating system name: " + sysname + "\n";
}

void fcloseOrDie(FILE * file)
{
    if(!fclose(file))
    {
        return;
    };
    /* error */
    perror("fclose error");
    exit(1);
}

FILE * fopenOrDie(const char * filename, const char * instruction)
{
    FILE * const inFile = fopen(filename, instruction);
    if(inFile) return inFile;
    /* error */
    perror(filename);
    exit(1);
}

size_t freadOrDie(void * buffer, size_t sizeToRead, FILE * file)
{
    size_t const readSize = fread(buffer, 1, sizeToRead, file);
    if(readSize == sizeToRead) return readSize; /* good */
    if(feof(file)) return readSize; /* good, reached end of file */
    /* error */
    perror("fread");
    exit(1);
}

size_t count_lines(const std::string & fpath)
{
    std::ifstream in(fpath);
    if(!in.is_open()) throw invalid_argument("can not open " + fpath);
    size_t count = 0;
    std::string line;
    while(getline(in, line))
    {
        count++;
    }
    return count;
}

std::string timestamp()
{
    auto t1 = std::chrono::system_clock::now();
    std::time_t tc = std::chrono::system_clock::to_time_t(t1);
    std::string str(std::ctime(&tc));
    str.pop_back(); // str[str.size() - 1] = '.';
    str = std::string("[") + str + std::string("] ");
    return str;
}

// Sign correction to ensure deterministic output from SVD.
// see https://www.kite.com/python/docs/sklearn.utils.extmath.svd_flip
void flip_UV(MyMatrix & U, MyMatrix & V, bool ubase)
{
    if(ubase)
    {
        Eigen::Index x, i;
        for(i = 0; i < U.cols(); ++i)
        {
            U.col(i).cwiseAbs().maxCoeff(&x);
            if(U(x, i) < 0)
            {
                U.col(i) *= -1;
                if(V.cols() == U.cols())
                {
                    V.col(i) *= -1;
                }
                else if(V.rows() == U.cols())
                {
                    V.row(i) *= -1;
                }
                else
                {
                    throw std::runtime_error("the dimention of U and V have different k ranks.\n");
                }
            }
        }
    }
    else
    {
        Eigen::Index x, i;
        for(i = 0; i < V.cols(); ++i)
        {
            if(V.cols() == U.cols())
            {
                V.col(i).cwiseAbs().maxCoeff(&x);
                if(V(x, i) < 0)
                {
                    U.col(i) *= -1;
                    V.col(i) *= -1;
                }
            }
            else if(V.rows() == U.cols())
            {
                V.row(i).cwiseAbs().maxCoeff(&x);
                if(V(i, x) < 0)
                {
                    U.col(i) *= -1;
                    V.row(i) *= -1;
                }
            }
            else
            {
                throw std::runtime_error("the dimention of U and V have different k ranks.\n");
            }
        }
    }
}

void flip_Omg(MyMatrix & Omg2, MyMatrix & Omg)
{
    for(Eigen::Index i = 0; i < Omg.cols(); ++i)
    {
        // if signs of half of values are flipped then correct signs.
        if((Omg2.col(i) - Omg.col(i)).array().abs().sum()
           > 2 * (Omg2.col(i) + Omg.col(i)).array().abs().sum())
        {
            Omg.col(i) *= -1;
        }
    }
    Omg2 = Omg;
}

void flip_Y(const MyMatrix & X, MyMatrix & Y)
{
    for(Eigen::Index i = 0; i < X.cols(); ++i)
    {
        // if signs of half of values are flipped then correct signs.
        if((X.col(i) - Y.col(i)).array().abs().sum() > 2 * (X.col(i) + Y.col(i)).array().abs().sum())
        {
            Y.col(i) *= -1;
        }
    }
}

double rmse(const MyMatrix & X, const MyMatrix & Y)
{
    MyMatrix Z = Y;
    flip_Y(X, Z);
    return sqrt((X - Z).array().square().sum() / (X.cols() * X.rows()));
}

// Y is the truth matrix, X is the test matrix
Eigen::VectorXd minSSE(const MyMatrix & X, const MyMatrix & Y)
{
    Eigen::Index w1, w2;
    Eigen::VectorXd res(X.cols());
    for(Eigen::Index i = 0; i < X.cols(); ++i)
    {
        // test against the original matrix to find the index with mincoeff
        ((-Y).colwise() + X.col(i)).array().square().colwise().sum().minCoeff(&w1);
        // test against the flipped matrix with the opposite sign
        (Y.colwise() + X.col(i)).array().square().colwise().sum().minCoeff(&w2);
        // get the minSSE value for X.col(i) against -Y.col(w1)
        auto val1 = (-Y.col(w1) + X.col(i)).array().square().sum();
        // get the minSSE value for X.col(i) against Y.col(w2)
        auto val2 = (Y.col(w2) + X.col(i)).array().square().sum();
        if(w1 != w2 && val1 > val2)
            res[i] = val2;
        else
            res[i] = val1;
    }
    return res;
}

double mev(const MyMatrix & X, const MyMatrix & Y)
{
    double res = 0;
    for(Eigen::Index i = 0; i < X.cols(); ++i)
    {
        res += (X.transpose() * Y.col(i)).norm();
    }
    return res / X.cols();
}

void mev_rmse_byk(const MyMatrix & X, const MyMatrix & Y, MyVector & Vm, MyVector & Vr)
{
    for(Eigen::Index i = 0; i < X.cols(); ++i)
    {
        Vm(i) = 1 - mev(X.leftCols(i + 1), Y.leftCols(i + 1));
        Vr(i) = rmse(X.leftCols(i + 1), Y.leftCols(i + 1));
    }
}

double get_median(std::vector<double> v)
{
    size_t n = v.size();
    if(n == 0)
    {
        return 0;
    }
    else
    {
        std::sort(v.begin(), v.end());
        if(n % 2 == 0)
        {
            return (v[n / 2 - 1] + v[n / 2]) / 2.0;
        }
        else
        {
            return v[n / 2];
        }
    }
}

std::vector<std::string> split_string(const std::string & s, const std::string & separators)
{
    std::vector<std::string> ret;
    bool is_seperator[256] = {false};
    for(auto & ch : separators)
    {
        is_seperator[(unsigned int)ch] = true;
    }
    int begin = 0;
    for(int i = 0; i <= (int)s.size(); i++)
    {
        if(is_seperator[(uint8_t)s[i]] || i == (int)s.size())
        {
            ret.push_back(std::string(s.begin() + begin, s.begin() + i));
            begin = i + 1;
        }
    }
    return ret;
}

// could use a new struct
void get_snp_pos_bim(const std::string & filebim,
                     Int1D & pos,
                     Int1D & chr_pos_end,
                     std::vector<std::string> & chrs)
{
    std::ifstream fin(filebim);
    if(!fin.is_open()) throw invalid_argument("can not open " + filebim);
    std::string line, chr_cur, chr_prev, sep{" \t"};
    int i = 0;
    while(getline(fin, line))
    {
        auto tokens = split_string(line, sep);
        chr_cur = tokens[0];
        if(chr_prev.empty()) chrs.push_back(chr_cur);
        if(!chr_prev.empty() && chr_prev != chr_cur)
        {
            chr_pos_end.push_back(i);
            chrs.push_back(chr_cur);
        }
        chr_prev = chr_cur;
        pos.push_back(std::stoi(tokens[3]));
        i++;
    }
    chr_pos_end.push_back(i); // add the last SNP
}

// given a list of snps, find its index per chr in the original pos
// assume chromosomes are continuous
Int2D get_target_snp_idx(const std::string & filebim,
                         const Int1D & pos,
                         const Int1D & chr_pos_end,
                         const std::vector<std::string> & chrs)
{
    Int1D t_pos, t_chr_pos_end, idx;
    std::vector<std::string> t_chrs;
    get_snp_pos_bim(filebim, t_pos, t_chr_pos_end, t_chrs);
    std::unordered_map<int, int> mpos;
    std::ifstream fin(filebim);
    if(!fin.is_open()) throw invalid_argument("can not open " + filebim);
    std::string line, chr_cur, chr_prev, sep{" \t"};
    int c, s, e, p, i;
    Int2D ret(t_chrs.size());
    for(int tc = 0; tc < (int)t_chrs.size(); tc++)
    {
        for(c = 0; c < (int)chrs.size(); c++)
            if(chrs[c] == t_chrs[tc]) break;
        e = chr_pos_end[c];
        s = c > 0 ? chr_pos_end[c - 1] : 0;
        for(i = s; i < e; i++) mpos[pos[i]] = i;
        e = t_chr_pos_end[tc];
        s = tc > 0 ? t_chr_pos_end[tc - 1] : 0;
        for(i = s; i < e; i++)
        {
            p = t_pos[i];
            if(mpos.count(p)) idx.push_back(mpos[p]);
        }
        ret[tc] = idx;
        idx.clear();
        mpos.clear();
    }
    return ret;
}

MyVector calc_sds(const MyMatrix & X)
{
    // compute degree of freedom
    const int df = X.rows() - 1; // N-1
    return (X.array().square().colwise().sum() / df).sqrt();
}

void calc_ld_metrics(std::string fileout,
                     MyMatrix & G,
                     const MyVector & F,
                     const Int1D & snp_pos,
                     const Int1D & chr_pos_end,
                     int ld_window_bp,
                     double r2_tol,
                     bool verbose = false)
{
    cao << tick.date() << "start calculating ld  metrics" << std::endl;
    G.rowwise() -= G.colwise().mean(); // Centering
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
                int o = F(k) > F(i) ? i : k;
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
                   MyMatrix & G,
                   const MyVector & F,
                   const Int1D & snp_pos,
                   const Int1D & chr_pos_end,
                   const std::vector<std::string> & chrs)
{
    cao << tick.date() << "start calculating pairwise ld r2 given a list of SNPs " << std::endl;
    G.rowwise() -= G.colwise().mean(); // Centering
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
            for(int k = j ; k < m; k++)
            {
                double r =
                    (G.col(idx[j]).array() * G.col(idx[k]).array() * (sds(idx[j]) * sds(idx[k]))).sum() * df;
                r *= r;
                ofs.write((char *)&r, sizeof(r));
            }
        }
    }
}
