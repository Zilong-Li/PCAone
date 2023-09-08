#include "Utils.hpp"

#include <stdexcept>
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

// assert(X.rows() == Y.rows())
// X is 2D, Y is 1D
// X and Y are centered beforehand
MyVector cor_cross(const MyMatrix & X, const MyVector & Y)
{
    if(X.rows() != Y.rows())
        throw std::runtime_error("X and Y has incompatible dimension. the rows are not the same!\n");
    // compute degree of freedom
    const int df = X.rows() - 1; // N-1
    // the standard deviations of X and Y
    double inv_sds_Y = std::sqrt(Y.array().square().sum() / df);
    MyVector r(X.cols());
    for(auto i = 0; i < X.cols(); i++)
    {
        double inv_sds_X = std::sqrt(X.col(i).array().square().sum() / df);
        r(i) = (X.col(i).array() * Y.array() / (inv_sds_X * inv_sds_Y)).sum() / df;
    }
    return r;
}

// X is nsamples x nsnps;
// ws stores the starts of windows
// we stores the number of sites in a window
// return cor matrix of (nsnps - w + 1, w);
void cor_by_window(const std::string & fileout,
                   MyMatrix & X,
                   const std::vector<int> & ws,
                   const std::vector<int> & we,
                   double r2_tol = 0.0)
{
    X.rowwise() -= X.colwise().mean(); // Centering
    std::ofstream ofs_out(fileout + ".ld.prune.out");
    std::ofstream ofs_in(fileout + ".ld.prune.in");
    vector<ArrayXb> flags(ws.size());
#pragma omp parallel for
    for(int i = 0; i < ws.size(); i++)
    {
        MyVector r2 = cor_cross(X.middleCols(ws[i], we[i]), X.col(ws[i])).array().square();
        ArrayXb ikeep = ArrayXb::Constant(we[i], true);
        for(int j = 1; j < we[i]; j++)
            if(r2[j] > r2_tol) ikeep(j) = false;
        flags[i] = ikeep;
    }
    ArrayXb keep = ArrayXb::Constant(X.cols(), true);
    for(int i = 0; i < ws.size(); i++)
        for(int j = 1; j < we[i]; j++)
            if(!flags[i][j]) keep(ws[i] + j) = false;
    std::ifstream fin(fileout + ".kept.bim");
    if(!fin.is_open()) throw invalid_argument("can not open " + fileout + ".kept.bim");
    std::string line, chr_cur, chr_prev, sep{" \t"};
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

void calc_ld_metrics(const std::string & fileout,
                     MyMatrix & G,
                     const MyMatrix & U,
                     const MyVector & S,
                     const MyMatrix & V,
                     const std::vector<int> & snp_pos,
                     const std::vector<int> & chr_pos_end,
                     int ld_window_bp,
                     double r2_tol = 0.5)
{
    cao << tick.date() << "start calculating ld  metrics" << std::endl;
    G -= U * S.asDiagonal() * V.transpose(); // get residuals matrix
#if defined(DEBUG)
    std::ofstream ofs_res(fileout + ".residuals");
    ofs_res.write((char *)G.data(), G.size() * sizeof(double));
#endif
    std::ofstream ofs_win(fileout + ".ld.window");
    ofs_win << "#window\tchr\tpos_start\tpos_end\tnsites" << std::endl;
    std::vector<int> ws, we;
    int nsnp = snp_pos.size();
    int i{0}, j{0}, w{0}, c{0}, pos_end, nsites;
    for(i = 0; i < nsnp; i++)
    {
        pos_end = snp_pos[chr_pos_end[c]];
        for(j = i; j < chr_pos_end[c]; j++)
            if(snp_pos[j] > snp_pos[i] + ld_window_bp) break;
        nsites = j - i + 1;
        ws.push_back(i); // start of the window
        we.push_back(nsites); // the number of sites
        ofs_win << w++ << "\t" << c + 1 << "\t" << snp_pos[i] << "\t" << snp_pos[j] << "\t" << nsites
                << std::endl;
        if(c == chr_pos_end.size() - 1 && pos_end < snp_pos[i] + ld_window_bp) break;
        if(c < chr_pos_end.size() - 1 && pos_end < snp_pos[i] + ld_window_bp) i = chr_pos_end[c++];
    }
    cor_by_window(fileout, G, ws, we, r2_tol);
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

// resize pos beforehand
void get_snp_pos_bim(const std::string & filebim, std::vector<int> & pos, std::vector<int> & chr_pos_end)
{
    std::ifstream fin(filebim);
    if(!fin.is_open()) throw invalid_argument("can not open " + filebim);
    std::string line, chr_cur, chr_prev, sep{" \t"};
    int i = 0;
    while(getline(fin, line))
    {
        auto tokens = split_string(line, sep);
        chr_cur = tokens[0];
        if(!chr_prev.empty() && chr_prev != chr_cur) chr_pos_end.push_back(i);
        chr_prev = chr_cur;
        pos[i++] = std::stoi(tokens[3]);
    }
    chr_pos_end.push_back(i); // add the last SNP
}
