#include "Utils.hpp"

#include <stdexcept>

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
// X and Y are centered beforehand
MyMatrix cor_cross(const MyMatrix & X, const MyMatrix & Y)
{
    if(X.rows() != Y.rows())
        throw std::runtime_error("X and Y has incompatible dimension. the rows are not the same!\n");
    // compute degree of freedom
    const int df = X.rows() - 1; // N-1
    MyMatrix cor = X.transpose() * Y / df;
    // Compute 1 over the standard deviations of X and Y
    Eigen::VectorXd inv_sds_X = (X.colwise().norm() / sqrt(df)).array().inverse();
    Eigen::VectorXd inv_sds_Y = (Y.colwise().norm() / sqrt(df)).array().inverse();
    // Scale the covariance matrix
    cor = cor.cwiseProduct(inv_sds_X * inv_sds_Y.transpose());
    return cor;
}

// X is nsamples x nsnps;
// wi stores the SNP index that all SNPs in a window compare to
// ws stores the starts of windows
// we stores the ends of windows
// return cor matrix of (nsnps - w + 1, w);
void cor_by_window(MyMatrix & X,
                   const std::vector<int> & wi,
                   const std::vector<int> & ws,
                   const std::vector<int> & we)
{
    X.rowwise() -= X.colwise().mean(); // Centering
    for(int i = 0; i < wi.size(); i++)
    {
        auto r2 = cor_cross(X.middleCols(ws[i], we[i]), X.col(wi[i])).array().square();
    }
}

void calc_ld_metrics(const std::string & fileout,
                     MyMatrix & G,
                     const MyMatrix & U,
                     const MyVector & S,
                     const MyMatrix & V)
{
    G -= U * S * V.transpose(); // get residuals matrix
    std::ofstream ofs_res(fileout + ".residuals");
    ofs_res.write((char *)G.data(), G.size() * sizeof(double));
}
