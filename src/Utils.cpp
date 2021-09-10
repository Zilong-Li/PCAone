#include "Utils.hpp"

size_t count_lines(const string& fpath)
{
    std::ifstream in(fpath);
    size_t count = 0;
    string line;
    while (getline(in, line)) {
        count++;
    }
    return count;
}

string timestamp()
{
    time_t t = time(NULL);
    char *s = asctime(localtime(&t));
    s[strlen(s) - 1] = '\0';
    string str(s);
    str = string("[") + str + string("] ");
    return str;
}

void permute_plink(string& fin)
{
    cout << timestamp() << "begin to permute plink data.\n";
    uint64 nsnps = count_lines(fin + ".bim");
    uint64 nsamples = count_lines(fin + ".fam");
    uint64 bed_bytes_per_snp = (nsamples+3)>>2;
    vector<uint64> seqs;
    seqs.reserve(nsnps);
    for(uint64 i = 0; i < nsnps; i++) { seqs.push_back(i); }
    auto rng = std::default_random_engine {};
    std::shuffle(std::begin(seqs), std::end(seqs), rng);

    string fout = fin + ".perm";
    std::ifstream in(fin + ".bed", std::ios::binary);
    std::ofstream out(fout + ".bed", std::ios::binary);
    if (!in.is_open()) {
        throw std::invalid_argument("ERROR: Cannot open bed file.\n");
    }
    uchar header[3];
    in.read(reinterpret_cast<char *> (&header[0]), 3);
    if ( (header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01) ) {
        throw std::invalid_argument("ERROR: Incorrect magic number in bed file.\n");
    }
    out.write(reinterpret_cast<char *> (&header[0]), 3);
    vector<uchar> inbed;
    inbed.resize(bed_bytes_per_snp);
    std::ifstream in_bim(fin + ".bim", std::ios::in);
    std::ofstream out_bim(fout + ".bim", std::ios::out);
    vector<string> bims(std::istream_iterator<Line>{in_bim},
                        std::istream_iterator<Line>{});
    uint64 idx;
    for(uint64 i = 0; i < nsnps; i++)
    {
        idx = 3 + seqs[i] * bed_bytes_per_snp;
        in.seekg(idx, std::ios_base::beg);
        in.read(reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_snp);
        out.write(reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_snp);
        out_bim << bims[seqs[i]] + "\n";
    }
    in.close(); out.close(); out_bim.close();
    std::ifstream in_fam(fin + ".fam");
    std::ofstream out_fam(fout + ".fam");
    out_fam << in_fam.rdbuf();
    // set bed_prefix to the permuted one;
    fin = fout;
}

// Sign correction to ensure deterministic output from SVD.
// see https://www.kite.com/python/docs/sklearn.utils.extmath.svd_flip
void flip_UV(MatrixXd& U, MatrixXd& V, bool ubase)
{
    if (ubase)
    {
        Eigen::Index x, i;
        for (i = 0; i < U.cols(); ++i)
        {
            U.col(i).cwiseAbs().maxCoeff(&x);
            if (U(x, i) < 0)
            {
                U.col(i) *= -1;
                if (V.cols() == U.cols())
                {
                    V.col(i) *= -1;
                } else if (V.rows() == U.cols()) {
                    V.row(i) *= -1;
                } else {
                    throw std::runtime_error("Error: the dimention of U and V have different k ranks.\n");
                }
            }

        }
    } else {
        Eigen::Index x, i;
        for (i = 0; i < V.cols(); ++i)
        {
            if (V.cols() == U.cols())
            {
                V.col(i).cwiseAbs().maxCoeff(&x);
                if (V(x, i) < 0)
                {
                    U.col(i) *= -1;
                    V.col(i) *= -1;
                }
            } else if (V.rows() == U.cols()) {
                V.row(i).cwiseAbs().maxCoeff(&x);
                if (V(i, x) < 0)
                {
                    U.col(i) *= -1;
                    V.row(i) *= -1;
                }
            } else {
                throw std::runtime_error("Error: the dimention of U and V have different k ranks.\n");
            }
        }
    }
}

void flip_Y(const MatrixXd& X, MatrixXd& Y)
{
    for (Eigen::Index i = 0; i < X.cols(); ++i)
    {
        // if signs of half of values are flipped then correct signs.
        if ((X.col(i) - Y.col(i)).array().abs().sum() > 2 * (X.col(i) + Y.col(i)).array().abs().sum())
        {
            Y.col(i) *= -1;
        }
    }
}

double rmse(const MatrixXd& X, const MatrixXd& Y)
{
    MatrixXd Z = Y;
    flip_Y(X, Z);
    return sqrt( (X - Z).array().square().sum() / (X.cols() * X.rows()) );
}

double mev(const MatrixXd& X, const MatrixXd& Y)
{
    double res = 0;
    for (Eigen::Index i = 0; i < X.cols(); ++i)
    {
        res += (X.transpose() * Y.col(i)).norm();
    }
    return res / X.cols();
}

void mev_rmse_byk(const MatrixXd& X, const MatrixXd& Y, VectorXd& Vm, VectorXd& Vr)
{
    for (Eigen::Index i = 0; i < X.cols(); ++i)
    {
        Vm(i) = 1 - mev(X.leftCols(i+1), Y.leftCols(i+1));
        Vr(i) = rmse(X.leftCols(i+1), Y.leftCols(i+1));
    }
}