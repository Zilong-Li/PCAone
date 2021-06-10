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
    cerr << timestamp() << "begin to permute plink data.\n";
    uint nsnps = count_lines(fin + ".bim");
    uint nsamples = count_lines(fin + ".fam");
    uint bed_bytes_per_snp = (nsamples+3)>>2;
    vector<uint> seqs;
    seqs.reserve(nsnps);
    for(uint i = 0; i < nsnps; i++) { seqs.push_back(i); }
    // std::random_device r;
    // auto rng = std::default_random_engine {r()};
    auto rng = std::default_random_engine {};
    std::shuffle(std::begin(seqs), std::end(seqs), rng);

    string fout = fin + ".perm";
    std::ifstream in(fin + ".bed", std::ios::in | std::ios::binary);
    std::ofstream out(fout + ".bed", std::ios::out | std::ios::binary);
    if (!in.is_open())
    {
        cerr << "ERROR: Cannot open bed file.\n";
        exit(EXIT_FAILURE);
    }
    uchar header[3];
    in.read(reinterpret_cast<char *> (&header[0]), 3);
    if ( (header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01) ) {
        cerr << "ERROR: Incorrect magic number in bed file.\n";
        exit(EXIT_FAILURE);
    }
    out.write(reinterpret_cast<char *> (&header[0]), 3);
    vector<uchar> inbed;
    inbed.resize(bed_bytes_per_snp);
    std::ifstream in_bim(fin + ".bim", std::ios::in);
    std::ofstream out_bim(fout + ".bim", std::ios::out);
    vector<string> bims(std::istream_iterator<Line>{in_bim},
                        std::istream_iterator<Line>{});
    uint idx;
    for(uint i = 0; i < nsnps; i++)
    {
        in.read(reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_snp);
        // get the position where to write
        idx = 3 + seqs[i] * bed_bytes_per_snp;
        out.seekp(idx);
        out.write(reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_snp);
        out_bim << bims[seqs[i]] + "\n";
    }
    in.close(); out.close();
    std::ifstream in_fam(fin + ".fam");
    std::ofstream out_fam(fout + ".fam");
    out_fam << in_fam.rdbuf();
    // set bed_prefix to the new;
    fin = fout;
}

void flip_UV(MatrixXf& U, MatrixXf& V, bool ubase)
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
                    cerr << "Error: the dimention of U and V have different k ranks.\n";
                    exit(EXIT_FAILURE);
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
                cerr << "Error: the dimention of U and V have different k ranks.\n";
                exit(EXIT_FAILURE);
            }
        }
    }
}

void flip_Y(const MatrixXf& X, MatrixXf& Y)
{
    VectorXf a = (X.array().sign() - Y.array().sign()).abs().colwise().sum().matrix();
    for (Eigen::Index i = 0; i < X.cols(); ++i)
    {
        // if signs of half of values are flipped then correct signs.
        if (a(i) > X.rows())
        {
            Y.col(i) *= -1;
        }
    }
}

double rmse(const MatrixXf& X, const MatrixXf& Y)
{
    MatrixXf Z = Y;
    flip_Y(X, Z);
    double diff = sqrt((X - Z).array().square().sum() / (X.cols() * X.rows()));
    return diff;
}