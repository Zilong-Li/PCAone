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

// structured permutation with cached buffer
void permute_plink2(string& fin, uint gb)
{
    uint   nbands = 64;
    uint64 nsnps = count_lines(fin + ".bim");
    uint64 nsamples = count_lines(fin + ".fam");
    uint64 bed_bytes_per_snp = (nsamples+3)>>2;

    // calculate the readin number of snps of certain big buffer like 2GB.
    // must be a multiple of nbands.
    uint64 twoGB = (uint64)1073741824 * gb;
    uint64 twoGB_snps = (uint64)floor((double) twoGB / bed_bytes_per_snp);
    if (twoGB_snps > nsnps) twoGB_snps = nsnps;
    uint64 bufsize = (uint64)floor((double) twoGB_snps / nbands);
    twoGB_snps = bufsize * nbands;     // initially twoGB_snps is a multiple of nbands
    assert(nsnps >= twoGB_snps);
    uint64 nblocks = (uint64)ceil((double) nsnps / twoGB_snps);
    uint modr2 = nsnps % twoGB_snps;
    uint64 bed_bytes_per_block = bed_bytes_per_snp * twoGB_snps;
    vector<uchar> inbed;                      // keep the input buffer
    inbed.resize(bed_bytes_per_block);
    vector<uchar> outbed;                     // keep the output buffer
    uint64 out_bytes_per_block = bed_bytes_per_snp * bufsize;
    outbed.resize(out_bytes_per_block);

    // get index of first snp of each band
    vector<uint64> bandidx;
    bandidx.resize(nbands);
    uint modr = nsnps % nbands;
    uint64 bandsize = (uint64)ceil((double) nsnps / nbands);
    if (modr == 0) {
        for (uint i = 0; i < nbands; ++i) {
            bandidx[i] = i * bandsize;
        }
    } else {
        for (uint i = 0; i < nbands; ++i) {
            if (i < modr) {
                bandidx[i] = i * bandsize;
            } else {
                bandidx[i] = modr * bandsize + (bandsize - 1) * (i - modr);
            }
        }
    }

    ios_base::sync_with_stdio(false);
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
    std::ifstream in_bim(fin + ".bim", std::ios::in);
    std::ofstream out_bim(fout + ".bim", std::ios::out);
    vector<string> bims(std::istream_iterator<Line>{in_bim},
                        std::istream_iterator<Line>{});
    vector<string> bims2;
    bims2.resize(nsnps);
    uint64 b, i, j, twoGB_snps2, idx, bufidx=bufsize;
    for(i = 0; i < nblocks; i++) {
        if (i == nblocks - 1 && modr2 != 0) {
            twoGB_snps2 = nsnps - (nblocks - 1) * twoGB_snps;
            bed_bytes_per_block = bed_bytes_per_snp * twoGB_snps2;
            inbed.resize(bed_bytes_per_block);
            // in last block, twoGB_snps is not neccessary a multiple of nbands and smaller than the previous
            bufsize = (uint64)ceil((double) twoGB_snps2 / nbands);
            modr2 = twoGB_snps2 % nbands;
            out_bytes_per_block = bed_bytes_per_snp * bufsize;
            outbed.resize(out_bytes_per_block);
        }
        in.read(reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_block);
        for (b = 0; b < nbands; b++) {
            idx = 3 + (i * bufidx + bandidx[b]) * bed_bytes_per_snp;
            for (j = 0; j < bufsize - 1; j++) {
                std::copy(inbed.begin()+(j * nbands + b) * bed_bytes_per_snp, inbed.begin()+(j * nbands + b + 1) * bed_bytes_per_snp, outbed.begin() + j * bed_bytes_per_snp);
                // cout << i * twoGB_snps + j * nbands + b << endl;
                bims2[i * bufidx + bandidx[b] + j] = bims[i * twoGB_snps + j * nbands + b];
            }
            if (i != nblocks - 1 || (i == nblocks - 1 && b < modr2) || modr2 == 0) {
                std::copy(inbed.begin()+(j * nbands + b) * bed_bytes_per_snp, inbed.begin()+(j * nbands + b + 1) * bed_bytes_per_snp, outbed.begin() + j * bed_bytes_per_snp);
                bims2[i * bufidx + bandidx[b] + j] = bims[i * twoGB_snps + j * nbands + b];
            } else {
                out_bytes_per_block = bed_bytes_per_snp * (bufsize - 1);

            }
            out.seekp(idx, std::ios_base::beg);
            out.write(reinterpret_cast<char *> (&outbed[0]), out_bytes_per_block);
        }

    }
    for(auto b : bims2) {out_bim << b << "\n";}
    in.close(); out.close(); out_bim.close();

    std::ifstream in_fam(fin + ".fam");
    std::ofstream out_fam(fout + ".fam");
    out_fam << in_fam.rdbuf();
    fin = fout;
}

void permute_plink(string& fin, uint blocksize)
{
    cout << timestamp() << "begin to permute plink data.\n";
    uint64 nsnps = count_lines(fin + ".bim");
    uint64 nsamples = count_lines(fin + ".fam");
    uint64 bed_bytes_per_snp = (nsamples+3)>>2;
    uint64 bed_bytes_per_block = bed_bytes_per_snp * blocksize;
    uint64 nblocks = (unsigned int)ceil((double)nsnps / blocksize);

    setlocale(LC_ALL,"C");
    ios_base::sync_with_stdio(false);
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
    inbed.resize(bed_bytes_per_block);
    std::ifstream in_bim(fin + ".bim", std::ios::in);
    std::ofstream out_bim(fout + ".bim", std::ios::out);
    vector<string> bims;
    bims.resize(nblocks);
    vector<uint64> seqs;
    seqs.reserve(nblocks);
    string line;
    uint64 i, j;
    for(i = 0; i < nblocks; i++) {
        seqs.push_back(i);
        j = 0;
        while(getline(in_bim, line)) {
            if (j < blocksize) {
                bims[i] += line + "\n";
                j++;
                if (j >= blocksize) break;
            }
        }
    }
    auto rng = std::default_random_engine {};
    std::shuffle(std::begin(seqs), std::end(seqs), rng);
    uint64 idx;
    uint64 len = 3 + bed_bytes_per_snp * nsnps;
    string sites = "";
    for(i = 0; i < nblocks; i++)
    {
        idx = 3 + (seqs[i] + 1) * bed_bytes_per_snp * blocksize;
        if (idx >= len -1) {
            bed_bytes_per_block = len - 3 - (nblocks - 1) * bed_bytes_per_snp * blocksize;
        } else {
            bed_bytes_per_block = blocksize * bed_bytes_per_snp;
        }
        idx = 3 + seqs[i] * bed_bytes_per_snp * blocksize;
        in.read(reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_block);
        out.seekp(idx, std::ios_base::beg);
        out.write(reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_block);
        sites += bims[seqs[i]];
    }
    out_bim << sites;
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

void flip_Omg(MatrixXd& Omg2, MatrixXd& Omg)
{
    for (Eigen::Index i = 0; i < Omg.cols(); ++i)
    {
        // if signs of half of values are flipped then correct signs.
        if ((Omg2.col(i) - Omg.col(i)).array().abs().sum() > 2 * (Omg2.col(i) + Omg.col(i)).array().abs().sum())
        {
            Omg.col(i) *= -1;
        }
    }
    Omg2 = Omg;
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
