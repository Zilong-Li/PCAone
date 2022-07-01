#include "FileBeagle.hpp"

using namespace std;

// inspired by Angsd
int tgets(gzFile gz, char** buf, uint64* l)
{
    int rlen = 0;
neverUseGoto:
    char* tok = gzgets(gz, *buf + rlen, *l - rlen);
    if (!tok)
        return rlen;
    int tmp = tok ? strlen(tok) : 0;
    if (tok[tmp - 1] != '\n')
    {
        rlen += tmp;
        *l *= 2;
        *buf = (char*)realloc(*buf, *l);
        goto neverUseGoto;
    }
    rlen += tmp;
    return rlen;
}


// read all data and estimate F
void FileBeagle::read_all()
{
    // read all GP data into P. parsing in parallel to speedup
#pragma omp parallel for
    for (uint64 j = 0; j < nsnps; j++) {
        size_t p;
        for(uint t = 0; t < 3; t++)
        {
            p = bgls[j].find("\t");
            bgls[j].erase(0, p + 1);
        }
        for (uint64 i = 0; i < nsamples; i++) {
            for (uint t = 0; t < 2; t++)
            {
                p = bgls[j].find("\t");
                if(p != std::string::npos) P(2 * i + t, j) = std::stod(bgls[j].substr(0, p));
                bgls[j].erase(0, p + 1);
            }
            p = bgls[j].find("\t");
            if(p != std::string::npos) bgls[j].erase(0, p + 1);
        }
    }

    bgls.clear();
    bgls.shrink_to_fit();
    llog << timestamp() << "begin to estimate allele frequencies" << endl;
    F = MyVector::Constant(nsnps, 0.25);
    { // out of scope: eigen object will be released;
        MyVector Ft = MyVector::Zero(nsnps);
        double diff;
        // run EM to estimate allele frequencies
        for (uint it = 0; it < params.maxiter; it++)
        {
#pragma omp parallel for
            for (uint j = 0; j < nsnps; j++)
            {
                Ft(j) = F(j);
                double p0, p1, p2, pt = 0.0;
                for (uint i = 0; i < nsamples; i++)
                {
                    p0 = P(2 * i + 0, j) * (1.0 - F(j)) * (1.0 - F(j));
                    p1 = P(2 * i + 1, j) * 2 * F(j) * (1.0 - F(j));
                    p2 = (1 - P(2 * i + 0, j) - P(2 * i + 1, j)) * F(j) * F(j);
                    pt += (p1 + 2 * p2) / (2 * (p0 + p1 + p2));
                }
                F(j) = pt / (double)nsamples;
            }
            // calculate differences between iterations
            diff = sqrt((F - Ft).array().square().sum() / nsnps);
            // Check for convergence
            if (diff < params.tolmaf)
            {
                llog << timestamp() << "EM (MAF) converged at iteration: " << it + 1 << endl;
                break;
            }
            else if (it == (params.maxiter - 1))
            {
                llog << timestamp() << "EM (MAF) did not converge.\n";
            }
        }
    }
    // filter snps and resize G;
    filterSNPs_resizeF();
    // resize P, only keep columns matching the indecis in idx;
    // P = P(Eigen::all, idx).eval(); // aliasing issue!!!
    G = MyMatrix::Zero(nsamples, nsnps); // initial E which is G
#pragma omp parallel for
    for (uint j = 0; j < nsnps; j++)
    {
        double p0, p1, p2;
        for (uint i = 0; i < nsamples; i++)
        {
            p0 = P(2 * i + 0, keepSNPs[j]) * (1.0 - F(j)) * (1.0 - F(j));
            p1 = P(2 * i + 1, keepSNPs[j]) * 2 * F(j) * (1.0 - F(j));
            p2 = (1 - P(2 * i + 0, keepSNPs[j]) - P(2 * i + 1, keepSNPs[j])) * F(j) * F(j);
            G(i, j) = (p1 + 2 * p2) / (p0 + p1 + p2) - 2.0 * F(j);
        }
    }
}
