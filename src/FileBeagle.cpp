#include "FileBeagle.hpp"

using namespace std;

// inspired by Angsd
int tgets(gzFile gz,char**buf,int *l)
{
    int rlen = 0;
neverUseGoto:
    char *tok = gzgets(gz,*buf+rlen,*l-rlen);
    if(!tok) return rlen;
    int tmp = tok?strlen(tok):0;
    if(tok[tmp-1]!='\n'){
        rlen += tmp;
        *l *= 2;
        *buf = (char*) realloc(*buf,*l);
        goto neverUseGoto;
    }
    rlen += tmp;
    return rlen;
}


// read all data and estimate F
void FileBeagle::read_all_and_centering()
{
    cout << timestamp() << "begin to read whole data" << endl;
    fp = gzopen(params.beagle.c_str(), "r");
    tgets(fp, &buffer, &bufsize);
    char *tok;
    uint i=0, j=0;
    // read all GP data into P
    while(tgets(fp, &buffer, &bufsize)) {
        if(buffer!=original){
            original=buffer;
        }
        tok = strtok_r(buffer,delims,&buffer);
        tok = strtok_r(NULL,delims,&buffer);
        tok = strtok_r(NULL,delims,&buffer);
        for (i = 0; i < nsamples * 3; i++) {
            tok = strtok_r(NULL,delims,&buffer);
            assert(tok!=NULL);
            P(j, i) = std::stod(std::string(tok));
        }
        buffer = original;
        j++;
    }
    gzclose(fp);
    assert(j==nsnps);

    cout << timestamp() << "begin to estimate allele frequencies" << endl;
    MyVector Ft(nsnps);
    F = MyVector::Constant(nsnps, 0.25);
    // run EM to estimate allele frequencies
    double diff;
    for (uint it = 0; it < params.maxiter; it++)
    {
        #pragma omp parallel for
        for (uint j = 0; j < nsnps; j++) {
            Ft(j) = F(j);
            double p0, p1, p2, pt = 0.0;
            for (uint i = 0; i < nsamples; i++) {
                p0 = P(j, 3 * i + 0) * (1.0 - F(j)) * (1.0 - F(j));
                p1 = P(j, 3 * i + 1) * 2 * F(j) * (1.0 - F(j));
                p2 = P(j, 3 * i + 2) * F(j) * F(j);
                pt += (p1 + 2 * p2) / (2 * (p0 + p1 + p2));
            }
            F(j) = pt / (double) nsamples;
        }
        // calculate differences between iterations
        diff = sqrt((F - Ft).array().square().sum() / nsnps);
        // Check for convergence
        if (diff < params.tolmaf) {
            cout << timestamp() << "EM (MAF) converged at iteration: " << it+1 << endl;
            break;
        } else if (it == (params.maxiter-1)) {
            cout << timestamp() << "EM (MAF) did not converge.\n";
        }
    }
    // initial E which is G
    G = MyMatrix(nsamples, nsnps);
    #pragma omp parallel for
    for (uint j = 0; j < nsnps; j++) {
        double p0, p1, p2;
        for (uint i = 0; i < nsamples; i++) {
            p0 = P(j, 3 * i + 0) * (1.0 - F(j)) * (1.0 - F(j));
            p1 = P(j, 3 * i + 1) * 2 * F(j) * (1.0 - F(j));
            p2 = P(j, 3 * i + 2) * F(j) * F(j);
            G(i, j) = (p1 + 2 * p2)/(p0 + p1 + p2) - 2.0 * F(j);
        }
    }
}
