#include "FileBgen.hpp"


void FileBgen::get_matrix_dimensions()
{
    // open file and print file info
    bg = new Bgen(params.bgen, "", true);
    // get nsamples and nsnps;
    nsamples = bg->header.nsamples;
    nsnps = bg->header.nvariants;
    cout << timestamp() << "the layout of bgen file is " << bg->header.layout << ". N samples is " << nsamples << ". M snps is " << nsnps << endl;
}

void FileBgen::read_all_and_centering()
{
    uint c, i, j;
    Variant var;
    cout << timestamp() << "begin to parse the bgen file.\n";
    if (!params.pcangsd)
    {
        F = VectorXf::Zero(nsnps);
        G = MatrixXf(nsamples, nsnps);
        if (params.maxiter > 0) C.resize(nsnps * nsamples);
        float* dosages;
        for (j = 0; j < nsnps; j++) {
            try {
                var = bg->next_var();
                dosages = var.minor_allele_dosage();
                c = 0;
                // calculate allele frequency
// #pragma omp parallel for private(i) reduction(+:c) reduction(+:VectorXf)
                for (i = 0; i < nsamples; i++) {
                    if (std::isnan(dosages[i])) {
                        if (params.maxiter > 0) C[j * nsamples + i] = 1;
                    } else {
                        if (params.maxiter > 0) C[j * nsamples + i] = 0;
                        G(i, j) = dosages[i] / 2.0; // map to [0, 1];
                        F(j) += G(i, j);
                        c += 1;
                    }
                }
                if (c==0) throw std::runtime_error("Error: the allele frequency should not be 0. should do filtering first.");
                F(j) /= c;
                // do centering and initialing
                #pragma omp parallel for
                for (i = 0; i < nsamples; i++) {
                    if (std::isnan(dosages[i])) {
                        G(i, j) = 0;
                    } else {
                        G(i, j) -= F(j);
                    }
                }
            } catch (const std::out_of_range & e) {
                break;
            }
        }
    } else {
        // read all GP data into P;
        float* probs1d;
        for (j = 0; j < nsnps; j++) {
            try {
                var = bg->next_var();
                probs1d = var.probs_1d();
                #pragma omp parallel for
                for (i = 0; i < nsamples; i++) {
                    P(j, i * 3 + 0) = probs1d[i * 3 + 0];
                    P(j, i * 3 + 1) = probs1d[i * 3 + 1];
                    P(j, i * 3 + 2) = probs1d[i * 3 + 2];
                }
            } catch (const std::out_of_range & e) {
                break;
            }
        }
        assert( j == nsnps );
        cout << timestamp() << "begin to estimate allele frequencies using GP" << endl;
        VectorXf Ft(nsnps);
        F = VectorXf::Constant(nsnps, 0.25);
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
                cout << "EM (MAF) converged at iteration: " << it+1 << endl;
                break;
            } else if (it == (params.maxiter-1)) {
                cerr << "EM (MAF) did not converge.\n";
            }
        }
        // initial E which is G
        G = MatrixXf(nsamples, nsnps);
        #pragma omp parallel for
        for (j = 0; j < nsnps; j++) {
            double p0, p1, p2;
            for (i = 0; i < nsamples; i++) {
                p0 = P(j, 3 * i + 0) * (1.0 - F(j)) * (1.0 - F(j));
                p1 = P(j, 3 * i + 1) * 2 * F(j) * (1.0 - F(j));
                p2 = P(j, 3 * i + 2) * F(j) * F(j);
                G(i, j) = (p1 + 2 * p2)/(p0 + p1 + p2) - 2.0 * F(j);
            }
        }
    }
}