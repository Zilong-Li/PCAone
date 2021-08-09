#include "FileBgen.hpp"


void FileBgen::get_matrix_dimensions()
{
    // open file and print file info
    bgen.open( params.bgen );
    params.verbose && bgen.summarise( cout );
    // get nsamples and nsnps;
    nsamples = bgen.number_of_samples();
    nsnps = bgen.number_of_variants();
    cout << timestamp() << "N samples is " << nsamples << ". M snps is " << nsnps << endl;
}

void FileBgen::read_all_and_centering()
{
    uint c, i, j, m = 0;
    uint position;
    std::string chromosome, rsid;
    std::vector< std::string > alleles ;
    std::vector< std::vector< double > > probs ;
    cout << timestamp() << "begin to parse the bgen file.\n";
    if (!params.pcangsd)
    {
        F = VectorXf::Zero(nsnps);
        G = MatrixXf(nsamples, nsnps);
        if (params.maxiter > 0) C.resize(nsnps * nsamples);
        while(  bgen.read_variant( &chromosome, &position, &rsid, &alleles ) )
        {
            bgen.read_probs( &probs );
            assert( probs[0].size() == 3 ); // only support haploidy
            c = 0;
            #pragma omp parallel for private(i, j) reduction(+ : c)
            for (i = 0; i < probs.size(); ++i)
            {
                // get GT or missing;
                for (j = 0; j < 3; ++j)
                {
                    if (probs[i][j] > GENOTYPE_THRESHOLD) break;
                }
                if (j < 3)
                {
                    // j = 0, 1, 2
                    G(i, m) = BGEN2GENO[j];
                    if (params.maxiter > 0) C[m * nsamples + i] = 0;
                    F(m) += G(i, m);
                    c += 1;
                } else {
                    // no allele with GP > GENOTYPE_THRESHOLD, so set it as missing
                    G(i, m) = BGEN2GENO[j]; // j = 3;
                }
            }

            if (c == 0)
            {
                F(m) = 0.0;
            } else {
                F(m) /= c;
            }
            if (F(m) == 0)
            {
                cerr << "Warning: the allele frequency should not be 0. should do filtering first.\n";
                exit(EXIT_FAILURE);
            }
            // do centering and initialing
            for(i=0; i<nsamples; ++i)
            {
                if (G(i, m) == BGEN_MISSING_VALUE) {
                    G(i, m) = 0.0;
                } else {
                    G(i, m) -= F(m);
                }
            }
            m++;
        }
        assert( m == nsnps );
    } else {
        // cerr << "read GP instead of GT for pcangsd\n";
        // read all GP data into P;
        while(  bgen.read_variant( &chromosome, &position, &rsid, &alleles ) )
        {
            bgen.read_probs( &probs );
            assert( probs[0].size() == 3 ); // only support haploidy
            #pragma omp parallel for private(i, j)
            for (i = 0; i < probs.size(); ++i)
            {
                for (j = 0; j < 3; ++j)
                {
                    P(m, i * 3 + j) = probs[i][j];
                }
            }
            m++;
        }
        assert( m == nsnps );
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
}