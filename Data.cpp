#include "Data.hpp"
#include "SVD.hpp"

// #include <Spectra/SymEigsSolver.h>
#include <Spectra/contrib/PartialSVDSolver.h>
using namespace Spectra;

void Data::run()
{
    if (params.intype == "bfile") {
        if (params.batch) {
            read_bed_whole_genomat();
            init_whole_E();
            PartialSVDSolver< double, MatrixXd > svds(G, params.k, params.k * 2 + 1);
            cerr << timestamp() << "begin to do svds\n";
            uint nconv = svds.compute(params.maxiter, params.tol);
            if (nconv != params.k) {
                cerr << "warning: something wrong\n";
            }
            VectorXd svals = svds.singular_values();
            MatrixXd U = svds.matrix_U(params.k);
            MatrixXd V = svds.matrix_V(params.k);
            cerr << timestamp() << "Singular values found:\n" << svals << endl;
            // run_whole_EM();
        } else {
            cout << "blockwise\n";
        }
    } else if (params.intype == "pfile") {
        cout << "begin to read pfiles\n";
    } else if (params.intype == "bgen") {
        cout << "begin to read bgen\n";
    } else {
        cerr << "ERROR: You must use either --bfile, --pfile or --bgen.\n";
        exit(EXIT_SUCCESS);
    }
}

void Data::read_bed_whole_genomat()
{
    string fbed = params.bed_prefix + ".bed";
    string fbim = params.bed_prefix + ".bim";
    string ffam = params.bed_prefix + ".fam";
    nsamples = count_lines(ffam);
    nsnps = count_lines(fbim);
    cout << "nsamples is " << nsamples << "\nnsnps is " << nsnps << endl;
    G = MatrixXd(nsnps, nsamples);
    F = VectorXd(nsnps);

    std::ifstream bed_ifstream;
    bed_ifstream.open(fbed, std::ios::in | std::ios::binary);
    if (!bed_ifstream.is_open())
    {
        cerr << "ERROR: Cannot open bed file.\n";
        exit(EXIT_FAILURE);
    }
    uchar header[3];
    bed_ifstream.read( reinterpret_cast<char *> (&header[0]), 3);
    // check magic number of bed file
    if ( (header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01) ) {
        cerr << "ERROR: Incorrect magic number in bed file.\n";
        exit(EXIT_FAILURE);
    }
    // begin to decode the plink bed
    uint64 bed_bytes_per_snp = (nsamples+3)>>2; // get ceiling(nsamples/4)
    inbed.reserve(bed_bytes_per_snp);
    uchar buf;
    vector<double> geno_no_missing;
    cerr << "begin to decode the plink bed.\n";
    C.reserve(nsnps * nsamples);
    for(uint64 i=0; i<nsnps; ++i)
    {
        bed_ifstream.read( reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_snp);
        for (size_t b=0, j=0; b < inbed.size(); ++b)
        {
            buf = inbed[b];
            for (uint k=0; k<4; ++k, ++j)
            {
                if (j < nsamples)
                {
                    G(i, j) = MAP2GENO[buf & 3];
                    if (G(i,j) != MISS_INTEGER) {
                        geno_no_missing.push_back(G(i,j));
                        C.push_back(0); // 0 indicate G(i,j) don't need to be predicted.
                    } else {
                        C.push_back(1); // 1 indicate G(i,j) need to be predicted and updated.
                    }
                    buf = buf >> 2;  // shift packed data and throw away genotype just processed.
                } else {
                    // when j is out of range, we're in the padding data now
                    // as an extra sanity check, the remaining data should be all zero (that's how the encoding is supposed to work)
                    if (buf != 0)
                    {
                        cerr << "Row " << i + 1 << "padding is non-zero. Either the specified number of individuals is incorrect or the input file is corrupt!\n";
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
        F(i) = 0.0;
        for (auto g: geno_no_missing) F(i) += g;
        F(i) /= 2 * geno_no_missing.size();
        geno_no_missing.clear();
    }
    // read bed to MatrixXd done.
}

void Data::init_whole_E()
{
    // #pragma omp parallel for
    for(uint64 i=0; i<nsnps; ++i)
    {
        for(uint64 j=0; j<nsamples; ++j)
        {
            if (G(i,j) == MISS_INTEGER) {
                G(i,j) = 0;
            } else {
                G(i,j) -= 2 * F(i);
            }
        }
    }
}

void Data::run_whole_EM()
{
    // SvdOp op(G);
    // SymEigsSolver< double, LARGEST_ALGE, SvdOp > eigs(&op, params.k, params.k * 2 + 1);
    // eigs.init();
    // cerr << timestamp() << "begin to do eigs\n";
    // eigs.compute(params.maxiter, params.tol);
    // VectorXd evalues;
    // if(eigs.info() == Spectra::SUCCESSFUL)
    //     evalues = eigs.eigenvalues();

    // cout << timestamp() << "Eigenvalues found:\n" << evalues << endl;
    PartialSVDSolver< double, MatrixXd > svds(G, params.k, params.k * 2 + 1);
    cerr << timestamp() << "begin to do svds\n";
    uint nconv = svds.compute(params.maxiter, params.tol);
    if (nconv != params.k) {
        cerr << "warning: something wrong\n";
    }
    VectorXd svals = svds.singular_values();
    MatrixXd U = svds.matrix_U(params.k);
    MatrixXd V = svds.matrix_V(params.k);
    cerr << timestamp() << "Singular values found:\n" << svals << endl;
}

void Data::update_whole_E()
{

}