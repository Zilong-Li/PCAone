#include "Arnoldi.hpp"
#include <Spectra/SymEigsSolver.h>
#include <Spectra/contrib/PartialSVDSolver.h>

using namespace Spectra;

void ArnoldiOpData::perform_op(const float *x_in, float* y_out)
{
   Map<const VectorXf> x(x_in, n);
   Map<VectorXf> y(y_out, n);
   // open bed
   data->open_check_file();
   if (update) {
       data->read_snp_block_update(data->start[0], data->stop[0], U, S, VT, standardize);
   } else {
       data->read_snp_block_initial(data->start[0], data->stop[0], standardize);
   }
   y.noalias() = data->G * (data->G.transpose() * x);
   for (uint k=1; k < data->nblocks; ++k)
   {
       if (update) {
           data->read_snp_block_update(data->start[k], data->stop[k], U, S, VT, standardize);
       } else {
           data->read_snp_block_initial(data->start[k], data->stop[k], standardize);
       }
       //TODO: Kahan summation
       y.noalias() = y + data->G * (data->G.transpose() * x);
   }
   // close bed
   data->close_check_file();
   nops++;
   data->params.verbose && cerr << "Arnoldi Op=" << nops << ".\n";
}


void run_pca_with_arnoldi(Data* data, const Param& params)
{
    cout << timestamp() << "begin to do svds\n";
    VectorXf svals, evals;
    uint nconv, nu;
    double diff;
    if (params.batch)
    {
        MatrixXf U, V, V2;
        // SpMatrix sG = data->G.sparseView();
        // PartialSVDSolver< float, SpMatrix > svds(sG, params.k, params.ncv);
        PartialSVDSolver< float, MatrixXf > svds(data->G, params.k, params.ncv);
        if (params.maxiter == 0)
        {
            data->standardize_E();
        }
        nconv = svds.compute(params.imaxiter, params.itol);
        if (nconv != params.k) {
            params.verbose && cerr << "Warning: the nconv is not equal to k.\n";
            exit(EXIT_FAILURE);
        }
        U = svds.matrix_U(params.k);
        svals = svds.singular_values();
        if (params.maxiter == 0)
        {
            cout << timestamp() << "Final SVD done!\n";
            evals = svals.array().square() / data->nsnps;
            data->write_eigs_files(evals, U);
            return;
        }
        V = svds.matrix_V(params.k);
        flip_UV(U, V);
        cout << timestamp() << "begin to do EM iteration.\n";
        for (uint i = 1; i <= params.maxiter; ++i)
        {
            data->update_batch_E(U, svals, V.transpose());
            nconv = svds.compute(params.imaxiter, params.itol);
            svals = svds.singular_values();
            U = svds.matrix_U(params.k);
            V2 = svds.matrix_V(params.k);
            flip_UV(U, V2);
            diff = rmse(V2, V);
            params.verbose && cout << timestamp() << "Individual allele frequencies estimated (iter=" << i << "), RMSE=" << diff <<".\n";
            V = V2;
            if (diff < params.tol)
            {
                cout << timestamp() << "Come to convergence!\n";
                break;
            }
        }

        cout << timestamp() << "begin to standardize the matrix\n";
        if (params.pcangsd)
        {
            data->pcangsd_standardize_E(U, svals, V.transpose());
            MatrixXf C = data->G * data->G.transpose();
            C.array() /= (double) data->nsnps;
            C.diagonal() = data->Dc.array() / (double) data->nsnps;
            std::ofstream out_cov(params.outfile + ".cov");
            if (out_cov.is_open()) {
                out_cov << C << "\n";
            }
        } else {
            data->standardize_E();
        }
        // PartialSVDSolver< float, MatrixXf > svd(data->G, params.k, params.ncv);
        svds.compute(params.imaxiter, params.itol);
        svals = svds.singular_values();
        U = svds.matrix_U(params.k);
        V = svds.matrix_V(params.k);
        flip_UV(U, V);
        cout << timestamp() << "Final SVD done!\n";
        evals = svals.array().square() / data->nsnps;
        // write to files;
        data->write_eigs_files(evals, U);

    } else {
        // for blockwise
        MatrixXf T, VT;
        ArnoldiOpData *op = new ArnoldiOpData(data);
        SymEigsSolver< float, LARGEST_ALGE, ArnoldiOpData > *eigs = new SymEigsSolver< float, LARGEST_ALGE, ArnoldiOpData >(op, params.k, params.ncv);
        if (params.maxiter == 0) op->setFlags(false, true, false);
        eigs->init();
        nconv = eigs->compute(params.imaxiter, params.itol);
        if (nconv < params.k) {

            params.verbose && cerr << "Warning: the nconv is not equal to k.\n";
        }
        nu = min(params.k, nconv);
        if(eigs->info() == Spectra::SUCCESSFUL)
        {
            op->U = eigs->eigenvectors().leftCols(nu);
            if (params.maxiter == 0)
            {
                evals = eigs->eigenvalues() / data->nsnps;
                data->write_eigs_files(evals, op->U);
                cerr << timestamp() << "Final SVD done!\n";
                return;
            }
            op->S = eigs->eigenvalues().cwiseSqrt();
            // V = G' * U / s
            // VT = (U' / s) * G
            // T = U' / s
            // T' = (U.array().rowwise() / eigs.eigenvalues().head(nv).transpose().array().sqrt()).matrix();
            T = (eigs->eigenvectors().leftCols(nu).transpose().array().colwise() / eigs->eigenvalues().head(nu).array().sqrt()).matrix();
            op->VT = MatrixXf::Zero(T.rows(), data->nsnps);
            data->calcu_vt_initial(T, op->VT);
            flip_UV(op->U, op->VT);
            cerr << timestamp() << "Begin to do EM\n";
            op->setFlags(true, false, params.pcangsd);
            for (uint i = 1; i <= params.maxiter; ++i)
            {
                VT = op->VT;
                eigs->init();
                nconv = eigs->compute(params.imaxiter, params.itol);
                if (nconv < params.k) {
                    params.verbose && cerr << "Warning: the nconv is not equal to k.\n";
                }
                if(eigs->info() == Spectra::SUCCESSFUL)
                {
                    cout << op->S.transpose() << endl;
                    nu = min(params.k, nconv);
                    T = (eigs->eigenvectors().leftCols(nu).transpose().array().colwise() / eigs->eigenvalues().head(nu).array().sqrt()).matrix();

                    data->calcu_vt_update(T, op->U, op->S, op->VT, false);
                    op->S = eigs->eigenvalues().cwiseSqrt();
                    op->U = eigs->eigenvectors().leftCols(nu);
                    flip_UV(op->U, op->VT);
                    diff = rmse(op->VT, VT);
                    params.verbose && cout << timestamp() << "Individual allele frequencies estimated (iter=" << i << "), RMSE=" << diff <<".\n";
                    if (diff < params.tol)
                    {
                        cout << timestamp() << "Come to convergence!\n";
                        break;
                    }

                } else {
                    cerr << "Error: something wrong with Spectra SymEigsSolver\n";
                    exit(EXIT_FAILURE);
                }
            }

            cout << timestamp() << "Begin to standardize the matrix\n";
            op->setFlags(true, true, params.pcangsd);
            eigs->init();
            nconv = eigs->compute(params.imaxiter, params.itol);
            if (nconv < params.k) {
                params.verbose && cerr << "Warning: the nconv is not equal to k.\n";
            }
            if(eigs->info() == Spectra::SUCCESSFUL)
            {
                nu = min(params.k, nconv);
                T = (eigs->eigenvectors().leftCols(nu).transpose().array().colwise() / eigs->eigenvalues().head(nu).array().sqrt()).matrix();

                data->calcu_vt_update(T, op->U, op->S, op->VT, true);
                op->U = eigs->eigenvectors().leftCols(nu);
                flip_UV(op->U, op->VT);
                evals = eigs->eigenvalues() / data->nsnps;
                data->write_eigs_files(evals, op->U);
            } else {
                cerr << "Error: something wrong with Spectra SymEigsSolver\n";
                exit(EXIT_FAILURE);
            }

        } else {
            cerr << "Error: something wrong with Spectra SymEigsSolver\n";
            exit(EXIT_FAILURE);
        }
        delete op;
        delete eigs;
    }
    cout << timestamp() << "eigenvecs and eigenvals are saved. have a nice day. bye!\n";
}
