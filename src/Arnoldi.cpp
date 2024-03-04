#include "Arnoldi.hpp"

#include "LD.hpp"
#include "Utils.hpp"
#include <Spectra/SymEigsSolver.h>
#include <Spectra/contrib/PartialSVDSolver.h>

using namespace std;
using namespace Spectra;

void ArnoldiOpData::perform_op(const double * x_in, double * y_out) const
{
    if(data->params.verbose) cao << tick.date() << "Arnoldi Matrix Operation = " << data->nops << endl;
    Eigen::Map<const MyVector> x(x_in, n);
    Eigen::Map<MyVector> y(y_out, n);
    auto t1 = std::chrono::high_resolution_clock::now();
    data->check_file_offset_first_var();
    if(update)
    {
        data->read_block_update(data->start[0], data->stop[0], U, S, VT, standardize);
    }
    else
    {
        data->read_block_initial(data->start[0], data->stop[0], standardize);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    data->readtime += std::chrono::duration<double>(t2 - t1).count()
                      * std::chrono::duration<double>::period::num
                      / std::chrono::duration<double>::period::den;

    y.noalias() = data->G * (data->G.transpose() * x);
    for(uint k = 1; k < data->nblocks; ++k)
    {
        auto t1 = std::chrono::high_resolution_clock::now();
        if(update)
        {
            data->read_block_update(data->start[k], data->stop[k], U, S, VT, standardize);
        }
        else
        {
            data->read_block_initial(data->start[k], data->stop[k], standardize);
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        data->readtime += std::chrono::duration<double>(t2 - t1).count()
                          * std::chrono::duration<double>::period::num
                          / std::chrono::duration<double>::period::den;
        // TODO: Kahan summation
        // optimal evaluation see https://eigen.tuxfamily.org/dox/TopicWritingEfficientProductExpression.html
        y.noalias() += data->G * (data->G.transpose() * x);
    }
    data->nops++;
}

void run_pca_with_arnoldi(Data * data, const Param & params)
{
    if(params.out_of_core)
        cao << tick.date() << "running PCAone IRAM with out-of-core mode" << endl;
    else
        cao << tick.date() << "running PCAone IRAM with in-core mode" << endl;
    MyMatrix U, V, V2;
    MyVector svals, evals;
    uint nconv, nu;
    double diff;
    if(!params.out_of_core)
    {
        // SpMatrix sG = data->G.sparseView();
        PartialSVDSolver<MyMatrix> svds(data->G, params.k, params.ncv);
        if(!params.ld && !params.runem
           && (params.file_t == FileType::PLINK || params.file_t == FileType::BGEN))
            data->standardize_E();
        nconv = svds.compute(params.imaxiter, params.itol);
        if(nconv != params.k) cao.error("the nconv is not equal to k.");
        U = svds.matrix_U(params.k);
        V = svds.matrix_V(params.k);
        svals = svds.singular_values();
        evals.noalias() = svals.array().square().matrix() / data->nsnps;
        if(params.runem)
        {
            flip_UV(U, V);
            cao << tick.date() << "begin to do EM iteration.\n";
            for(uint i = 1; i <= params.maxiter; ++i)
            {
                data->update_batch_E(U, svals, V.transpose());
                nconv = svds.compute(params.imaxiter, params.itol);
                svals = svds.singular_values();
                U = svds.matrix_U(params.k);
                V2 = svds.matrix_V(params.k);
                flip_UV(U, V2);
                diff = rmse(V2, V);
                if(params.verbose)
                    cao << tick.date() << "individual allele frequencies estimated (iter=" << i
                        << "), RMSE=" << diff << ".\n";
                V = V2;
                if(diff < params.tolem)
                {
                    cao << tick.date() << "come to convergence!\n";
                    break;
                }
            }

            cao << tick.date() << "begin to standardize the matrix\n";
            if(params.pcangsd)
            {
                data->pcangsd_standardize_E(U, svals, V.transpose());
                MyMatrix C = data->G * data->G.transpose();
                C.array() /= (double)data->nsnps;
                C.diagonal() = data->Dc.array() / (double)data->nsnps;
                std::ofstream out_cov(params.fileout + ".cov");
                if(out_cov.is_open()) out_cov << C << "\n";
                // calculate eigenvectors
                DenseSymMatProd<double> op2(C);
                SymEigsSolver<DenseSymMatProd<double>> eigs2(op2, params.k, params.ncv);
                eigs2.init();
                nconv = eigs2.compute(SortRule::LargestAlge, params.imaxiter, params.itol);
                assert(eigs2.info() == CompInfo::Successful);
                nu = min(params.k, nconv);
                U = eigs2.eigenvectors().leftCols(nu);
                V = eigs2.eigenvectors().leftCols(nu);
                evals = eigs2.eigenvalues();
            }
            else
            {
                data->standardize_E();
                svds.compute(params.imaxiter, params.itol);
                svals = svds.singular_values();
                U = svds.matrix_U(params.k);
                V = svds.matrix_V(params.k);
                flip_UV(U, V);
                evals.noalias() = svals.array().square().matrix() / data->nsnps;
            }
        }
        // write to files;
        data->write_eigs_files(evals, U, V);
        // NOTE: pcangsd only gives us evals of covariance matrix
        if(params.ld && !params.pcangsd) data->write_residuals(svals, U, V);
    }
    else
    {
        // for blockwise
        ArnoldiOpData * op = new ArnoldiOpData(data);
        // SymEigsSolver< double, LARGEST_ALGE, ArnoldiOpData > *eigs = new SymEigsSolver< double,
        // LARGEST_ALGE, ArnoldiOpData >(op, params.k, params.ncv);
        SymEigsSolver<ArnoldiOpData> * eigs = new SymEigsSolver<ArnoldiOpData>(*op, params.k, params.ncv);
        if(!params.runem) op->setFlags(false, true, false);
        eigs->init();
        nconv = eigs->compute(SortRule::LargestAlge, params.imaxiter, params.itol);
        if(nconv < params.k) cao.error("the nconv is not equal to k");
        nu = min(params.k, nconv);
        assert(eigs->info() == CompInfo::Successful);
        op->U = eigs->eigenvectors().leftCols(nu);
        op->S = eigs->eigenvalues().cwiseSqrt();
        // V = G' * U / s
        // VT = (U' / s) * G
        // T = U' / s
        // T' = (U.array().rowwise() / eigs.eigenvalues().head(nv).transpose().array().sqrt()).matrix();
        // reuse MyMatrix U = T, V = VT;
        U = (eigs->eigenvectors().leftCols(nu).transpose().array().colwise()
             / eigs->eigenvalues().head(nu).array().sqrt())
                .matrix();
        op->VT = MyMatrix::Zero(U.rows(), data->nsnps);
        data->calcu_vt_initial(U, op->VT, true);
        evals.noalias() = eigs->eigenvalues() / data->nsnps;
        if(params.runem)
        {
            data->calcu_vt_initial(U, op->VT, false);
            flip_UV(op->U, op->VT);
            cao << tick.date() << "begin to do EM iteration.\n";
            op->setFlags(true, false, params.pcangsd);
            for(uint i = 1; i <= params.maxiter; ++i)
            {
                V = op->VT;
                eigs->init();
                nconv = eigs->compute(SortRule::LargestAlge, params.imaxiter, params.itol);
                if(nconv < params.k) cao.error("the nconv is not equal to k.");
                assert(eigs->info() == CompInfo::Successful);
                nu = min(params.k, nconv);
                U = (eigs->eigenvectors().leftCols(nu).transpose().array().colwise()
                     / eigs->eigenvalues().head(nu).array().sqrt())
                    .matrix();

                data->calcu_vt_update(U, op->U, op->S, op->VT, false);
                op->S = eigs->eigenvalues().cwiseSqrt();
                op->U = eigs->eigenvectors().leftCols(nu);
                flip_UV(op->U, op->VT);
                diff = rmse(op->VT, V);
                if(params.verbose)
                    cao << tick.date() << "individual allele frequencies estimated (iter=" << i
                        << "), RMSE=" << diff << ".\n";
                if(diff < params.tol)
                    {
                        cao << tick.date() << "come to convergence!\n";
                        break;
                    }
            }

            cao << tick.date() << "begin to standardize the matrix\n";
            op->setFlags(true, true, params.pcangsd);
            eigs->init();
            nconv = eigs->compute(SortRule::LargestAlge, params.imaxiter, params.itol);
            if(nconv < params.k) cao.error("the nconv is not equal to k.");
            assert(eigs->info() == CompInfo::Successful);
            nu = min(params.k, nconv);
            U = (eigs->eigenvectors().leftCols(nu).transpose().array().colwise()
                 / eigs->eigenvalues().head(nu).array().sqrt())
                    .matrix();

            data->calcu_vt_update(U, op->U, op->S, op->VT, true);
            op->U = eigs->eigenvectors().leftCols(nu);
            flip_UV(op->U, op->VT);
            evals.noalias() = eigs->eigenvalues() / data->nsnps;
        }
        data->write_eigs_files(evals, op->U, op->VT.transpose());
        delete op;
        delete eigs;
    }
    
    cao << tick.date() << "run_pca_with_arnoldi done\n";

    return;
}
