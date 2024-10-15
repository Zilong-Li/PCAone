/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Arnoldi.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Arnoldi.hpp"

#include <Spectra/SymEigsSolver.h>
#include <Spectra/contrib/PartialSVDSolver.h>

#include "Utils.hpp"

using namespace std;
using namespace Spectra;

void ArnoldiOpData::perform_op(const double* x_in, double* y_out) const {
  if (data->params.verbose) cao.print(tick.date(), "Arnoldi Matrix Operation =", data->nops);
  Eigen::Map<const Mat1D> x(x_in, n);
  Eigen::Map<Mat1D> y(y_out, n);
  tick.clock();
  data->check_file_offset_first_var();
  if (update) {
    data->read_block_update(data->start[0], data->stop[0], U, S, VT, standardize);
  } else {
    data->read_block_initial(data->start[0], data->stop[0], standardize);
  }
  data->readtime += tick.reltime();

  y.noalias() = data->G * (data->G.transpose() * x);
  for (uint k = 1; k < data->nblocks; ++k) {
    tick.clock();
    if (update) {
      data->read_block_update(data->start[k], data->stop[k], U, S, VT, standardize);
    } else {
      data->read_block_initial(data->start[k], data->stop[k], standardize);
    }
    data->readtime += tick.reltime();
    // TODO: Kahan summation
    // optimal evaluation see
    // https://eigen.tuxfamily.org/dox/TopicWritingEfficientProductExpression.html
    y.noalias() += data->G * (data->G.transpose() * x);
  }
  data->nops++;
}

void run_pca_with_arnoldi(Data* data, const Param& params) {
  if (params.out_of_core)
    cao.print(tick.date(), "running IRAM SVD with out-of-core mode.");
  else
    cao.print(tick.date(), "running IRAM SVD with in-core mode.");
  Mat2D U, V, V2;
  Mat1D svals, evals;
  uint nconv, nu;
  double diff;
  if (!params.out_of_core) {
    // SpMatrix sG = data->G.sparseView();
    PartialSVDSolver<Mat2D> svds(data->G, params.k, params.ncv);
    if (!params.ld && !params.impute && (params.file_t == FileType::PLINK || params.file_t == FileType::BGEN))
      data->standardize_E();
    nconv = svds.compute(params.imaxiter, params.itol);
    if (nconv != params.k) cao.error("the nconv is not equal to k.");
    U = svds.matrix_U(params.k);
    V = svds.matrix_V(params.k);
    svals = svds.singular_values();
    evals.noalias() = svals.array().square().matrix() / data->nsnps;
    if (params.impute) {
      flip_UV(U, V);
      cao.print(tick.date(), "do EM-PCA algorithms for data with uncertainty.");
      for (uint i = 1; i <= params.maxiter; ++i) {
        data->update_batch_E(U, svals, V.transpose());
        nconv = svds.compute(params.imaxiter, params.itol);
        svals = svds.singular_values();
        U = svds.matrix_U(params.k);
        V2 = svds.matrix_V(params.k);
        flip_UV(U, V2);
        diff = rmse(V2, V);
        if (params.verbose)
          cao.print(tick.date(), "individual allele frequencies estimated (iter =", i, "), RMSE =", diff);
        V = V2;
        if (diff < params.tolem) {
          cao.print(tick.date(), "come to convergence!");
          break;
        }
      }

      if (params.pcangsd) {
        data->pcangsd_standardize_E(U, svals, V.transpose());
        evals.noalias() = svals.array().square().matrix() / data->nsnps;
        Mat2D C = data->G * data->G.transpose();
        C.array() /= (double)data->nsnps;
        C.diagonal() = data->Dc.array() / (double)data->nsnps;
        std::ofstream fcov(params.fileout + ".cov");
        if (fcov.is_open()) fcov << C << "\n";
        Eigen::JacobiSVD<Mat2D> svd(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
        // output real eigenvectors of covariance in eigvecs2
        write_eigvecs2_beagle(svd.matrixU(), params.filein, params.fileout + ".eigvecs2");

        // DenseSymMatProd<double> op2(C);
        // SymEigsSolver<DenseSymMatProd<double>> eigs2(op2, params.k, params.ncv);
        // eigs2.init();
        // nconv = eigs2.compute(SortRule::LargestAlge, params.imaxiter, params.itol);
        // assert(eigs2.info() == CompInfo::Successful);
        // output real eigenvectors of covariance in eigvecs2
        // nu = min(params.k, nconv);
        // U = eigs2.eigenvectors().leftCols(nu);
        // V = eigs2.eigenvectors().leftCols(nu);
        // evals = eigs2.eigenvalues();
        // write_eigvecs2_beagle(eigs2.eigenvectors(), params.filein, params.fileout + ".eigvecs2");
      } else {
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
    if (params.ld && !params.pcangsd) data->write_residuals(svals, U, V.transpose());
  } else {
    // for blockwise
    ArnoldiOpData* op = new ArnoldiOpData(data);
    // SymEigsSolver< double, LARGEST_ALGE, ArnoldiOpData > *eigs = new
    // SymEigsSolver< double, LARGEST_ALGE, ArnoldiOpData >(op, params.k,
    // params.ncv);
    SymEigsSolver<ArnoldiOpData>* eigs = new SymEigsSolver<ArnoldiOpData>(*op, params.k, params.ncv);
    if (!params.impute) op->setFlags(false, true, false);
    eigs->init();
    nconv = eigs->compute(SortRule::LargestAlge, params.imaxiter, params.itol);
    if (nconv < params.k) cao.error("the nconv is not equal to k");
    nu = min(params.k, nconv);
    assert(eigs->info() == CompInfo::Successful);
    op->U = eigs->eigenvectors().leftCols(nu);
    op->S = eigs->eigenvalues().cwiseSqrt();
    // V = G' * U / s
    // VT = (U' / s) * G
    // T = U' / s
    // T' = (U.array().rowwise() /
    // eigs.eigenvalues().head(nv).transpose().array().sqrt()).matrix(); reuse
    // MyMatrix U = T, V = VT;
    U = (eigs->eigenvectors().leftCols(nu).transpose().array().colwise() /
         eigs->eigenvalues().head(nu).array().sqrt())
            .matrix();
    op->VT = Mat2D::Zero(U.rows(), data->nsnps);
    data->calcu_vt_initial(U, op->VT, true);
    evals.noalias() = eigs->eigenvalues() / data->nsnps;
    if (params.impute) {
      cao.print(tick.date(), "starts EM iteration");
      data->calcu_vt_initial(U, op->VT, false);
      flip_UV(op->U, op->VT);
      op->setFlags(true, false, params.pcangsd);
      for (uint i = 1; i <= params.maxiter; ++i) {
        V = op->VT;
        eigs->init();
        nconv = eigs->compute(SortRule::LargestAlge, params.imaxiter, params.itol);
        if (nconv < params.k) cao.error("the nconv is not equal to k.");
        assert(eigs->info() == CompInfo::Successful);
        nu = min(params.k, nconv);
        U = (eigs->eigenvectors().leftCols(nu).transpose().array().colwise() /
             eigs->eigenvalues().head(nu).array().sqrt())
                .matrix();

        data->calcu_vt_update(U, op->U, op->S, op->VT, false);
        op->S = eigs->eigenvalues().cwiseSqrt();
        op->U = eigs->eigenvectors().leftCols(nu);
        flip_UV(op->U, op->VT);
        diff = rmse(op->VT, V);
        if (params.verbose)
          cao.print(tick.date(), "individual allele frequencies estimated (iter =", i, "), RMSE =", diff);
        if (diff < params.tol) {
          cao.print(tick.date(), "come to convergence!");
          break;
        }
      }

      op->setFlags(true, true, params.pcangsd);
      eigs->init();
      nconv = eigs->compute(SortRule::LargestAlge, params.imaxiter, params.itol);
      if (nconv < params.k) cao.error("the nconv is not equal to k.");
      assert(eigs->info() == CompInfo::Successful);
      nu = min(params.k, nconv);
      U = (eigs->eigenvectors().leftCols(nu).transpose().array().colwise() /
           eigs->eigenvalues().head(nu).array().sqrt())
              .matrix();

      data->calcu_vt_update(U, op->U, op->S, op->VT, true);
      op->U = eigs->eigenvectors().leftCols(nu);
      flip_UV(op->U, op->VT);
      evals.noalias() = eigs->eigenvalues() / data->nsnps;
    }

    data->write_eigs_files(evals, op->U, op->VT.transpose());
    if (params.ld && !params.pcangsd) data->write_residuals(op->S, op->U, op->VT);

    delete op;
    delete eigs;
  }

  cao.print(tick.date(), "PCAone - IRAM SVD done");
  return;
}
