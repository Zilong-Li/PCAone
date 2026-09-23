/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Arnoldi.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Arnoldi.hpp"

#include <Spectra/contrib/PartialSVDSolver.h>
#include <Spectra/SymEigsSolver.h>

#include "Cmd.hpp"
#include "Utils.hpp"

using namespace std;
using namespace Spectra;

void ArnoldiOpData::perform_op(const double* x_in, double* y_out) const {
  if (data->params.verbose > 1) cao.print(tick.date(), "Arnoldi Matrix Operation =", data->nops);
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
    bool standardized = !(params.missme || params.ld);
    if (standardized) data->standardize_E();
    nconv = svds.compute(params.imaxiter, params.itol);
    if (nconv != params.k) cao.error("the nconv is not equal to k.");
    U = svds.matrix_U(params.k);
    V = svds.matrix_V(params.k);
    svals = svds.singular_values();
    evals.noalias() = svals.array().square().matrix() / data->nsnps;
    // impute information via EM-PCA
    if (params.missme) {
      if (data->p_miss == 0.0) cao.warn("there is no missing values");
      flip_UV(U, V);
      cao.print(tick.date(), "starts EM iteration. maxiter =", params.maxiter);
      for (uint i = 1; i <= params.maxiter; ++i) {
        data->fit_with_pi(U, svals, V.transpose());
        nconv = svds.compute(params.imaxiter, params.itol);
        svals = svds.singular_values();
        U = svds.matrix_U(params.k);
        V2 = svds.matrix_V(params.k);
        flip_UV(U, V2);
        // Same measure as run_pca_with_halko() and the FULL solver, so that a
        // given --tol-em means the same thing whichever -d the user picked.
        // This used to be rmse(), which divides by sqrt(nsnps * k) and so hit
        // any fixed tolerance far earlier -- IRAM stopped after 4 EM
        // iterations where the other solvers took 10, and settled short of the
        // EM fixed point despite being the more accurate decomposition.
        if (params.mev)
          diff = 1.0 - mev(V2, V);
        else
          diff = minSSE(V2, V).sum() / V.cols();
        if (params.verbose)
          cao.print(tick.date(), "individual allele frequencies estimated (iter =", i, "), diff =", diff);
        V = V2;
        if (diff < params.tolem) {
          cao.print(tick.date(), "come to convergence!");
          break;
        }
      }

      if (params.pcangsd) {
        cao.print(tick.date(), "estimate GRM for pcangsd");
        data->pcangsd_standardize_E(U, svals, V.transpose());
        evals.noalias() = svals.array().square().matrix() / data->nsnps;
        if (params.file_t == FileType::BEAGLE) {
          Mat2D C = data->G * data->G.transpose();
          C.array() /= (double)data->nsnps;
          C.diagonal() = data->Dc.array() / (double)data->nsnps;
          std::ofstream fcov(params.fileout + ".cov");
          if (fcov.is_open()) fcov << C << "\n";
          Eigen::JacobiSVD<Mat2D> svd(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
          // output real eigenvectors of covariance in eigvecs2
          write_eigvecs2_beagle(svd.matrixU(), params.filein, params.fileout + ".eigvecs2");
        }
      }

      if (params.emu) {
        cao.print(tick.date(), "standardize the final matrix");
        standardized = true;
        data->standardize_E();
        svds.compute(params.imaxiter, params.itol);
        svals = svds.singular_values();
        U = svds.matrix_U(params.k);
        V = svds.matrix_V(params.k);
        flip_UV(U, V);
        evals.noalias() = svals.array().square().matrix() / data->nsnps;
      }
    }
    // write to files; NOTE: pcangsd only gives us evals of covariance matrix
    if (params.ld && !params.pcangsd) data->write_residuals(svals, U, V.transpose());
    data->set_svd_transform(standardized);
    data->write_eigs_files(evals, svals, U, V);
  } else {
    // for blockwise
    ArnoldiOpData* op = new ArnoldiOpData(data);
    // SymEigsSolver< double, LARGEST_ALGE, ArnoldiOpData > *eigs = new
    // SymEigsSolver< double, LARGEST_ALGE, ArnoldiOpData >(op, params.k,
    // params.ncv);
    SymEigsSolver<ArnoldiOpData>* eigs = new SymEigsSolver<ArnoldiOpData>(*op, params.k, params.ncv);
    // write_residuals() re-reads the blocks unstandardized, so under --ld the
    // decomposition it subtracts has to be unstandardized too. This is the
    // convention in-core IRAM (Arnoldi.cpp above) and Halko already follow.
    bool standardized = !(params.missme || params.ld);
    op->setFlags(false, standardized);
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
    U = (eigs->eigenvectors().leftCols(nu).transpose().array().colwise() / eigs->eigenvalues().head(nu).array().sqrt())
            .matrix();
    op->VT = Mat2D::Zero(U.rows(), data->nsnps);
    data->calcu_vt_initial(U, op->VT, standardized);
    evals.noalias() = eigs->eigenvalues() / data->nsnps;
    // impute information via EM-PCA
    if (params.missme) {
      if (data->p_miss == 0.0) cao.warn("there is no missing values");
      cao.print(tick.date(), "starts EM iteration. maxiter =", params.maxiter);
      data->calcu_vt_initial(U, op->VT, false);
      flip_UV(op->U, op->VT);
      op->setFlags(true, false);
      for (uint i = 1; i <= params.maxiter; ++i) {  // in-core runs maxiter, not maxiter - 1
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
        // VT is k x nsnps here, so transpose before measuring like the in-core
        // loop above does. The stopping threshold is --tol-em; this compared
        // against params.tol (--tol-rsvd, and ten times looser by default),
        // which left --tol-em with no effect at all on out-of-core IRAM.
        if (params.mev)
          diff = 1.0 - mev(op->VT.transpose(), V.transpose());
        else
          diff = minSSE(op->VT.transpose(), V.transpose()).sum() / op->VT.rows();
        if (params.verbose)
          cao.print(tick.date(), "individual allele frequencies estimated (iter =", i, "), diff =", diff);
        if (diff < params.tolem) {
          cao.print(tick.date(), "come to convergence!");
          break;
        }
      }

      standardized = !params.ld;
      op->setFlags(true, standardized);
      eigs->init();
      nconv = eigs->compute(SortRule::LargestAlge, params.imaxiter, params.itol);
      if (nconv < params.k) cao.error("the nconv is not equal to k.");
      assert(eigs->info() == CompInfo::Successful);
      nu = min(params.k, nconv);
      U = (eigs->eigenvectors().leftCols(nu).transpose().array().colwise() /
           eigs->eigenvalues().head(nu).array().sqrt())
              .matrix();

      data->calcu_vt_update(U, op->U, op->S, op->VT, standardized);
      // S comes from this final solve as well. op->U, op->VT and evals all do,
      // so keeping the previous iteration's S left .sigvals disagreeing with
      // .eigvals and with the U/V written beside it.
      op->S = eigs->eigenvalues().cwiseSqrt();
      op->U = eigs->eigenvectors().leftCols(nu);
      flip_UV(op->U, op->VT);
      evals.noalias() = eigs->eigenvalues() / data->nsnps;
    }

    if (params.ld && !params.pcangsd) data->write_residuals(op->S, op->U, op->VT);
    // `standardized` is the flag the solve that produced op->U/S/VT ran with
    data->set_svd_transform(standardized);
    data->write_eigs_files(evals, op->S, op->U, op->VT.transpose());

    delete op;
    delete eigs;
  }

  cao.print(tick.date(), "PCAone - IRAM SVD done");
  return;
}
