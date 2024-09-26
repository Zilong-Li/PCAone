// clang-off
#define _DECLARE_TOOLBOX_HERE

#include <omp.h>

#include "Arnoldi.hpp"
#include "Cmd.hpp"
#include "Data.hpp"
#include "FileBeagle.hpp"
#include "FileBgen.hpp"
#include "FileBinary.hpp"
#include "FileCsv.hpp"
#include "FilePlink.hpp"
#include "Halko.hpp"
#include "LD.hpp"

#ifdef WITH_OPENBLAS
#include "lapacke.h"
#elif defined WITH_MKL
#include "mkl_lapacke.h"
#endif

// clang-on
using namespace std;

int main(int argc, char* argv[]) {
  Param params(argc, argv);
  cao.cao.open(params.fileout + ".log");
  string commandargs = params.ss.str();
  cao << get_machine() << commandargs << endl;
  // set number of threads
  // openblas_set_num_threads(params.threads);
  omp_set_num_threads(params.threads);
  Data* data = nullptr;

  // particular case for LD
  if ((params.file_t == FileType::BINARY) &&
      (params.print_r2 || params.ld_r2 > 0 || !params.clump.empty())) {
    data = new FileBin(params);
    data->prepare();
    run_ld_stuff(params, data);
    delete data;
    cao.print(tick.date(), "total elapsed wall time:", tick.abstime(),
              "seconds");
    return 0;
  }

  if (params.perm && params.out_of_core) {
    tick.clock();
    if (params.file_t == FileType::PLINK) {
      auto perm = permute_plink(params.filein, params.fileout, params.buffer,
                                params.bands);
      data = new FileBed(params);
      data->perm = perm;
    } else if (params.file_t == FileType::BGEN) {
      auto perm = permute_bgen(params.filein, params.fileout, params.threads);
      data = new FileBgen(params);
      data->perm = perm;
    } else if (params.file_t == FileType::CSV) {
      auto perm = shuffle_csvzstd_to_bin(params.filein, params.fileout,
                                         params.buffer, params.scale);
      params.file_t = FileType::BINARY;
      data = new FileBin(params);
      data->perm = perm;
    } else {
      cao.error("wrong file type used!");
    }
    cao.print(tick.date(), "elapsed time of permuting data:", tick.reltime(),
              "seconds");
  } else {
    if (params.file_t == FileType::PLINK) {
      data = new FileBed(params);
    } else if (params.file_t == FileType::BGEN) {
      data = new FileBgen(params);
    } else if (params.file_t == FileType::BEAGLE) {
      data = new FileBeagle(params);
    } else if (params.file_t == FileType::BINARY) {
      data = new FileBin(params);
    } else if (params.file_t == FileType::CSV) {
      data = new FileCsv(params);
    } else {
      cao.error("invalid input files!");
    }
  }

  // be prepared for run
  data->prepare();

  // begin to run PCA
  if (params.svd_t == SvdType::IRAM) {
    run_pca_with_arnoldi(data, params);
  } else if (params.svd_t == SvdType::PCAoneAlg1 ||
             params.svd_t == SvdType::PCAoneAlg2) {
    run_pca_with_halko(data, params);
  } else if (params.svd_t == SvdType::FULL) {
    cao.print(tick.date(), "running the Full SVD with in-core mode.");
    if (params.file_t == FileType::PLINK || params.file_t == FileType::BGEN)
      data->standardize_E();
    Eigen::JacobiSVD<MyMatrix> svd(data->G,
                                   Eigen::ComputeThinU | Eigen::ComputeThinV);
    data->write_eigs_files(svd.singularValues().array().square() / data->nsnps,
                           svd.matrixU(), svd.matrixV());
  } else {
    cao.error("unsupported PCA method!");
  }

  delete data;

  if (params.file_t == FileType::PLINK)
    make_plink2_eigenvec_file(params.k, params.fileout + ".eigvecs2",
                              params.fileout + ".eigvecs",
                              params.filein + ".fam");

  cao.print(tick.date(), "total elapsed reading time: ", data->readtime,
            "seconds");
  cao.print(tick.date(), "total elapsed wall time:", tick.abstime(), "seconds");
  cao.print(tick.date(),
            "eigenvecs and eigenvals are saved. have a nice day. bye!");

  return 0;
}
