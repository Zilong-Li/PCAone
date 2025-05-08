// clang-off
#define _DECLARE_TOOLBOX_HERE

#include <omp.h>
#include <thread>

#include "Arnoldi.hpp"
#include "Cmd.hpp"
#include "Common.hpp"
#include "Data.hpp"
#include "FileBeagle.hpp"
#include "FileBgen.hpp"
#include "FileBinary.hpp"
#include "FileCsv.hpp"
#include "FilePlink.hpp"
#include "FileUSV.hpp"
#include "Halko.hpp"
#include "Inbreeding.hpp"
#include "LD.hpp"
#include "Projection.hpp"

#ifdef WITH_OPENBLAS
#include "lapacke.h"
#elif defined WITH_MKL
#include "mkl_lapacke.h"
#endif

// clang-on
using namespace std;

static int bye() {
  cao.print(tick.date(), "total elapsed wall time:", tick.abstime(), " seconds");
  cao.print(tick.date(), "have a nice day. bye!");
  return 0;
}

static void setEnvironmentVariables(int threadCount) {
  std::string countStr = std::to_string(threadCount);
  setenv("MKL_NUM_THREADS", countStr.c_str(), 1);
  setenv("OMP_NUM_THREADS", countStr.c_str(), 1);
  setenv("OPENBLAS_NUM_THREADS", countStr.c_str(), 1);
}

int main(int argc, char* argv[]) {
  Param params(argc, argv);
  cao.cao.open(params.fileout + ".log");
  if (params.verbose > 0) cao.is_screen = true;
  cao.print(get_machine(), params.ss.str());
  // limit the number of threads
  uint max_threads = std::thread::hardware_concurrency();
  max_threads = params.threads > max_threads ? max_threads : params.threads;
  setEnvironmentVariables(max_threads);
  cao.print(tick.date(), "program started");
  Data* data = nullptr;

  // particular case for inbreeding
  if (params.inbreed > 0) {
    data = new FileUSV(params);
    run_inbreeding(data, params);
    delete data;
    return bye();
  }

  // particular case for LD
  if (((params.file_t == FileType::BINARY || ((params.file_t == FileType::PLINK) && !params.fileU.empty())) &&
       (params.print_r2 || params.ld_r2 > 0 || !params.clump.empty()))) {
    if (params.filebim.empty()) params.filebim = params.filein + ".bim";
    if (params.file_t == FileType::BINARY)
      data = new FileBin(params);
    else {
      params.memory = 0, params.out_of_core = false;
      data = new FileBed(params);
    }
    run_ld_stuff(data, params);
    delete data;
    return bye();
  }

  // particular case for projection
  if ((params.project > 0) && (params.file_t == FileType::PLINK)) {
    data = new FileBed(params);
    run_projection(data, params);
    delete data;
    return bye();
  }

  if (params.perm && params.out_of_core) {
    tick.clock();
    if (params.file_t == FileType::PLINK) {
      auto perm = permute_plink(params.filein, params.fileout, params.buffer, params.bands);
      data = new FileBed(params);
      data->perm = perm;
    } else if (params.file_t == FileType::BGEN) {
      auto perm = permute_bgen(params.filein, params.fileout, params.threads);
      data = new FileBgen(params);
      data->perm = perm;
    } else if (params.file_t == FileType::CSV) {
      auto perm = shuffle_csvzstd_to_bin(params.filein, params.fileout, params.buffer, params.scale);
      params.file_t = FileType::BINARY;
      data = new FileBin(params);
      data->perm = perm;
    } else {
      cao.error("wrong file type used!");
    }
    cao.print(tick.date(), "elapsed time of permuting data:", tick.reltime(), " seconds");
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
  } else if (params.svd_t == SvdType::PCAoneAlg1 || params.svd_t == SvdType::PCAoneAlg2) {
    run_pca_with_halko(data, params);
  } else if (params.svd_t == SvdType::FULL) {
    if (params.file_t == FileType::PLINK || params.file_t == FileType::BGEN) data->standardize_E();
    cao.print(tick.date(), "running the Full SVD with in-core mode.");
    Eigen::JacobiSVD<Mat2D> svd(data->G, Eigen::ComputeThinU | Eigen::ComputeThinV);
    data->write_eigs_files(svd.singularValues().array().square() / data->nsnps, svd.matrixU(), svd.matrixV());
  } else {
    cao.error("unsupported PCA method!");
  }

  cao.print(tick.date(), "total elapsed reading time: ", data->readtime, " seconds");

  delete data;

  if (params.file_t == FileType::PLINK)
    make_plink2_eigenvec_file(params.k, params.fileout + ".eigvecs2", params.fileout + ".eigvecs",
                              params.filein + ".fam");

  // remove temp files if verbose < 3
  if (params.perm && params.out_of_core && params.verbose < 3) {
    if (params.file_t == FileType::PLINK) {
      for (auto suf : std::vector<std::string>{".bed", ".bim", ".fam"}) {
        std::filesystem::path tmpfile{params.filein + suf};
        std::filesystem::remove(tmpfile);
      }
    }
    if ((params.file_t == FileType::BGEN) || (params.file_t == FileType::CSV)) {
      std::filesystem::path tmpfile{params.filein};
      std::filesystem::remove(tmpfile);
    }
  }
  return bye();
}
