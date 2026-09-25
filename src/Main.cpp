// clang-off
#define _DECLARE_TOOLBOX_HERE

#include <omp.h>

#include <thread>

#include "Arnoldi.hpp"
#include "Cmd.hpp"
#include "Common.hpp"
#include "Data.hpp"
#include "EvalAdmix.hpp"
#include "Exact.hpp"
#include "FileBeagle.hpp"
#include "FileBgen.hpp"
#include "FileBinary.hpp"
#include "FileCsv.hpp"
#include "FilePgen.hpp"
#include "FilePlink.hpp"
#include "FileUSV.hpp"
#include "Halko.hpp"
#include "InbredSamples.hpp"
#include "InbredSites.hpp"
#include "LD.hpp"
#include "Projection.hpp"
#include "Selection.hpp"

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
  // setenv alone does not bound OpenMP: the runtime reads OMP_NUM_THREADS when
  // it initializes, before main() runs, so setting it here is too late and -n
  // was silently ignored -- every parallel region used the whole machine. MKL
  // and OpenBLAS read their variables lazily, so those two do take effect.
  omp_set_num_threads(threadCount > 0 ? threadCount : 1);
}

static int run(int argc, char* argv[]) {
  Param params(argc, argv);
  cao.cao.open(params.fileout + ".log");
  if (params.verbose > 0) cao.is_screen = true;
  cao.print(get_machine(), params.ss.str());
  // limit the number of threads
  uint max_threads = std::thread::hardware_concurrency();
  max_threads = params.threads > max_threads ? max_threads : params.threads;
  setEnvironmentVariables(max_threads);
  cao.print(tick.date(), "program started with " + std::to_string(max_threads) + " threads");
  Data* data = nullptr;

  // particular case for inbreeding sites
  if (params.inbreed == 1) {
    data = new FileUSV(params);
    run_inbred_sites(data, params);
    delete data;
    return bye();
  }

  // particular case for inbreeding samples
  // if (params.inbreed == 2) {
  //   data = new FileUSV(params);
  //   run_inbreed_coef_sample(data, params);
  //   delete data;
  //   return bye();
  // }

  // particular case for LD: R2, pruning and clumping on the genotypes, with the
  // PCs of -P/--USV removed as they are read. Cmd.cpp allows only PLINK/PGEN
  // here and turns off the shuffling, as LD walks the sites in file order.
  if (params.ld) {
    if (params.file_t == FileType::PLINK)
      data = new FileBed(params);
    else
      data = new FilePgen(params);
    run_ld_stuff(data, params);
    delete data;
    return bye();
  }

  // particular case for projection
  if (params.project > 0 &&
      (params.file_t == FileType::PLINK || params.file_t == FileType::BEAGLE || params.file_t == FileType::PGEN)) {
    if (params.file_t == FileType::PLINK)
      data = new FileBed(params);
    else if (params.file_t == FileType::BEAGLE)
      data = new FileBeagle(params);
    else
      data = new FilePgen(params);
    run_projection(data, params);
    delete data;
    return bye();
  }

  // particular case for Selection
  if ((params.selection > 0) && (params.file_t == FileType::PLINK || params.file_t == FileType::PGEN)) {
    // perm must be off, and here it is not merely safe but required.
    //
    // params.perm is on by default, since the default --svd 2 is winSVD and
    // --no-shuffle is off. But --selection decomposes nothing: it reads a
    // reference U and passes over the genotypes once. The permutation this flag
    // announces is only built further down, past this early return, so nothing
    // initializes it here. For PGEN the permutation is *logical* -- FilePgen
    // maps every read through perm.indices() -- so read_block_initial() indexes
    // an empty permutation and segfaults out-of-core. PLINK escaped it only
    // because its permutation lives in a temp .bed that --selection never
    // writes, leaving it to read the original file in stored order.
    //
    // Unlike --evaladmix, whose statistic is a sum over sites, site order
    // matters here: the rows of .zscore / .galinsky / .pcadapt* are positional
    // against the .bim/.pvar. So the fix is to clear the flag rather than to
    // build a permutation.
    params.perm = false;
    if (params.file_t == FileType::PLINK)
      data = new FileBed(params);
    else
      data = new FilePgen(params);
    run_selection(data, params);
    delete data;
    return bye();
  }

  const bool ooc_permutation = params.perm && params.out_of_core;

  if (ooc_permutation) {
    tick.clock();
    if (params.file_t == FileType::PLINK) {
      auto perm = permute_plink(params.filein, params.fileout, params.buffer, params.bands);
      data = new FileBed(params);
      data->perm = perm;
    } else if (params.file_t == FileType::PGEN) {
      // Logical permutation is initialized after prepare(), when blocksize is known.
      data = new FilePgen(params);
    } else if (params.file_t == FileType::BGEN) {
      auto perm = permute_bgen(params.filein, params.fileout, params.threads, params.seed);
      data = new FileBgen(params);
      data->perm = perm;
    } else if (params.file_t == FileType::CSV) {
      auto perm = shuffle_csvzstd_to_bin(params.filein, params.fileout, params.buffer, params.scale, params.scaleFactor,
                                         params.seed);
      params.file_t = FileType::BINARY;
      data = new FileBin(params);
      data->perm = perm;
    } else {
      cao.error("wrong file type used!");
    }
    if (params.file_t != FileType::PGEN)
      cao.print(tick.date(), "elapsed time of permuting data:", tick.reltime(), " seconds");
  } else {
    if (params.file_t == FileType::PLINK) {
      data = new FileBed(params);
    } else if (params.file_t == FileType::PGEN) {
      data = new FilePgen(params);
    } else if (params.file_t == FileType::BGEN) {
      data = new FileBgen(params);
    } else if (params.file_t == FileType::BEAGLE) {
      data = new FileBeagle(params);
    } else if (params.file_t == FileType::CSV) {
      data = new FileCsv(params);
    } else {
      cao.error("invalid input files!");
    }
  }

  // --svd 3 with no more samples than sites needs only the N x N GRM, so it
  // streams the genotypes block by block instead of holding G (N x M doubles):
  // 400 x 656,281 took 2.1 GB, now 0.11 GB. -m sets the block size. --maf
  // filters only in-core, and N > M decomposes the M x M G'G in-core.
  if (params.svd_t == SvdType::FULL) {
    if (!params.filterSNP && data->nsamples <= data->nsnps)
      params.out_of_core = true;
    else if (params.memory > 0)
      cao.warn("-m is ignored: --svd 3 runs in-core with --maf, or with more samples than sites");
  }

  // be prepared for run
  data->prepare();
  // every SVD below indexes the first sample and site
  if (data->nsamples == 0 || data->nsnps == 0)
    cao.error("the input has " + std::to_string(data->nsamples) + " samples and " + std::to_string(data->nsnps) +
              " sites; nothing to decompose");
  // At most min(N, M) - 1 PCs of centred data carry anything: more came back as
  // arbitrary vectors from the null space, and different on every run.
  const uint rank_max = std::min(data->nsamples, data->nsnps);
  if (params.k >= rank_max)
    cao.error("-k/--pc " + std::to_string(params.k) + " must be smaller than the number of samples and of sites (" +
              std::to_string(data->nsamples) + " and " + std::to_string(data->nsnps) + ")");
  // the RSVD cannot use more test vectors, nor the IRAM more Lanczos vectors, than that
  if (params.k + params.oversamples > rank_max) {
    params.oversamples = rank_max - params.k;
    cao.warn("--oversamples reduced to " + std::to_string(params.oversamples) +
             " so that k + oversamples <= " + std::to_string(rank_max));
  }
  if (params.ncv > rank_max) {
    params.ncv = rank_max;
    cao.warn("--ncv reduced to " + std::to_string(params.ncv) + ", the size of the problem");
  }
  if (ooc_permutation && params.file_t == FileType::PGEN) {
    // The permutation draws each block from a fixed number of shuffled source
    // regions, one per reading thread at that many threads. It used the thread
    // count, so the same --seed gave different PCs with a different -n.
    const uint regions = 16;
    data->perm = compute_pgen_perm(data->nsnps, params.bands, data->blocksize, regions, params.seed);
    cao.print(tick.date(), "initialized logical PGEN permutation. blocksize:", data->blocksize,
              ", batches:", params.bands, ", source regions:", regions);
  }

  // begin to run PCA
  if (params.svd_t == SvdType::IRAM) {
    run_pca_with_arnoldi(data, params);
  } else if (params.svd_t == SvdType::PCAoneAlg1 || params.svd_t == SvdType::PCAoneAlg2) {
    run_pca_with_halko(data, params);
  } else if (params.svd_t == SvdType::FULL) {
    run_pca_exact(data, params);
  } else {
    cao.error("unsupported PCA method!");
  }

  cao.print(tick.date(), "total elapsed reading time: ", data->readtime, " seconds");

  delete data;

  // evalAdmix: correlation of residuals given the PCs just computed.
  // needs a second pass over the raw (uncentered, unstandardized) genotypes.
  if (params.evaladmix) {
    if (params.file_t != FileType::PLINK && params.file_t != FileType::PGEN)
      cao.error("--evaladmix currently supports PLINK bed/pgen input only");
    // Second pass over the genotypes: centred, never standardized. Param is
    // not copyable, so mutate in place -- the PCA is finished by now.
    //
    // center must stay TRUE. With center=false, read_all() leaves missing
    // calls as BED_MISSING_VALUE (-9), because the imputation in FileBed /
    // FilePgen sits inside `if (params.center)`. Those -9s poison G*G' and
    // drive the heterozygosity accumulator negative, so the statistic comes
    // out saturated at the +-1 clamp. Centring costs nothing: the statistic is
    // invariant to it (see run_evaladmix), and missing calls are then imputed
    // to the site mean, the convention used everywhere else in PCAone.
    //
    // dopca stays on so the same --maf filter selects the sites the PCA used,
    // and read_all() / read_block_initial(.., false) never standardize.
    params.center = true;
    // perm must be off. For PGEN the permutation is *logical*: FilePgen maps
    // every read through perm.indices(), and Main only initializes that on the
    // first Data object (above, once its blocksize is known). d2 would index an
    // empty permutation and segfault. For PLINK the permutation is already
    // baked into the temp .bed on disk, so clearing the flag simply reads that
    // file in stored order. Either way the site order does not matter here: A,
    // b and d are all sums over sites.
    params.perm = false;
    Data* d2 = (params.file_t == FileType::PLINK) ? (Data*)new FileBed(params) : (Data*)new FilePgen(params);
    d2->prepare();
    run_evaladmix(d2, params);
    delete d2;
  }

  if (params.file_t == FileType::PLINK)
    make_plink2_eigenvec_file(params.k, params.fileout + ".eigvecs2", params.fileout + ".eigvecs",
                              params.filein + ".fam");
  else if (params.file_t == FileType::PGEN)
    make_plink2_eigenvec_from_psam(params.k, params.fileout + ".eigvecs2", params.fileout + ".eigvecs",
                                   params.filein + ".psam");

  // remove temp files if verbose < 3
  if (ooc_permutation && params.verbose < 3) {
    if (params.file_t == FileType::PLINK) {
      for (auto suf : std::vector<std::string>{".bed", ".bim", ".fam"}) {
        std::filesystem::path tmpfile{params.filein + suf};
        std::filesystem::remove(tmpfile);
      }
    }
    // a shuffled CSV is read back as BINARY, from <out>.perm.bin
    if ((params.file_t == FileType::BGEN) || (params.file_t == FileType::CSV) ||
        (params.file_t == FileType::BINARY)) {
      std::filesystem::path tmpfile{params.filein};
      std::filesystem::remove(tmpfile);
    }
  }
  return bye();
}

// cao.error() throws, and so do Eigen (std::bad_alloc) and the standard library.
// Uncaught, each of them aborted the run with "terminate called after throwing
// ..." and exit code 134, which looks like a crash to the user and to workflow
// managers. Report the message and exit with 1 instead.
int main(int argc, char* argv[]) {
  try {
    return run(argc, argv);
  } catch (const std::bad_alloc&) {
    std::cerr << "Error: out of memory (std::bad_alloc). consider -m/--memory for the out-of-core mode\n";
  } catch (const std::exception& e) {
    std::string msg = e.what();
    if (msg.empty() || msg.back() != '\n') msg += '\n';
    std::cerr << "Error: " << msg;
  }
  return EXIT_FAILURE;
}
