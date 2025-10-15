#ifndef CMD_H_
#define CMD_H_

#include <iostream>
#include <iterator>
#include <sstream>
#include <cstdint>

using uint = std::uint32_t;
using uint64 = std::uint64_t;

enum class FileType { PLINK, CSV, BEAGLE, BINARY, BGEN };

enum class SvdType { IRAM, PCAoneAlg1, PCAoneAlg2, FULL };

class Param {
 public:
  Param(int argc, char** argv);
  ~Param();

  FileType file_t;  // PLINK, CSV, BEAGLE, BGEN
  SvdType svd_t;
  std::string fileU, fileS, fileV;
  std::string filein;
  std::string fileout = "pcaone";
  double memory = 0;  // 0 for disable
  uint nsamples = 0;
  uint nsnps = 0;
  uint k = 10;
  uint maxp = 20;  // maximum number of power iterations
  uint threads = 12;
  uint bands = 64;
  bool pca = true;
  // for normalization
  double scaleFactor = 1;
  // for emu iteration
  uint maxiter = 100;
  double alpha = 0.001;
  // can be tol_emu or tol_pcangsd
  double tolem = 1e-5;
  double tolmaf = 1e-6;
  double maf = 0.0;
  // for arnoldi
  uint ncv = 20;  // max(20, 2*k + 1)
  uint imaxiter = 1000;
  double itol = 1e-6;
  // for halko
  uint oversamples = 10;
  double tol = 1e-4;
  uint buffer = 2;
  uint rand = 1;
  bool perm = false;  // wheather data is permuted
  // for ld stuff
  bool print_r2 = false;
  // for ld pruning
  std::string filebim;  // the 7-th column can be MAF
  int ld_stats = 0;     // 0: adj; 1: std
  double ld_r2 = 0;
  bool ld = false;  // true if tolld > 0
  uint ld_bp = 0;   // base pairs not number of snps
  // for clumping
  std::string clump;
  std::string assoc_colnames;
  double clump_p1 = 0.0001;
  double clump_p2 = 0.01;
  double clump_r2 = 0.5;
  uint clump_bp = 250000;
  // for projection
  int project = 0;
  // for inbreeding
  int inbreed = 0;

  // general
  uint verbose = 1; // verbose level. 0: no verbose; 1: output to screen; 2: debug info 
  uint scale = 0;  // do scaling. 0: just centering. 1: log scaling. 2: cpmed
  bool groff = false;
  bool cpmed = false;
  bool printv = false;
  bool impute = false;    // enable EM-PCA for data with uncertainty. impute information!
  bool noshuffle = false;
  bool emu = false;
  bool pcangsd = false;  // enable pcangsd procedure
  // bool fancyem = false;  // true if emu/pcangsd + svd 2
  bool mev = true;
  bool out_of_core = false;  // otherwise load all matrix into RAM.
  int ploidy = 2;
  int seed = 101; // seeding 
  bool keepsnp = false;
  bool center = true; // false if G is raw data likelihood
  bool estaf = true; // false if project,inbreed enabled
  // bool estpi = false; // true if output pi is needed

  std::ostringstream ss;
};

#endif  // CMD_H_
