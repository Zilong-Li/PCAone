#ifndef CMD_H_
#define CMD_H_

#include <iostream>
#include <iterator>
#include <sstream>

typedef unsigned int uint;
typedef unsigned long long uint64;

enum class FileType { PLINK, CSV, BEAGLE, BINARY, BGEN };

enum class SvdType { IRAM, PCAoneAlg1, PCAoneAlg2, FULL };

class Param {
 public:
  Param(int argc, char** argv);
  ~Param();

  FileType file_t;  // PLINK, CSV, BEAGLE, BGEN
  SvdType svd_t;
  std::string filein;
  std::string fileout = "pcaone";
  double memory = 0;  // 0 for disable
  uint nsamples = 0;
  uint nsnps = 0;
  uint k = 10;
  uint maxp = 40;  // maximum number of power iterations
  uint threads = 10;
  uint blocksize = 0;
  uint bands = 64;
  // for emu iteration
  uint maxiter = 100;
  double alpha = 0.001;
  // can be tol_emu or tol_pcangsd
  double tolem = 1e-6;
  double tolmaf = 1e-4;
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
  // for ld pruning
  std::string filebim;
  int ld_stats = 0;  // 0: adj; 1: std
  double ld_r2 = 0;
  bool ld = false;       // true if tolld > 0
  uint ld_bp = 1000000;  // base pairs not number of snps
  // for clumping
  std::string clump;
  std::string assoc_colnames;
  double clump_p1 = 0.0001;
  double clump_p2 = 0.01;
  double clump_r2 = 0.5;
  uint clump_bp = 250000;

  // general
  uint scale = 0;  // do scaling. 0: just centering. 1: log scaling. 2: cpmed
  bool groff = false;
  bool cpmed = false;
  bool printv = false;
  bool printu = false;
  bool runem = false;
  bool noshuffle = false;
  bool emu = false;
  bool pcangsd = false;  // read GP field for PCAngsd instead of GT.
  bool verbose = false;
  bool mev = true;
  bool out_of_core = false;  // otherwise load all matrix into RAM.
  bool haploid = false;
  bool diploid = false;
  bool keepsnp = false;

  std::ostringstream ss;
};

#endif  // CMD_H_
