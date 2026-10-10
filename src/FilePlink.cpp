/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/FilePlink.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "FilePlink.hpp"
#include "BedShuffle.hpp"
#include <climits>
#include <fstream>
#include <sys/stat.h>
#ifdef __linux__
#include <sys/sysmacros.h>
#endif
#include <filesystem>

using namespace std;

void FileBed::check_file_offset_first_var() {
  setlocale(LC_ALL, "C");
  ios_base::sync_with_stdio(false);
  long long offset = 3 + nsnps * bed_bytes_per_snp;
  if (bed_ifstream.tellg() == offset) {
    // reach the end of bed, reset the position to the first variant;
    bed_ifstream.seekg(3, std::ios_base::beg);
  } else if (bed_ifstream.tellg() == 3) {
    ;
  } else {
    bed_ifstream.seekg(3, std::ios_base::beg);
    if (params.verbose) cao.warn("confirm you are running the window-based RSVD (algorithm2)");
  }
}

FileBed::~FileBed() {
  if (prefetcher) {
    prefetcher->cancel();
    if (params.verbose > 1)
      cao.print(tick.date(), "BED blocks read in the background:", prefetcher->predicted(), ", waited for reads",
                prefetcher->wait_seconds(), "seconds");
  }
}

const uchar* FileBed::read_records(uint64 start_idx, uint64 count) {
  if (!params.noprefetch || !bed_order.empty()) {
    try {
      if (!prefetcher) {
        prefetcher = std::make_unique<PCAone::RecordPrefetcher>(bed_path, 3, bed_bytes_per_snp, nsnps);
        if (!bed_order.empty()) prefetcher->set_order(bed_order);
        if (params.noprefetch) prefetcher->set_read_ahead(false);
      }
      return prefetcher->get(start_idx, count);
    } catch (const std::exception& e) {
      cao.error(std::string("cannot read ") + bed_path + ": " + e.what());
    }
  }
  if (bed_ifstream.tellg() != (long long)(3 + start_idx * bed_bytes_per_snp)) cao.error("read_records: offset wrong!");
  if (inbed.size() < bed_bytes_per_snp * count) inbed.resize(bed_bytes_per_snp * count);
  bed_ifstream.read(reinterpret_cast<char*>(inbed.data()), bed_bytes_per_snp * count);
  if (!bed_ifstream) cao.error("cannot read " + bed_path);
  return inbed.data();
}

void FileBed::read_all() {
  check_file_offset_first_var();
  // Begin to decode the plink bed
  inbed.resize(bed_bytes_per_snp * nsnps);
  bed_ifstream.read(reinterpret_cast<char*>(inbed.data()), bed_bytes_per_snp * nsnps);
  uint64 c, i, j, b, k;
  uchar buf;
  // sometimes no need to estimate AF, e.g. projection
  if (params.dopca && !frequency_was_estimated) {
    F = Mat1D::Zero(nsnps);
    uint64 nmono = 0;
    // estimate allele frequency first
#pragma omp parallel for private(i, j, b, c, k, buf) reduction(+ : nmono)
    for (i = 0; i < nsnps; ++i) {
      for (b = 0, c = 0, j = 0; b < bed_bytes_per_snp; ++b) {
        buf = inbed[i * bed_bytes_per_snp + b];
        for (k = 0; k < 4; ++k, ++j) {
          if (j < nsamples) {
            if (BED2GENO[buf & 3] != BED_MISSING_VALUE) {
              F(i) += BED2GENO[buf & 3];
              c++;
            }
            buf >>= 2;
          }
        }
      }
      if (c == 0)
        F(i) = 0;
      else
        F(i) /= c;
      // should remove sites with F=0 and 1.0. counted, not logged here: the
      // logger is not thread-safe, and one line per site garbled the log
      if (F(i) == 0.0 || F(i) == 1.0) ++nmono;
    }
    warn_monomorphic(nmono);
    filter_snps_resize_F();  // filter and resize nsnps
  }

  const bool filter = !keepSNPs.empty();
  if (filter) nsnps = keepSNPs.size();
  G = Mat2D::Zero(nsamples, nsnps);  // fill in G with new size after filtering

  if (params.missme) C = ArrBool::Zero((uint64)nsnps * nsamples);
  if (params.pcangsd) P = Mat2D::Zero(nsamples * 2, nsnps);

#pragma omp parallel for private(i, j, b, c, k, buf)
  for (i = 0; i < nsnps; ++i) {
    uint s = filter ? keepSNPs[i] : i;
    for (b = 0, c = 0, j = 0; b < bed_bytes_per_snp; ++b) {
      buf = inbed[s * bed_bytes_per_snp + b];
      for (k = 0; k < 4; ++k, ++j) {
        if (j < nsamples) {
          G(j, i) = BED2GENO[buf & 3];
          buf >>= 2;
          if (params.missme && (G(j, i) == BED_MISSING_VALUE)) C[i * nsamples + j] = 1;
        }
      }
    }
    if (params.pcangsd) {
      for (j = 0; j < nsamples; ++j) {
        if (G(j, i) == BED2GENO[3]) {  // RR
          P(2 * j + 0, i) = 1.00;
          P(2 * j + 1, i) = 0.00;
        } else if (G(j, i) == BED2GENO[0]) {  // AA
          P(2 * j + 0, i) = 0.00;
          P(2 * j + 1, i) = 0.00;
        } else if (G(j, i) == BED2GENO[2]) {  // RA
          P(2 * j + 0, i) = 0.00;
          P(2 * j + 1, i) = 1.00;
        } else {
          P(2 * j + 0, i) = 0.333333;
          P(2 * j + 1, i) = 0.333333;
        }
      }
    }
    // do centering and initialing
    if (params.center) {
      for (j = 0; j < nsamples; ++j) {
        if (G(j, i) == BED_MISSING_VALUE)
          G(j, i) = 0.0;  // impute to mean
        else
          G(j, i) -= F(i);
      }
    }
  }

  if (params.missme) {
    p_miss = (double)C.count() / (double)C.size();
    cao.print(tick.date(), "the proportion of missingness  is", p_miss);
  }

  // read bed to matrix G done and close bed_ifstream
  inbed.clear();
  inbed.shrink_to_fit();
}

void FileBed::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize) {
  uint actual_block_size = stop_idx - start_idx + 1;
  // if G is not initial then initial it
  // if actual_block_size is smaller than blocksize, don't resize G;
  if (G.cols() < blocksize || (actual_block_size < blocksize)) {
    G = Mat2D::Zero(nsamples, actual_block_size);
  }
  uint64 c, b, i, j, k, snp_idx;
  uchar buf;
  const uchar* bed = read_records(start_idx, actual_block_size);
  if (!params.dopca) frequency_was_estimated = true;  // read AF from external
  if (frequency_was_estimated) {
#pragma omp parallel for private(i, j, b, k, snp_idx, buf)
    for (i = 0; i < actual_block_size; ++i) {
      snp_idx = start_idx + i;
      const double f = F(snp_idx);
      double scale_factor = 1.0;
      if (standardize && params.scale == SCALE_STANDARDIZE_GENETIC) {
        double sd = sqrt(f * (1.0 - f));
        if (sd > VAR_TOL) scale_factor = sqrt((double)params.ploidy) / sd;
      }
      for (b = 0, j = 0; b < bed_bytes_per_snp; ++b) {
        buf = bed[i * bed_bytes_per_snp + b];
        for (k = 0; k < 4; ++k, ++j) {
          if (j < nsamples) {
            if (params.center) {
              G(j, i) = centered_geno_lookup(buf & 3, snp_idx);
            } else {
              G(j, i) = BED2GENO[buf & 3];
            }
            G(j, i) *= scale_factor;
            buf >>= 2;
          }
        }
      }
    }
  } else {
    // estimate allele frequencies
    if (start_idx == 0) nmono_seen = 0;  // a pass over the blocks restarted before F was complete
    uint64 nmono = 0;
#pragma omp parallel for private(c, i, j, b, k, snp_idx, buf) reduction(+ : nmono)
    for (i = 0; i < actual_block_size; ++i) {
      snp_idx = start_idx + i;
      c = 0;
      // start from 0: the LD code may read a block again before the pass that
      // estimates F has reached the last block, e.g. clumping several files
      F(snp_idx) = 0.0;
      for (b = 0, j = 0; b < bed_bytes_per_snp; ++b) {
        buf = bed[i * bed_bytes_per_snp + b];
        for (k = 0; k < 4; ++k, ++j) {
          if (j < nsamples) {
            if ((buf & 3) != 1) {
              // g is {0, 0.5, 1}
              F(snp_idx) += BED2GENO[buf & 3];
              c++;
            }
            buf >>= 2;
          }
        }
      }
      // calculate F and centered_geno_lookup
      if (c == 0) {
        F(snp_idx) = 0;
      } else {
        F(snp_idx) /= c;
      }
      // should remove sites with F=0 and 1.0
      if (F(snp_idx) == 0.0 || F(snp_idx) == 1.0) ++nmono;
      // do centering and initialing
      centered_geno_lookup(1, snp_idx) = 0.0;                       // missing
      centered_geno_lookup(0, snp_idx) = BED2GENO[0] - F(snp_idx);  // minor hom
      centered_geno_lookup(2, snp_idx) = BED2GENO[2] - F(snp_idx);  // het
      centered_geno_lookup(3, snp_idx) = BED2GENO[3] - F(snp_idx);  // major hom
      // get centered and standardized G
      const double f = F(snp_idx);
      double scale_factor = 1.0;
      if (standardize && params.scale == SCALE_STANDARDIZE_GENETIC) {
        double sd = sqrt(f * (1.0 - f));
        if (sd > VAR_TOL) scale_factor = sqrt((double)params.ploidy) / sd;
      }
      for (b = 0, j = 0; b < bed_bytes_per_snp; ++b) {
        buf = bed[i * bed_bytes_per_snp + b];
        for (k = 0; k < 4; ++k, ++j) {
          if (j < nsamples) {
            G(j, i) = centered_geno_lookup(buf & 3, snp_idx) * scale_factor;
            buf >>= 2;
          }
        }
      }
    }
    nmono_seen += nmono;
  }

  if (stop_idx + 1 == nsnps && !frequency_was_estimated) {
    frequency_was_estimated = true;
    warn_monomorphic(nmono_seen);
  }
}

void FileBed::read_block_update(
    uint64 start_idx, uint64 stop_idx, const Mat2D& U, const Mat1D& svals, const Mat2D& VT, bool standardize) {
  uint actual_block_size = stop_idx - start_idx + 1;
  if (G.cols() < blocksize || (actual_block_size < blocksize)) {
    G = Mat2D::Zero(nsamples, actual_block_size);
  }
  const uchar* bed = read_records(start_idx, actual_block_size);
  uint64 b, i, j, snp_idx;
  uint ks = svals.rows();
  uint ki, k;
  uchar buf;
#pragma omp parallel for private(i, j, b, ki, k, snp_idx, buf)
  for (i = 0; i < actual_block_size; ++i) {
    snp_idx = start_idx + i;
    const double f = F(snp_idx);
    double scale_factor = 1.0;
    if (standardize && params.scale == SCALE_STANDARDIZE_GENETIC) {
      double sd = sqrt(f * (1.0 - f));
      if (sd > VAR_TOL) scale_factor = sqrt((double)params.ploidy) / sd;
    }
    for (b = 0, j = 0; b < bed_bytes_per_snp; ++b) {
      buf = bed[i * bed_bytes_per_snp + b];
      for (ki = 0; ki < 4; ++ki, ++j) {
        if (j < nsamples) {
          G(j, i) = centered_geno_lookup(buf & 3, snp_idx);
          // EMU procedure
          if (params.emu && ((buf & 3) == 1)) {
            G(j, i) = 0.0;
            for (k = 0; k < ks; ++k) {
              G(j, i) += U(j, k) * svals(k) * VT(k, snp_idx);
            }
            // map to domain(0,1)
            G(j, i) = fmin(fmax(G(j, i), -F(snp_idx)), 1 - F(snp_idx));
          }
          if (params.pcangsd) {
            double pt, p1, p2, p0;
            pt = 0.0;
            for (k = 0; k < ks; ++k) {
              pt += U(j, k) * svals(k) * VT(k, snp_idx);
            }
            pt = (pt + 2.0 * F(snp_idx)) / 2.0;
            pt = fmin(fmax(pt, 1e-4), 1.0 - 1e-4);
            if ((buf & 3) == 3) {
              p0 = 1.00;
              p1 = 0.00;
              p2 = 0.00;
            } else if ((buf & 3) == 0) {
              p0 = 0.00;
              p1 = 0.00;
              p2 = 1.00;
            } else if ((buf & 3) == 2) {
              p0 = 0.00;
              p1 = 1.00;
              p2 = 0.00;
            } else {
              p0 = 0.333333;
              p1 = 0.333333;
              p2 = 0.333333;
            }
            p0 *= (1.0 - pt) * (1.0 - pt);
            p1 *= 2.0 * pt * (1.0 - pt);
            p2 *= pt * pt;
            G(j, i) = (p1 + 2.0 * p2) / (p0 + p1 + p2) - 2.0 * F(snp_idx);
          }
          // scale G
          G(j, i) *= scale_factor;
          // shift packed data and throw away genotype just processed.
          buf >>= 2;
        }
      }
    }
  }
}

bool FileBed::apply_permutation(Param& config) {
  // winSVD changes Omg only between -w bands of bandFactor blocks, so the
  // order within a band is free (BedShuffle.hpp). One bucket per band keeps
  // the writes large; one per block would shrink them as the file grows.
  const uint64 bucket = (uint64)blocksize * bandFactor;
  if (!config.bedcopy) {
    // The order of the copy, read from the input: each block's SNPs are a
    // sorted subset of it, read in file order (RecordPrefetcher::set_order).
    if (nsnps > INT_MAX) cao.error("too many SNPs for the BED permutation.");
    struct stat st;
    if (::stat(bed_path.c_str(), &st) != 0 || (uint64)st.st_size != 3 + (uint64)nsnps * bed_bytes_per_snp)
      cao.error("BED size does not match BIM/FAM dimensions.");
    bed_order = bed_bucket_order(nsnps, bucket, config.seed);
    Eigen::VectorXi indices(nsnps);
    for (uint64 d = 0; d < nsnps; ++d) indices(d) = bed_order[d];
    perm = PermMat(indices);
    prefetcher.reset();
    cao.print(tick.date(), "shuffle SNPs into random -w bands, read from the input; seed:", config.seed,
              ", SNPs per band:", bucket);
    if (on_rotating_disk(bed_path) == 1)
      cao.warn("the BED is on a spinning disk, where reading its shuffled SNPs takes a seek per SNP once the file is "
               "not in the memory. if the BED (" + std::to_string(st.st_size >> 20) + " MiB) is larger than the free "
               "memory, --bed-copy, which writes a shuffled copy first, is usually much faster");
    return false;
  }
  perm = permute_plink(config.filein, config.fileout, config.buffer, bucket, config.seed);
  bed_ifstream.close();
  bed_ifstream.clear();
  bed_ifstream.open(config.filein + ".bed", std::ios::binary);
  bed_ifstream.seekg(3);
  bed_path = config.filein + ".bed";
  prefetcher.reset();
  if (!bed_ifstream) throw std::runtime_error("Cannot reopen permuted BED file");
  return true;
}

std::vector<uint32_t> bed_bucket_order(uint64 nsnps, uint64 bucket, int seed) {
  // the shuffle of permute_plink()
  std::vector<uint32_t> order(nsnps);
  std::iota(order.begin(), order.end(), 0);
  PortableRng rng(seed);
  portable_shuffle(order.begin(), order.end(), rng);
  // rewrite_bed_buckets() keeps the source order within a bucket
  for (uint64 b = 0; b < nsnps; b += bucket) std::sort(order.begin() + b, order.begin() + std::min(nsnps, b + bucket));
  return order;
}

int on_rotating_disk(const std::string& path) {
#ifdef __linux__
  struct stat st;
  if (::stat(path.c_str(), &st) != 0 || major(st.st_dev) == 0) return -1;
  // a whole disk has queue/, a partition takes its disk's
  const std::string dev =
      "/sys/dev/block/" + std::to_string(major(st.st_dev)) + ":" + std::to_string(minor(st.st_dev));
  for (const char* queue : {"/queue/rotational", "/../queue/rotational"}) {
    std::ifstream f(dev + queue);
    int r;
    if (f >> r) return r;
  }
#else
  (void)path;
#endif
  return -1;
}

namespace {
// true if `path` exists and is the same file as `input`, through any spelling,
// symlink or hard link. "-b ./x.perm -o x" truncated the input before.
bool same_file(const std::string& input, const std::string& path) {
  std::error_code ec;
  return std::filesystem::exists(path, ec) && std::filesystem::equivalent(input, path, ec);
}

// Removes the permuted files this run created unless it finished writing them.
struct PermOutputGuard {
  std::vector<std::string> created;
  bool done = false;
  ~PermOutputGuard() {
    if (done) return;
    std::error_code ec;
    for (const auto& path : created) std::filesystem::remove(path, ec);
  }
};
}  // namespace

PermMat permute_plink(std::string& fin, const std::string& fout, uint gb, uint64 bucket, int seed) {
  const uint64 nsnps = count_lines(fin + ".bim");
  const uint64 nsamples = count_lines(fin + ".fam");
  if (!nsnps || !nsamples || nsnps > INT_MAX || !gb || !bucket)
    cao.error("Invalid dimensions or buffer for BED permutation.");
  const uint64 width = (nsamples + 3) / 4;
  const uint64 budget = uint64(gb) * 1073741824ULL;
  if (budget / 2 < width) cao.error("--buffer must hold two SNP records for the BED permutation.");
  const uint64 bytes = nsnps * width + 3;
  std::ifstream in(fin + ".bed", std::ios::binary | std::ios::ate);
  if (!in || in.tellg() != static_cast<std::streamoff>(bytes))
    cao.error("BED size does not match BIM/FAM dimensions.");
  in.seekg(0);
  char header[3];
  in.read(header, 3);
  if (!in || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01)
    cao.error("Incorrect magic number in plink bed file.");
  const std::string prefix = fout + ".perm";
  for (const char* suffix : {".bed", ".bim", ".fam"})
    if (same_file(fin + suffix, prefix + suffix))
      cao.error("the permuted " + prefix + suffix + " would overwrite the input. please use another -o.");

  // the same generator as the in-core, PGEN, BGEN and CSV shuffles: std::shuffle
  // gave a different order for the same --seed on libc++ (macOS) and libstdc++
  std::vector<uint32_t> order(nsnps);
  std::iota(order.begin(), order.end(), 0);
  PortableRng rng(seed);
  portable_shuffle(order.begin(), order.end(), rng);
  cao.print(tick.date(), "shuffle SNPs into random -w bands; seed:", seed, ", SNPs per band:", bucket);

  PermOutputGuard guard;
  std::ofstream out(prefix + ".bed", std::ios::binary);
  if (!out) throw std::runtime_error("Cannot open " + prefix + ".bed");
  guard.created.push_back(prefix + ".bed");
  out.write(header, 3);
  PCAone::rewrite_bed_buckets(in, out, order, width, budget, bucket);
  out.close();
  if (!out) throw std::runtime_error("Cannot finish permuted BED file");

  std::ifstream bim(fin + ".bim");
  std::vector<std::string> lines(std::istream_iterator<Line>{bim}, std::istream_iterator<Line>{});
  if (lines.size() != nsnps) throw std::runtime_error("Cannot read BIM during permutation");
  std::ofstream out_bim(prefix + ".bim");
  if (!out_bim) throw std::runtime_error("Cannot open " + prefix + ".bim");
  guard.created.push_back(prefix + ".bim");
  Eigen::VectorXi indices(nsnps);
  for (uint64 d = 0; d < nsnps; ++d) {
    indices(d) = order[d];
    out_bim << lines[order[d]] << "\n";
  }
  out_bim.close();
  std::ifstream fam(fin + ".fam", std::ios::binary);
  std::ofstream out_fam(prefix + ".fam", std::ios::binary);
  if (!out_fam) throw std::runtime_error("Cannot open " + prefix + ".fam");
  guard.created.push_back(prefix + ".fam");
  out_fam << fam.rdbuf();
  out_fam.close();
  if (!out_bim || !fam || !out_fam) throw std::runtime_error("Cannot write permuted BED metadata");
  guard.done = true;
  fin = prefix;
  return PermMat(indices);
}
