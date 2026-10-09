#include "FileBinary.hpp"

using namespace std;

void FileBin::check_file_offset_first_var() {
  setlocale(LC_ALL, "C");
  ios_base::sync_with_stdio(false);
  // magic += missing_points.size() * sizeof(uint64);
  long long offset = ibyte * 2 + nsnps * bytes_per_snp;
  if (ifs_bin.tellg() == offset) {
    // reach the end of bed, reset the position to the first variant;
    ifs_bin.seekg(ibyte * 2, std::ios_base::beg);
  } else if (ifs_bin.tellg() == ibyte * 2) {
    ;
  } else {
    ifs_bin.seekg(ibyte * 2, std::ios_base::beg);
    if (params.verbose) cao.warn("confirm you are running the window-based RSVD (algorithm2)");
  }
}

void FileBin::read_all() {
  check_file_offset_first_var();
  G = Mat2D::Zero(nsamples, nsnps);
  Eigen::VectorXf fg(nsamples);
  for (Eigen::Index i = 0; i < G.cols(); i++) {
    ifs_bin.read((char*)fg.data(), bytes_per_snp);
    G.col(i) = fg.cast<double>();
    // --scale 1-4 standardized the columns before they were written, so this
    // only removes the float rounding of the mean. --scale 0 (and the default,
    // which is no transform for CSV) must stay as it is: centring here made the
    // -m run a different, centred PCA from the in-core one.
    if (params.scale >= 1) G.col(i).array() -= G.col(i).mean();
  }
}

// TODO : can standardize
void FileBin::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize) {
  // magic += missing_points.size() * sizeof(uint64);
  // check where we are
  long long offset = ibyte * 2 + start_idx * bytes_per_snp;
  if (ifs_bin.tellg() != offset) cao.error("something wrong with read_snp_block!\n");
  uint actual_block_size = stop_idx - start_idx + 1;
  G = Mat2D(nsamples, actual_block_size);
  Eigen::VectorXf fg(nsamples);
  for (Eigen::Index i = 0; i < G.cols(); i++) {
    ifs_bin.read((char*)fg.data(), bytes_per_snp);
    G.col(i) = fg.cast<double>();
    // --scale 1-4 standardized the columns before they were written, so this
    // only removes the float rounding of the mean. --scale 0 (and the default,
    // which is no transform for CSV) must stay as it is: centring here made the
    // -m run a different, centred PCA from the in-core one.
    if (params.scale >= 1) G.col(i).array() -= G.col(i).mean();
  }
  // While this block is used, the kernel reads the next one into the page
  // cache (after the last block, the first one again). No buffer is held, so
  // the memory is that of -m.
  if (start_idx == 0) first_block_size = actual_block_size;
  if (params.noprefetch) return;
  const uint64 next = stop_idx + 1 < nsnps ? stop_idx + 1 : 0;
  const uint64 count = next ? std::min<uint64>(actual_block_size, nsnps - next) : first_block_size;
  try {
    if (!readahead) readahead = std::make_unique<PCAone::ReadAhead>(params.filein);
    readahead->hint(ibyte * 2 + next * bytes_per_snp, count * bytes_per_snp);
  } catch (const std::exception&) {
    // only a hint
  }
}
