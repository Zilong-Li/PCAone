# v0.9.0: full change notes

The detailed notes behind the v0.9.0 entry of the
[change log](../docs/changelog.md): what changed in the code, why, and how it
was measured. The change log itself lists only what users see. v0.8.1 was never
tagged, so this covers every commit since v0.8.0: 4012509, 05d8458, 885e58a,
ede0ee2, 274c429 and fb15458.

## PCA of a BED on the 2-bit calls (fb15458)

- `--svd 1` and `--svd 2` of a BED compute both products of each power iteration from the 2-bit
  calls, in-core and out-of-core, instead of decoding them to doubles. `Data::packed_pca()`
  decides: BED input, a PCA with `--svd 1` or `2`, no EM (`--emu`, `--pcangsd`), centred, with
  `--scale -9` (the default) or `0`, and none of LD, `--evaladmix`, `--project`, `--selection` or
  `--inbreed`. `--no-packed` (advanced) turns it off. The log line `initialize ... on the 2-bit
  calls` says when it is on.
- `PackedGenotypes` (`src/PackedPca.{hpp,cpp}`) holds the records, `(N + 3) / 4` bytes a site,
  and per site `a_j`, `b_j` and `m_j`: an entry of `X` is `a_j + b_j h[c]` for a call of code `c`
  and `m_j` for a missing one. For the RSVD `a_j = -b_j f_j`, `m_j = 0` (mean imputation), and
  `b_j = sqrt(ploidy / (f_j (1 - f_j)))` with the default scaling, 1 with `--scale 0`.
- `X * Y` (`mul_X_Y`) takes the sites 4 at a time with `packed_product_add` of the projection
  (`src/PackedScore.cpp`): the 256 sums of the 4 rows of `Y` for the codes the 4 sites can have,
  each row weighted by the value of its code at its site, are tabulated once, and each sample
  adds one row per 4 sites, read by its byte of codes. Samples in panels of at least 2,048, sites
  in chunks of at least 64, for about 64 tasks.
- `X' * Y` (`mul_Xt_Y`) is the new `packed_tproduct_add`: the samples are taken 4 at a time, one
  byte of codes at a site, and the 256 sums of their rows of `Y` times `h[c]` are tabulated once
  for all the sites, so each site adds one row per 4 samples. The centring is added as `a_j`
  times the column sums of `Y`, and the missing calls as `(m_j - a_j)` times the sum of the rows
  of `Y` of the samples missing at the site, found 32 samples a 64-bit word from the codes
  (`missing_word`). The columns of `Y` go 32 at a time, the widest sums kept in registers.
- Both products are split into about 64 tasks whatever `-n` is, each adding into its own partial
  sums, which are then added in a fixed order. So the result does not depend on `-n`.
- In-core, `FileBed::read_all_packed()` reads the BED in 16 MB chunks with `pread` on every
  thread and copies each record kept straight to its place in the winSVD shuffle, drawn as
  `permute_matrix()` draws it, so the records are held once and the `N x M` doubles are never
  formed. With `--maf`, a first pass finds `F` and so the sites kept. The blocks of a band, which
  share the same test matrix, are multiplied in one call.
- Out-of-core, `RsvdOpData::add_block()` takes the records of the block from
  `FileBed::read_block_packed()` (with the same read-ahead), estimates `F` on the first pass, and
  hands them to the products undecoded.
- `TallQR` (`src/LinAlg.{hpp,cpp}`) replaces Eigen's `HouseholderQR` of `G` and `H` in every RSVD,
  whatever the input. It factors panels of at least 8,192 rows (and of `4 s` for `s` columns, at
  most 1,024 panels) in parallel, then the stacked `R` once more (TSQR). Eigen's QR of a matrix of
  a few columns runs on one thread and took as long as an epoch over a million sites. The panels
  have the same height whatever `-n` is; with one panel it is the `HouseholderQR` of before.
- `Data::save_snps_in_mbim()` reads the `.bim` of a BED at once, finds its lines with `memchr` and
  writes the `.mbim` on all threads with `write_lines`. Reading it a line at a time into strings
  took as long as an epoch of the winSVD. A PGEN's `.pvar` goes the old way.
- Measured on 1000 Genomes, 3,202 x 6.84M: in-core 347 s and 177 GB before, 24 s and 8.2 GB
  after; `-m 16`: 244 s before, 119 s after. The PCs match those of the decoded path and do not
  depend on `-n`.
- On `example/plink` (400 x 656,281, `-k 10`, `-n 2`, cloud x86_64 with 2 cores), against the
  v0.8.0 binary:

  | | v0.8.0 | v0.9.0, `--no-packed` | v0.9.0 |
  |---|---|---|---|
  | winSVD in-core | 9.6 to 10.0 s, 2.29 GB | 6.9 s, 2.19 GB | 3.6 to 3.9 s, 0.35 GB |
  | winSVD `-m 0.05` | 12.7 s, 0.36 GB | 9.1 s, 0.26 GB | 8.5 s, 0.24 GB |
  | sSVD in-core | 18.6 s, 2.29 GB | 13.9 s, 2.19 GB | 7.8 s, 0.40 GB |

  The `.eigvecs`, `.eigvecs2`, `.eigvals`, `.sigvals` and `.mbim` are byte-identical to v0.8.0
  in all three. 4 (winSVD) and 3 (sSVD) of the 6.56M in-core `.loadings` differ by one in the
  last printed digit, the same with `--no-packed`, so from the panel QR. With `-n 1` the output is
  byte-identical to `-n 2`, and with `--maf 0.05` or `-S` the 2-bit and the decoded paths give
  the same `.eigvecs`.
- `tests/test_packed_pca.cpp` (`make test_packed_pca`, run by `make test_prefetch`) checks
  `X' * Y`, `X * Y` and `X * Y` added against the decoded matrix: numbers of samples that do and
  do not fill the last byte, missing calls imputed with 0 and with another value, one and several
  slabs of columns, sub-ranges of sites, the records of a block given apart (out-of-core), and
  several numbers of threads.

## Projection on the 2-bit calls (274c429)

- `--project 1` and `2` on a BED, or on the hard calls of a PGEN, read the blocks as 2-bit calls
  (`Data::packed_calls()`, `read_block_packed()`) and score them with `packed_product_add`, as
  above: 4 sites at a time, one row of the table per sample and 4 sites. The values of the codes
  at each site carry the centring by the reference `.mbim` frequencies and the scaling of the
  reference `.sigvals`; a missing call is 0.
- For `--project 2`, the missing calls of each sample are found from the codes, 32 samples at a
  time, and only they change its normal equations: `add_gram()` adds or removes `v_j v_j'` for
  them, two rows of the packed `K x K` sums in registers.
- PLINK and PGEN targets of `--project 1` and `2` now stream by default, within 2 GB unless `-m`
  is given; `-m 0` reads the target in-core as before. `--project-bootstrap` and PGEN dosages
  decode the blocks.
- 200,000 simulated samples, 98,000 variants with 1% missing calls, 20 PCs (a 4.6 GB BED in the
  page cache): `--project 1` 2.8 s and 0.5 GB, `--project 2` 4.0 s and 1.9 GB on an Apple M3 Pro
  (11 threads), against 12 s and 16 s decoding the calls (`-m 4`) and 11 s for plink2 `--score
  variance-standardize`. On a shared Xeon E5-2699 v4 (16 threads), 3 to 4.4 s and 5 to 6.7 s,
  against 37 s and 42 s decoding and 44 to 72 s for plink2. The scores are the same as when
  decoded. `scripts/benchmark/projection_vs_plink2.sh` reproduces the comparison.

## Out-of-core projection (885e58a)

- `run_projection_blockwise()` (`src/Projection.cpp`) runs `--project 1` and `2` with `-m` on a
  BED or PGEN in one pass over the matched sites, centred by the reference frequencies
  (`use_reference_af()`), holding sums per sample instead of the `9 x N x M` bytes of the in-core
  path.
- `--project 1` adds `G_b V_b S^-1` per block. `--project 2` solves the least squares per sample
  from its normal equations in `V`, `A z = b` with `u = z / S`: `b_i += G_b V_b`, and `A_i`,
  packed to `K(K+1)/2` numbers, is the sum over the block's sites less those the sample misses,
  or the sum over its called sites, whichever touches fewer sites. Each thread owns a chunk of
  samples, so there are no races and the result does not depend on `-n`. A sample with fewer
  than `K` called sites is NA.
- The bootstrap uses the in-core draws (`PortableRng(seed)`, a site-major byte array) and keeps
  `A`, `b` and the number of distinct called sites drawn per replicate; `V` weighted by the draws
  gives `b` of all replicates in one product per block.
- The memory plan splits the samples into groups, one pass each, only when the sums exceed half
  of what `-m` leaves; at least 256 samples a group and 64 sites a block, with a warning. With
  several passes `.proj.bootstrap.eigvecs` is written per pass and merged in replicate-major
  order. `add_X_Y_panels` multiplies row panels of about 8 MB, since Eigen's GEMM packed half of a
  tall block outside `-m`.
- `-m` with `--project 3`, or with BEAGLE input, stops with an error; it was ignored.
- Checked: in-core and `-m` outputs (`.eigvecs`, `.proj.bootstrap.{tsv,cov.tsv,eigvecs}`) agree
  in every digit written, for BED and PGEN, hard calls and dosages, with and without
  `--no-prefetch`, standardized and unscaled references, 0 to 99.9% missing calls, reference
  sites absent, swapped alleles, and 1 to 3 passes. A 60,000 x 20,000 target (10 GB in-core):
  `--project 1 -m 2` 4.9 s, `--project 2 -m 2` 7.8 s, 50 bootstrap replicates 43.5 s in 2 passes,
  peak RSS at most 2.04 GB (cloud x86_64, 2 cores). `tests/test_project_out_of_core.py`
  (`make test_project_out_of_core`).
- **Fix**: `--project-bootstrap` counted the draws, not the distinct sites, towards the `K` sites
  a replicate needs: a sample with 2 called sites each drawn twice passed with `K = 4`, and its
  replicate came from a singular system. It now counts distinct sites, in-core and with `-m`; in
  the test, 401 of 14,000 replicates of near-empty samples became NA.

## Sample IDs in the projection outputs (ede0ee2)

- `--project` (every mode, in-core and `-m`) also writes `.eigvecs2`, as a PCA does: a
  `#FID IID PC1 ...` header and the IDs of `read_sample_ids()`: the `.fam`, the `.psam` (FID `0`
  without a FID column) or the BEAGLE header for `--project 3` (FID = IID). A sample without a
  score is `NA`. The plain `.eigvecs` is unchanged.
- **breaking**: the bootstrap files have `FID` and `IID` after `sample`: `.proj.bootstrap.tsv`
  (`sample FID IID pc ...`), `.proj.bootstrap.cov.tsv` (`sample FID IID pc_x ...`) and
  `.proj.bootstrap.eigvecs` (`replicate sample FID IID PC1 ...`). If the IDs cannot be read, the
  columns are NA, with a warning.

## Text output and the readers of `-P` (4012509, 05d8458)

- One writer (`src/TextIO.{hpp,cpp}`) for every large text output: `put_g6` is printf `%.6g` and
  `put_f6` `%.6f`, with exact rounding (ties to even from the exact binary value), about 35 ns a
  value against 150 ns for `snprintf` and 230 ns through an ostream. `write_lines()` formats lines
  on all threads in about 4 MB per thread, and a `std::thread` writes each block while the next
  is formatted. Same bytes as before, except that a NaN is always `nan` (glibc printed `-nan` for
  some). 1M x 10 `.loadings` on an Apple Silicon Linux VM: 1.70 s before, 0.30 s on one
  thread, 0.14 s on four.
- Converted: `.sigvals`, `.eigvals`, `.eigvecs`, `.eigvecs2` (now written from `U` in memory, not
  re-read from `.eigvecs`), `.loadings`, `.mbim`, the PCAngsd `.cov`, the selection scans, the
  projection `.eigvecs` and bootstrap files, `.hwe`, `.ld.gz`, `.ld.prune.*`, clumping, and the
  evalAdmix matrices and `.kin0`.
- **Fix**: writing `.eigvecs2` from memory fixes an empty (header-only) `.eigvecs2` for a `.psam`
  without a header line.
- The `.loadings`, `.eigvecs` and `.mbim` of `-P/--USV` are read by `LineReader` in blocks, the
  next on a thread while this one is parsed on all threads, and the numbers by `parse_double()`:
  `std::from_chars` with libstdc++, and with libc++ (Apple Clang, which has no `from_chars` for
  double) an exact fast path for up to 19 digits and decimal exponents up to 22, else `strtod`.
  Every value is the double `strtod` gives. 5M x 10 `.loadings` on cloud x86_64: 5.17 s before,
  1.59 s on one thread, 0.85 s on two.
- **Fix**: `read_usv` read a field that is not a number, such as the `NA` of a projected sample,
  as 0. LD and `--evaladmix` now stop with an error naming the row.
- `Utils` is split into `LinAlg`, `Stats`, `TextIO` and `Variants`; dead and duplicated helpers
  are removed (four sample-ID readers become one, three `.mbim` matchers one).
- Checked: `tests/test_textio.cpp` (`make test_textio`) compares the formatters with printf on
  10M values, `write_matrix` with Eigen's output at 1 to 8 threads, and `parse_double` with
  `strtod` bit for bit on 6.2M strings. 97 output files of 25 runs (every PCA method, `-m`, PGEN,
  BEAGLE, `--maf`, selection, inbreeding, evalAdmix, LD, projection and bootstrap) were
  byte-identical to the v0.8.0 binary.
