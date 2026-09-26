#!/usr/bin/env Rscript
# Write copies of a PLINK fileset with genotypes set to missing.
#
#   Rscript make_missing.R <in prefix> <out prefix> <design> <seed>
#
# designs: mcar<rate> (e.g. mcar0.2: every call missing with that probability),
# varying (each sample's rate drawn from U(0, 0.4)), batch (samples split at
# random into two batches, each missing its own random 25% of the sites; the
# batch of each sample is written to <out>.batch).
args <- commandArgs(TRUE)
stopifnot(length(args) == 4)
src <- args[1]; out <- args[2]; design <- args[3]; set.seed(as.integer(args[4]))
n <- nrow(read.table(paste0(src, ".fam"))); m <- nrow(read.table(paste0(src, ".bim")))
nb <- ceiling(n / 4)
bed <- readBin(paste0(src, ".bed"), "raw", n = 3 + nb * m)
stopifnot(bed[1] == as.raw(0x6c), bed[2] == as.raw(0x1b), bed[3] == as.raw(0x01))
x <- matrix(as.integer(bed[-(1:3)]), nb, m)              # one column per site
slot <- lapply(0:3, function(k) (x %/% 4^k) %% 4)        # 2-bit code of sample 4*(row-1)+k
ind <- lapply(0:3, function(k) outer(4 * (seq_len(nb) - 1) + k + 1, rep(1, m)))  # 1-based sample
if (grepl("^mcar", design)) {
  rate <- as.numeric(sub("mcar", "", design))
  miss <- function(i) matrix(runif(nb * m) < rate, nb, m)
} else if (design == "varying") {
  r <- runif(n, 0, 0.4)
  miss <- function(i) matrix(runif(nb * m), nb, m) < matrix(r[pmin(i, n)], nb, m)
} else if (design == "batch") {
  b <- sample(0:1, n, TRUE); SA <- runif(m) < 0.25; SB <- runif(m) < 0.25
  write.table(b, paste0(out, ".batch"), row.names = FALSE, col.names = FALSE)
  miss <- function(i) { bi <- matrix(b[pmin(i, n)], nb, m); (bi == 0 & rep(SA, each = nb)) | (bi == 1 & rep(SB, each = nb)) }
} else stop("unknown design ", design)
y <- 0L; nmiss <- 0
for (k in 1:4) {
  s <- slot[[k]]; real <- ind[[k]] <= n                  # never touch the padding
  s[miss(ind[[k]]) & real] <- 1L                         # 01 = missing
  nmiss <- nmiss + sum(s[real] == 1L)
  y <- y + s * 4^(k - 1)
}
con <- file(paste0(out, ".bed"), "wb")
writeBin(c(as.raw(c(0x6c, 0x1b, 0x01)), as.raw(as.vector(y))), con); close(con)
invisible(file.copy(paste0(src, c(".bim", ".fam")), paste0(out, c(".bim", ".fam")), overwrite = TRUE))
cat(sprintf("%s: %.4f of the calls missing\n", out, nmiss / (n * m)))
