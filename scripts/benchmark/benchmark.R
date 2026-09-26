#!/usr/bin/env Rscript
# Reproduce every table in docs/evaladmix.md from the outputs of run_all.sh.
#
#   Rscript benchmark.R <workdir> <path/to/relateAdmix/data>
#
# GENESIS is optional: if installed, PC-Relate rows are included.
#   BiocManager::install("GENESIS")

args <- commandArgs(TRUE)
WORK <- if (length(args) >= 1) args[1] else stop("usage: benchmark.R <workdir> <datadir>")
DATA <- if (length(args) >= 2) args[2] else file.path(WORK, "src/relateAdmix/data")
OUT  <- file.path(WORK, "out")
SRC  <- file.path(WORK, "src")

## ---------------------------------------------------------------- helpers ---
read_bed <- function(prefix) {                    # -> n x m matrix of {0,1,2}
  fam <- read.table(paste0(prefix, ".fam")); n <- nrow(fam)
  bim <- read.table(paste0(prefix, ".bim")); m <- nrow(bim)
  bed <- readBin(paste0(prefix, ".bed"), "raw", n = 3 + ceiling(n/4) * m)
  stopifnot(bed[1] == as.raw(0x6c), bed[2] == as.raw(0x1b), bed[3] == as.raw(0x01))
  raw <- matrix(bed[-(1:3)], nrow = ceiling(n/4))
  lut <- c(2L, NA_integer_, 1L, 0L)               # 00=hom A1, 01=miss, 10=het, 11=hom A2
  G <- matrix(NA_integer_, n, m)
  bits <- function(v) { x <- as.integer(v); cbind(x %% 4, (x %/% 4) %% 4, (x %/% 16) %% 4, (x %/% 64) %% 4) }
  for (j in seq_len(m)) G[, j] <- lut[as.vector(t(bits(raw[, j]))) + 1L][1:n]
  dimnames(G) <- list(fam$V2, bim$V2)
  G
}
rmse <- function(x, y) sqrt(mean((x - y)^2))

## ------------------------------------------------------------------ truth ---
# The pedigree is visible in the data: individuals are ordered so that relatives
# are consecutive even-odd pairs, ten pairs per class. "Truth" below is the
# THEORETICAL kinship for the relationship, never any method's output.
G <- read_bed(file.path(DATA, "smallPlink"))
n <- nrow(G)
pairs <- t(combn(n, 2))                                   # 1-based, upper triangle
cls <- rep("unrelated", nrow(pairs))
designed <- (pairs[,2] == pairs[,1] + 1) & (pairs[,1] %% 2 == 1)   # 1-based even-odd
grp <- function(i0) {                                     # i0 = 0-based index
  if (i0 <= 5) "unrelated" else if (i0 <= 25) "duplicate" else if (i0 <= 45) "parent-offspring"
  else if (i0 <= 65) "full sib" else if (i0 <= 85) "half sib"
  else if (i0 <= 105) "first cousin" else "second cousin"
}
cls[designed] <- vapply(pairs[designed, 1] - 1L, grp, "")
theta <- c(unrelated = 0, duplicate = 0.5, "parent-offspring" = 0.25, "full sib" = 0.25,
           "half sib" = 0.125, "first cousin" = 0.0625, "second cousin" = 0.015625)
truth <- theta[cls]
CLS <- names(theta)[-1]
cat(sprintf("pairs: %d (%d related, %d unrelated)\n", nrow(pairs), sum(cls != "unrelated"), sum(cls == "unrelated")))

## ------------------------------------------------------------- estimates ----
est <- list()
gp  <- function(M) M[pairs]                               # matrix -> pair vector

k <- read.table(file.path(OUT, "relateadmix.k"), header = TRUE)
ra <- matrix(0, n, n); ra[cbind(k$ind1 + 1, k$ind2 + 1)] <- k$k1/4 + k$k2/2
est[["RelateAdmix"]] <- gp(ra + t(ra))

# evalAdmix returns the correlation of residuals, which estimates 2*phi
est[["evalAdmix (EM)"]]   <- gp(as.matrix(read.table(file.path(OUT, "evaladmix_em.corres"))))   / 2
est[["evalAdmix (proj)"]] <- gp(as.matrix(read.table(file.path(OUT, "evaladmix_proj.corres")))) / 2

pc1 <- as.matrix(read.table(file.path(OUT, "pcaone.kinship"), header = TRUE, check.names = FALSE))
est[["PCAone --evaladmix"]] <- gp(pc1)

## PCA + projection, R reference implementation
if (file.exists(file.path(SRC, "evalPopStructure/R/evalPCA.R"))) {
  source(file.path(SRC, "evalPopStructure/R/evalPCA.R"))
  p <- colMeans(G)/2; keep <- p > 0.05 & p < 0.95
  g <- t(G[, keep])                                       # SNPs x individuals
  pcaR <- makePCA(g, method = "standard", center = TRUE, scale = FALSE)
  est[["PCA + projection (R)"]] <- gp(evalPCA(pcaR, k = 1)$corres) / 2
}

## PC-Relate, if GENESIS is available
have_genesis <- requireNamespace("GENESIS", quietly = TRUE) &&
                requireNamespace("SNPRelate", quietly = TRUE) &&
                requireNamespace("GWASTools", quietly = TRUE)
pcrelate_npc <- NULL
if (have_genesis) {
  suppressMessages({ library(GENESIS); library(SNPRelate); library(GWASTools); library(gdsfmt) })
  gds <- file.path(OUT, "smallPlink.gds")
  if (!file.exists(gds))
    snpgdsBED2GDS(file.path(DATA, "smallPlink.bed"), file.path(DATA, "smallPlink.fam"),
                  file.path(DATA, "smallPlink.bim"), gds, verbose = FALSE)
  ff <- snpgdsOpen(gds)
  pca <- snpgdsPCA(ff, num.thread = 8, eigen.cnt = 8, maf = 0.05, missing.rate = 0.05, verbose = FALSE)
  snpgdsClose(ff)
  pcs <- pca$eigenvect; rownames(pcs) <- pca$sample.id
  run_pcrelate <- function(npc, correct = TRUE) {
    gd <- GenotypeData(GdsGenotypeReader(gds))
    r  <- pcrelate(GenotypeBlockIterator(gd, snpBlock = 20000), pcs = pcs[, 1:npc, drop = FALSE],
                   ibd.probs = TRUE, maf.thresh = 0.05, small.samp.correct = correct, verbose = FALSE)
    close(gd)
    kb <- as.data.frame(r$kinBtwn)
    i <- match(kb$ID1, rownames(G)); j <- match(kb$ID2, rownames(G))
    M <- matrix(0, n, n); M[cbind(pmin(i,j), pmax(i,j))] <- kb$kin
    list(kin = gp(M + t(M)), k0 = kb$k0, k2 = kb$k2, ID1 = kb$ID1, ID2 = kb$ID2)
  }
  est[["PC-Relate"]] <- run_pcrelate(1)$kin
  est[["PC-Relate (corr off)"]] <- run_pcrelate(1, correct = FALSE)$kin
  pcrelate_npc <- lapply(1:4, function(d) run_pcrelate(d)$kin)
} else {
  message("GENESIS not installed -- PC-Relate rows skipped. BiocManager::install(\"GENESIS\")")
}

## --------------------------------------------------------------- tables -----
u <- cls == "unrelated"
cat("\n=== mean kinship by relationship ===\n")
tb <- sapply(est, function(v) tapply(v, factor(cls, levels = c("unrelated", CLS)), mean))
print(round(cbind(truth = c(0, theta[CLS]), tb), 4))

cat("\n=== RMSE decomposition (the key table) ===\n")
dec <- t(sapply(est, function(v) c(
  `RMSE all` = rmse(v, truth), `RMSE unrelated` = rmse(v[u], truth[u]),
  `RMSE related` = rmse(v[!u], truth[!u]), `bias related` = mean(v[!u] - truth[!u]))))
print(round(dec[order(dec[, "RMSE related"]), ], 5))
cat("\n99.24% of pairs are unrelated, so 'RMSE all' mostly measures the unrelated\n",
    "column. Note the ranking differs between 'RMSE all' and 'RMSE related'.\n", sep = "")

cat("\n=== attenuation (estimate / truth) ===\n")
att <- t(sapply(est, function(v) sapply(CLS, function(c) unname(mean(v[cls == c]) / theta[c]))))
colnames(att) <- CLS
print(round(att, 3))

cat("\n=== effect of the [-1,1] clip ===\n")
clipped <- t(sapply(est, function(v) c(`as is` = rmse(v[!u], truth[!u]),
                                       clipped = rmse(pmin(v[!u], 0.5), truth[!u]))))
print(round(clipped, 5))

if (!is.null(pcrelate_npc)) {
  cat("\n=== sensitivity to the number of PCs (K=2, so 1 is correct) ===\n")
  sens <- sapply(1:4, function(d) c(
    `RMSE related` = rmse(pcrelate_npc[[d]][!u], truth[!u]),
    duplicates = mean(pcrelate_npc[[d]][cls == "duplicate"]),
    `cor with truth` = cor(pcrelate_npc[[d]], truth)))
  colnames(sens) <- paste(1:4, "PC"); cat("PC-Relate:\n"); print(round(sens, 5))
}

cat("\n=== implementation cross-checks (all should be 0) ===\n")
chk <- function(a, b) max(abs(a - b))
cat(sprintf("  --evaladmix-k 1 (4-col eigvecs) vs -k 1 : %.2e\n",
    chk(pc1, as.matrix(read.table(file.path(OUT, "pcaone_k4.kinship"), header = TRUE, check.names = FALSE)))))
cat(sprintf("  in-core vs out-of-core                  : %.2e\n",
    chk(as.matrix(read.table(file.path(OUT, "pcaone_ic.kinship"),  header = TRUE, check.names = FALSE)),
        as.matrix(read.table(file.path(OUT, "pcaone_ooc.kinship"), header = TRUE, check.names = FALSE)))))
if ("PCA + projection (R)" %in% names(est))
  cat(sprintf("  PCAone vs evalPCA() reference (r)       : %.6f\n",
      cor(est[["PCAone --evaladmix"]], est[["PCA + projection (R)"]])))

## ------------------------------------------------ detecting related pairs ---
# The tables above measure how close the estimates are. These ask the question
# a relatedness screen asks: is the pair related, and to what degree?
saveRDS(list(est = est, cls = cls, truth = truth), file.path(OUT, "estimates.rds"))

# KING degree bins (Manichaikul et al. 2010): degree d spans
# (2^-(d+1.5), 2^-(d+0.5)], extended here to 5th degree for second cousins
cut_deg <- 2^-(seq(0, 5) + 1.5)                      # 0.354 0.177 0.0884 0.0442 0.0221 0.0110
degree  <- function(v) findInterval(-v, -cut_deg)    # 0 = dup/MZ ... 5 = 5th, 6 = unrelated
true_deg <- c(unrelated = 6, duplicate = 0, "parent-offspring" = 1, "full sib" = 1,
              "half sib" = 2, "first cousin" = 3, "second cousin" = 5)[cls]

cat("\n=== degree called correctly (KING bins, 5th degree = second cousins) ===\n")
acc <- t(sapply(est, function(v) {
  dg <- degree(v); tapply(dg == true_deg, factor(cls, levels = c(CLS, "unrelated")), mean)
}))
print(round(acc, 2))

cat("\n=== unrelated pairs called related (false-positive rate) ===\n")
fp <- t(sapply(est, function(v) c(`>= 3rd degree` = mean(v[u] > cut_deg[4]),
                                  `>= 4th degree` = mean(v[u] > cut_deg[5]),
                                  `>= 5th degree` = mean(v[u] > cut_deg[6]),
                                  `max unrelated` = max(v[u]))))
print(signif(fp, 3))

cat("\n=== separation from the unrelated pairs ===\n")
# AUC: P(estimate of a related pair > estimate of an unrelated pair), per class.
# power: fraction of the class above the 99.9th percentile of the unrelated.
auc <- function(x, y) { r <- rank(c(x, y)); (sum(r[seq_along(x)]) - length(x) * (length(x) + 1) / 2) / (length(x) * length(y)) }
sep <- t(sapply(est, function(v) {
  t999 <- quantile(v[u], 0.999, names = FALSE)
  c(sapply(c("first cousin", "second cousin"), function(c) auc(v[cls == c], v[u])),
    sapply(c("first cousin", "second cousin"), function(c) mean(v[cls == c] > t999)))
}))
colnames(sep) <- c("AUC 1C", "AUC 2C", "power 1C @0.1%", "power 2C @0.1%")
print(round(sep, 3))
close <- sapply(est, function(v) min(sapply(CLS[1:5], function(c) auc(v[cls == c], v[u]))))
cat(sprintf("duplicates to half sibs: min AUC over the methods %.3f\n", min(close)))

## ------------------------------------------------------ missing genotypes ---
MISS <- file.path(OUT, "missing")
designs <- c(mcar0.05 = "5% at random", mcar0.1 = "10% at random", mcar0.2 = "20% at random",
             varying = "0-40% by sample", batch = "25% by batch")
designs <- designs[file.exists(file.path(MISS, paste0("pcaone_", names(designs), ".kinship")))]
if (length(designs)) {
  rd <- function(f, header) { M <- as.matrix(read.table(f, header = header, check.names = FALSE)); diag(M) <- 0; M }
  miss_est <- lapply(names(designs), function(d) list(
    `PCAone --evaladmix` = gp(rd(file.path(MISS, paste0("pcaone_", d, ".kinship")), TRUE)),
    `evalAdmix (EM)`     = gp(rd(file.path(MISS, paste0("evaladmix_em_", d, ".corres")), FALSE)) / 2))
  names(miss_est) <- names(designs)
  row <- function(v, keep = TRUE) c(`RMSE related` = rmse(v[!u & keep], truth[!u & keep]),
    `est/truth` = mean(v[!u & keep & cls != "second cousin"] / truth[!u & keep & cls != "second cousin"]),
    `degree right` = mean(degree(v[!u & keep]) == true_deg[!u & keep]),
    `RMSE unrelated` = rmse(v[u & keep], 0))
  cat("\n=== missing genotypes (complete data: see above) ===\n")
  tab <- do.call(rbind, lapply(names(designs), function(d) {
    t(sapply(miss_est[[d]], row))
  }))
  rownames(tab) <- paste(rep(designs, each = 2), rownames(tab), sep = " | ")
  print(round(tab, 5))
  if ("batch" %in% names(designs)) {
    b <- scan(file.path(MISS, "batch.batch"), quiet = TRUE)
    same <- b[pairs[, 1]] == b[pairs[, 2]]
    cat("\nbatch design, pairs in the same batch vs in different batches:\n")
    bt <- do.call(rbind, lapply(names(miss_est$batch), function(m) rbind(
      row(miss_est$batch[[m]], same), row(miss_est$batch[[m]], !same))))
    rownames(bt) <- paste(rep(names(miss_est$batch), each = 2), c("same batch", "different batches"), sep = " | ")
    print(round(bt, 5))
  }
}
