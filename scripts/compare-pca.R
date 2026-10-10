#!/usr/bin/env Rscript
# Compare sample eigenvectors, with O(N K^2) work and O(N K) memory.

read_pca_run <- function(file, k, format = "auto") {
  probe <- data.table::fread(file, nrows = 0, header = "auto")
  header <- "PC1" %in% names(probe)
  if (format == "auto") format <- if (header || grepl("\\.(eigvecs2|eigenvec)(\\.gz)?$", file)) "ids" else "matrix"
  if (!format %in% c("ids", "matrix")) stop("Format must be auto, ids or matrix")
  offset <- if (format == "ids") 2L else 0L
  if (ncol(probe) < k + offset) stop(file, ": fewer than ", k, " PCs")
  selected <- if (header) c(if (offset) names(probe)[1:2], paste0("PC", seq_len(k))) else seq_len(k + offset)
  d <- data.table::fread(file, select = selected, header = header,
                         colClasses = if (offset) list(character = selected[1:2]) else NULL)
  ids <- if (offset) d[, 1:2, with = FALSE] else NULL
  if (!is.null(ids)) {
    data.table::setnames(ids, c("FID", "IID"))
    if (anyNA(ids) || anyDuplicated(ids)) stop(file, ": missing or duplicate FID/IID pairs")
  }
  pc <- d[, seq_len(k) + offset, with = FALSE]
  if (!all(vapply(pc, is.numeric, logical(1)))) stop(file, ": PCs must be numeric")
  x <- as.matrix(pc)
  if (any(!is.finite(x))) stop(file, ": non-finite PC coordinates")
  list(ids = ids, x = x)
}

compare_pca <- function(file_a, file_b, k = 10L, format_a = "auto", format_b = "auto",
                        row_order = FALSE) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop('Install data.table: install.packages("data.table")')
  if (length(k) != 1 || !is.finite(k) || k < 1 || k != as.integer(k)) stop("k must be a positive integer")
  a <- read_pca_run(file_a, k, format_a); b <- read_pca_run(file_b, k, format_b)
  na <- nrow(a$x); nb <- nrow(b$x)
  if (!is.null(a$ids) && !is.null(b$ids)) {
    b$ids[, row_b := seq_len(.N)]
    index <- b$ids[a$ids, on = .(FID, IID), row_b]
    keep <- !is.na(index)
    a$x <- a$x[keep, , drop = FALSE]; b$x <- b$x[index[keep], , drop = FALSE]
  } else {
    if (!row_order) stop("Without IDs in both inputs, supply --row-order to confirm identical sample order")
    if (na != nb) stop("Row-order comparison requires equal sample counts")
  }
  n <- nrow(a$x)
  if (n <= k) stop("Need more shared samples than PCs")
  message("Shared samples: ", n, "; excluded from A: ", na - n, "; from B: ", nb - n)
  # Centre and unit-normalise each column: comparisons are scale independent.
  normalise <- function(x) {
    x <- sweep(x, 2, colMeans(x))
    norms <- sqrt(colSums(x * x))
    if (any(norms == 0)) stop("Constant PC column")
    sweep(x, 2, norms, "/")
  }
  x <- normalise(a$x); y <- normalise(b$x)
  cross <- crossprod(x, y)
  signs <- ifelse(diag(cross) < 0, -1, 1)
  # Solve min ||X - Y R||_F over orthogonal R (including reflections).
  fit <- svd(crossprod(y, x)); rotation <- fit$u %*% t(fit$v)
  aligned <- y %*% rotation
  metrics <- data.table::data.table(PC = seq_len(k), correlation = diag(cross),
    sign = signs, abs_correlation = abs(diag(cross)),
    sign_aligned_relative_error = sqrt(pmax(0, 2 - 2 * abs(diag(cross)))),
    rotation_aligned_relative_error = sqrt(colSums((x - aligned)^2)))
  # Whiten using small K x K Gram matrices, avoiding an N x N matrix.
  whitener <- function(z) {
    e <- eigen(crossprod(z), symmetric = TRUE)
    if (min(e$values) <= max(e$values) * 1e-10) stop("Selected PCs are rank deficient")
    sweep(e$vectors, 2, sqrt(e$values), "/")
  }
  cosines <- pmin(1, pmax(0, svd(t(whitener(x)) %*% cross %*% whitener(y), nu = 0, nv = 0)$d))
  angles <- acos(cosines) * 180 / pi
  summary <- data.table::data.table(shared_samples = n, only_a = na - n, only_b = nb - n,
    pcs = k, subspace_overlap = mean(cosines^2), max_principal_angle_deg = max(angles),
    procrustes_relative_error = sqrt(sum((x - aligned)^2) / k))
  list(metrics = metrics, summary = summary, correlation = cross,
       angles = data.table::data.table(direction = seq_len(k), angle_deg = angles))
}

plot_comparison <- function(result, labels = c("Run A", "Run B")) {
  old <- par(mfrow = c(1, 3), mar = c(5, 5, 3, 1)); on.exit(par(old))
  k <- nrow(result$metrics)
  image(seq_len(k), seq_len(k), abs(result$correlation), zlim = c(0, 1),
        col = grDevices::hcl.colors(100, "Blues 3"),
        xlab = paste(labels[1], "PC"), ylab = paste(labels[2], "PC"),
        main = "Absolute PC correlation (0-1)", axes = FALSE)
  axis(1, at = seq_len(k)); axis(2, at = seq_len(k)); box()
  matplot(result$metrics$PC, as.matrix(result$metrics[, c("sign_aligned_relative_error", "rotation_aligned_relative_error")]),
          type = "b", pch = c(16, 17), lty = 1, col = c("#0072B2", "#D55E00"),
          xlab = "PC", ylab = "Relative error (unit-norm PCs)", main = "Aligned coordinate error", ylim = c(0, 2))
  legend("topright", c("Sign only", "Orthogonal rotation"), col = c("#0072B2", "#D55E00"), pch = c(16, 17), bty = "n")
  plot(result$angles$direction, result$angles$angle_deg, type = "b", pch = 16,
       xlab = "Subspace direction (best to worst)", ylab = "Principal angle (degrees)",
       main = "Subspace agreement", ylim = c(0, 90))
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  help <- paste("Usage: Rscript scripts/compare-pca.R A.eigvecs2 B.eigvecs2 [options]",
    "  --k N              Compare first N PCs (default 10)",
    "  -o, --out PREFIX   Output prefix (default pca-comparison)",
    "  --labels A,B       Run names (default Run A,Run B)",
    "  --format-a FORMAT  auto, ids or matrix",
    "  --format-b FORMAT  auto, ids or matrix",
    "  --row-order        Confirm identical row order for inputs without IDs",
    "Writes PREFIX.pdf, .pcs.tsv, .angles.tsv and .summary.tsv", sep = "\n")
  if (!length(args) || any(args %in% c("-h", "--help"))) { cat(help, "\n"); return(invisible(NULL)) }
  if (length(args) < 2) stop(help)
  k <- 10; out <- "pca-comparison"; labels <- c("Run A", "Run B")
  fa <- fb <- "auto"; row_order <- FALSE; i <- 3L
  while (i <= length(args)) {
    opt <- args[i]
    if (opt == "--row-order") { row_order <- TRUE; i <- i + 1L; next }
    if (!opt %in% c("--k", "-o", "--out", "--labels", "--format-a", "--format-b")) stop("Unknown option: ", opt)
    if (i == length(args)) stop("Missing value for ", opt)
    v <- args[i + 1L]
    switch(opt, "--k" = { k <- as.numeric(v) }, "-o" = , "--out" = { out <- v },
           "--labels" = { labels <- strsplit(v, ",", fixed = TRUE)[[1]] },
           "--format-a" = { fa <- v }, "--format-b" = { fb <- v })
    i <- i + 2L
  }
  if (length(labels) != 2) stop("--labels needs two comma-separated names")
  r <- compare_pca(args[1], args[2], k, fa, fb, row_order)
  data.table::fwrite(r$metrics, paste0(out, ".pcs.tsv"), sep = "\t")
  data.table::fwrite(r$angles, paste0(out, ".angles.tsv"), sep = "\t")
  data.table::fwrite(r$summary, paste0(out, ".summary.tsv"), sep = "\t")
  grDevices::pdf(paste0(out, ".pdf"), width = 15, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_comparison(r, labels)
  print(r$summary)
}
if (sys.nframe() == 0L) main()
