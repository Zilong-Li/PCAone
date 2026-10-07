#!/usr/bin/env Rscript
# PCAone per-site HWE diagnostics, with bounded plotting work for large files.
# Public sourceable API: read_hwe(), plot_hwe(), save_hwe().

read_hwe <- function(file, bim = NULL, p_floor = 1e-300) {
  if (!requireNamespace("data.table", quietly = TRUE))
    stop('Install data.table: install.packages("data.table")')
  if (!is.finite(p_floor) || p_floor <= 0 || p_floor >= 1) stop("p_floor must be between 0 and 1")
  d <- data.table::fread(file, colClasses = list(character = 1L))
  data.table::setnames(d, sub("^#", "", names(d)))
  required <- c("ID", "HWE_P", "LRT", "Inbreeding_coefficient")
  if (!all(required %in% names(d))) stop("Expected PCAone .hwe columns: ", paste(required, collapse = ", "))
  if (!all(vapply(d[, required[-1], with = FALSE], is.numeric, logical(1))))
    stop("HWE_P, LRT and Inbreeding_coefficient must be numeric")
  if (!is.null(bim)) {
    b <- data.table::fread(bim, header = FALSE, select = c(1, 2, 4),
                           colClasses = list(character = c(1, 2)))
    if (nrow(b) != nrow(d) || !identical(b[[2]], d$ID))
      stop("BIM/MBIM must contain exactly the HWE variants in the same order")
    if (!is.numeric(b[[3]]) || any(!is.finite(b[[3]]) | b[[3]] < 0) || anyNA(b[[1]]))
      stop("Invalid chromosome or position in BIM/MBIM")
    d[, `:=`(chr = b[[1]], bp = b[[3]])]
  }
  valid <- is.finite(d$HWE_P) & d$HWE_P >= 0 & d$HWE_P <= 1
  if (any(!valid)) warning(sum(!valid), " variants with invalid HWE P excluded")
  d <- d[valid]
  if (!nrow(d)) stop("No valid HWE P values")
  clipped <- sum(d$HWE_P < p_floor)
  if (clipped) warning(clipped, " P values capped at ", p_floor, " for plotting (including zeros)")
  d[, logp := -log10(pmax(HWE_P, p_floor))]
  attr(d, "excluded") <- sum(!valid)
  attr(d, "clipped") <- clipped
  d
}

plot_hwe <- function(d, alpha = 0.05, title = "HWE diagnostics", max_points = 5000L) {
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  if (!is.finite(max_points) || max_points < 100 || max_points != as.integer(max_points))
    stop("max_points must be an integer >= 100")
  n <- nrow(d); cutoff <- alpha / n
  genomic <- all(c("chr", "bp") %in% names(d))
  old <- par(mfrow = if (genomic) c(2, 2) else c(1, 3), mar = c(4.5, 4.5, 3, 1), oma = c(0, 0, 2, 0))
  on.exit(par(old))
  # Exact sorted P values at logarithmically spaced ranks; retain extreme tail.
  observed <- sort(d$logp, decreasing = TRUE)
  rank <- unique(as.integer(round(exp(seq(0, log(n), length.out = min(n, max_points))))))
  expected <- -log10((rank - 0.5) / n)
  plot(expected, observed[rank], pch = 16, cex = 0.45, col = "#0072B2",
       xlim = c(0, max(expected)), ylim = c(0, max(expected, observed)),
       xlab = "Expected -log10(P)", ylab = "Observed -log10(P)", main = "HWE Q-Q")
  abline(0, 1, col = "grey40", lty = 2)
  hist(d$HWE_P, breaks = seq(0, 1, length.out = 51), col = "#0072B2", border = "white",
       main = "HWE P values", xlab = "P", ylab = "Variants")
  f <- d$Inbreeding_coefficient[is.finite(d$Inbreeding_coefficient)]
  if (length(f)) {
    hist(f, breaks = 80, col = "#D55E00", border = "white", main = "Per-site inbreeding",
         xlab = "Inbreeding coefficient F", ylab = "Variants")
    abline(v = 0, lty = 2)
  } else { plot.new(); title(main = "No finite inbreeding coefficients") }
  if (genomic) {
    # Chromosomes follow their order of first appearance in the input.
    chr <- unique(d$chr)
    lengths <- vapply(chr, function(ch) max(d$bp[d$chr == ch]), numeric(1))
    lengths <- pmax(lengths, 1)
    gap <- max(sum(lengths) * 0.005, 1)
    offsets <- c(0, head(cumsum(lengths + gap), -1))
    ci <- match(d$chr, chr)
    x <- d$bp + offsets[ci]
    # Keep the strongest signal in each position bin, separately by chromosome.
    bin <- floor(x / (sum(lengths) + gap * length(chr)) * max_points)
    b <- data.table::data.table(x = x, y = d$logp, chr = ci, bin = bin)
    peaks <- b[, .(x = x[which.max(y)], y = max(y)), by = .(chr, bin)]
    plot(peaks$x, peaks$y, pch = 16, cex = 0.5,
         col = c("#0072B2", "#666666")[(peaks$chr - 1) %% 2 + 1],
         xaxt = "n", xlab = "Chromosome", ylab = "-log10(HWE P)", main = "HWE across the genome",
         ylim = c(0, max(d$logp, -log10(cutoff))))
    axis(1, at = offsets + lengths / 2, labels = chr, cex.axis = 0.7)
    abline(h = -log10(cutoff), col = "#D55E00", lty = 2)
    mtext("Dashed: Bonferroni alpha / valid variants", side = 3, cex = 0.65)
  }
  mtext(sprintf("%s | %s valid variants | %s P < %.3g", title,
                format(n, big.mark = ","), format(sum(d$HWE_P < cutoff), big.mark = ","), cutoff),
        outer = TRUE, cex = 0.9)
  invisible(data.table::data.table(valid_variants = n, excluded_p = attr(d, "excluded"),
    capped_p = attr(d, "clipped"), bonferroni_cutoff = cutoff,
    below_cutoff = sum(d$HWE_P < cutoff), finite_f = length(f),
    median_f = if (length(f)) median(f) else NA_real_))
}

# Save diagnostics to a new device, close only that device, and return summary.
# d is the table returned by read_hwe(); summary_out = NULL disables the TSV.
save_hwe <- function(d, out = "hwe.png", alpha = 0.05,
                     title = "HWE diagnostics", max_points = 5000L,
                     width = NULL, height = NULL, res = 180,
                     summary_out = paste0(out, ".summary.tsv")) {
  ext <- tolower(tools::file_ext(out))
  if (!ext %in% c("png", "pdf")) stop("Output must be .png or .pdf")
  genomic <- all(c("chr", "bp") %in% names(d))
  if (is.null(width)) width <- if (genomic) 11 else 15
  if (is.null(height)) height <- if (genomic) 8 else 5
  if (ext == "png") grDevices::png(out, width = width, height = height, units = "in", res = res)
  else grDevices::pdf(out, width = width, height = height)
  device <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(device), add = TRUE)
  summary <- plot_hwe(d, alpha, title, max_points)
  if (!is.null(summary_out)) data.table::fwrite(summary, summary_out, sep = "\t")
  invisible(summary)
}

.hwe_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste("Usage: Rscript scripts/plot-hwe.R INPUT.hwe [options]",
    "  -o, --out FILE    PNG or PDF (default INPUT.hwe.png)",
    "  --bim FILE        Matching .bim/.mbim; adds a genome plot",
    "  --alpha A         Bonferroni family-wise threshold (default 0.05)",
    "  --p-floor P       Plotting floor for tiny/zero P (default 1e-300)",
    "  --max-points N    Q-Q ranks / genome bins (default 5000)",
    "  --title TEXT      Plot title", sep = "\n")
  if (!length(args) || any(args %in% c("-h", "--help"))) { cat(usage, "\n"); return(invisible(NULL)) }
  file <- args[1]; out <- paste0(file, ".png"); bim <- NULL
  alpha <- 0.05; floor <- 1e-300; points <- 5000; title <- "HWE diagnostics"; i <- 2L
  while (i <= length(args)) {
    opt <- args[i]
    if (!opt %in% c("-o", "--out", "--bim", "--alpha", "--p-floor", "--max-points", "--title")) stop("Unknown option: ", opt)
    if (i == length(args)) stop("Missing value for ", opt)
    v <- args[i + 1L]
    switch(opt, "-o" = , "--out" = { out <- v }, "--bim" = { bim <- v },
      "--alpha" = { alpha <- as.numeric(v) }, "--p-floor" = { floor <- as.numeric(v) },
      "--max-points" = { points <- as.numeric(v) }, "--title" = { title <- v })
    i <- i + 2L
  }
  d <- read_hwe(file, bim, floor)
  summary <- save_hwe(d, out, alpha, title, points)
  print(summary)
  message("Saved ", out)
}
if (sys.nframe() == 0L) .hwe_main()
