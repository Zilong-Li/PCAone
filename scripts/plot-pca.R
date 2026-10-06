#!/usr/bin/env Rscript
# Read only two PCs and aggregate every sample into a fixed-size raster.
# Source this file to use plot_pca() on the current graphics device.

plot_pca <- function(file, pcs = c(1L, 2L), format = c("auto", "ids", "matrix"),
                     bins = 700L, title = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE))
    stop('Install data.table first: install.packages("data.table")')
  format <- match.arg(format)
  if (length(pcs) != 2L || any(!is.finite(pcs)) || any(pcs < 1) ||
      any(pcs != as.integer(pcs)) || pcs[1] == pcs[2])
    stop("pcs must contain two distinct positive integers")
  if (length(bins) != 1L || !is.finite(bins) || bins < 2 || bins > 2000 ||
      bins != as.integer(bins)) stop("bins must be an integer from 2 to 2000")
  # PCAone .eigvecs2 has #FID IID PC1 ...; .eigvecs is a numeric matrix.
  probe <- data.table::fread(file, nrows = 0L, header = "auto")
  has_header <- all(paste0("PC", pcs) %in% names(probe))
  if (format == "auto") {
    format <- if (has_header || grepl("\\.(eigvecs2|eigenvec)(\\.gz)?$", file))
      "ids" else "matrix"
  }
  columns <- if (has_header) paste0("PC", pcs) else pcs + if (format == "ids") 2L else 0L
  d <- data.table::fread(file, select = columns, header = has_header,
                         showProgress = interactive())
  if (ncol(d) != 2L || !all(vapply(d, is.numeric, logical(1))))
    stop("Selected PC columns must be numeric; check --format and --pcs")
  keep <- is.finite(d[[1]]) & is.finite(d[[2]])
  excluded <- sum(!keep)
  if (excluded) warning(excluded, " samples with non-finite PCs excluded")
  x <- d[[1]][keep]; y <- d[[2]][keep]
  if (!length(x)) stop("No finite sample coordinates")
  bounds <- function(z) {
    r <- range(z)
    pad <- if (r[1] == r[2]) max(abs(r[1]), 1) * 0.01 else diff(r) * 0.02
    r + c(-pad, pad)
  }
  xr <- bounds(x); yr <- bounds(y)
  ix <- pmin(bins, pmax(1L, floor((x - xr[1]) / diff(xr) * bins) + 1L))
  iy <- pmin(bins, pmax(1L, floor((y - yr[1]) / diff(yr) * bins) + 1L))
  counts <- matrix(tabulate(ix + (iy - 1L) * bins, nbins = bins * bins), bins, bins)
  palette <- grDevices::hcl.colors(256L, "Inferno")
  intensity <- floor(log1p(counts) / log1p(max(counts)) * 255L) + 1L
  colors <- matrix(palette[intensity], bins, bins)
  colors[counts == 0L] <- "white"
  # Raster rows run top to bottom, whereas count matrix rows represent x.
  raster <- as.raster(t(colors[, bins:1L, drop = FALSE]))
  graphics::plot(NA_real_, NA_real_, xlim = xr, ylim = yr, xaxs = "i", yaxs = "i",
                 xlab = paste0("PC", pcs[1]), ylab = paste0("PC", pcs[2]),
                 main = if (is.null(title)) "PCA sample density" else title)
  graphics::rasterImage(raster, xr[1], yr[1], xr[2], yr[2], interpolate = FALSE)
  graphics::box()
  graphics::mtext(sprintf("%s samples | log-scaled count per bin | peak: %s",
                          format(length(x), big.mark = ","),
                          format(max(counts), big.mark = ",")), side = 3, cex = 0.7)
  # Compact colour key with counts, rather than transformed values.
  key <- unique(round(expm1(seq(0, log1p(max(counts)), length.out = 5))))
  graphics::legend("topright", legend = key,
                   fill = palette[floor(log1p(key) / log1p(max(counts)) * 255) + 1L],
                   title = "Samples/bin", bg = "white", cex = 0.75)
  invisible(list(samples = length(x), excluded = excluded, counts = counts,
                 xlim = xr, ylim = yr))
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste(
    "Usage: Rscript scripts/plot-pca.R INPUT [options]",
    "  INPUT                 PCAone .eigvecs2 or numeric .eigvecs file",
    "  -o, --out FILE        PNG or PDF (default: pca.png)",
    "  --pcs A,B             Two PCs (default: 1,2)",
    "  --format auto|ids|matrix  ids: first two columns are FID/IID",
    "  --bins N              Raster bins per axis, 2..2000 (default: 700)",
    "  --title TEXT          Plot title",
    "  -h, --help            Show this help", sep = "\n")
  if (!length(args) || any(args %in% c("-h", "--help"))) {
    cat(usage, "\n"); return(invisible(NULL))
  }
  file <- args[1]; out <- "pca.png"; pcs <- c(1, 2); fmt <- "auto"
  bins <- 700; title <- NULL
  i <- 2L
  while (i <= length(args)) {
    opt <- args[i]
    if (!opt %in% c("-o", "--out", "--pcs", "--format", "--bins", "--title"))
      stop("Unknown option: ", opt)
    if (i == length(args)) stop("Missing value for ", opt)
    value <- args[i + 1L]
    switch(opt, "-o" = , "--out" = { out <- value },
           "--pcs" = { pcs <- as.numeric(strsplit(value, ",", fixed = TRUE)[[1]]) },
           "--format" = { fmt <- value }, "--bins" = { bins <- as.numeric(value) },
           "--title" = { title <- value })
    i <- i + 2L
  }
  ext <- tolower(tools::file_ext(out))
  if (!ext %in% c("png", "pdf")) stop("Output must end in .png or .pdf")
  if (ext == "png") grDevices::png(out, width = 1800, height = 1500, res = 200)
  else grDevices::pdf(out, width = 9, height = 7.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  result <- plot_pca(file, pcs, fmt, bins, title)
  message("Plotted ", result$samples, " samples to ", out)
}

if (sys.nframe() == 0L) main()
