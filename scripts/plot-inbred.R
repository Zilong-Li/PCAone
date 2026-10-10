#!/usr/bin/env Rscript
# Public API: read_inbred(), plot_inbred(), save_inbred().
# Per-sample diagnostics for PCAone --inbreed 2.

read_inbred <- function(file) {
  if (!requireNamespace("data.table", quietly = TRUE))
    stop('Install data.table: install.packages("data.table")')
  d <- data.table::fread(file, colClasses = list(character = 1:2))
  data.table::setnames(d, sub("^#", "", names(d)))
  required <- c("FID", "IID", "N_SITES", "O_HET", "E_HET", "F")
  if (!all(required %in% names(d))) stop("Expected .inbred columns: ", paste(required, collapse = ", "))
  if (!nrow(d)) stop("No samples in input")
  # fread infers an all-NA F column as logical when every sample has no data.
  if (is.logical(d$F) && all(is.na(d$F))) d[, F := as.numeric(F)]
  if (!all(vapply(d[, required[3:6], with = FALSE], is.numeric, logical(1))))
    stop("N_SITES, O_HET, E_HET and F must be numeric")
  if (anyNA(d[, .(FID, IID)]) || anyDuplicated(d[, .(FID, IID)]))
    stop("Missing or duplicate FID/IID pairs")
  if (any(!is.finite(d$N_SITES) | d$N_SITES < 0 | d$N_SITES != floor(d$N_SITES)))
    stop("N_SITES must be finite nonnegative integers")
  if (any(is.finite(d$F) & abs(d$F) > 1)) warning("Some F values are outside [-1, 1]")
  if (any(!is.finite(d$F))) warning(sum(!is.finite(d$F)), " samples with non-finite F omitted from F plots")
  d
}

# A fixed-resolution density raster avoids drawing hundreds of thousands of points.
.inbred_density <- function(x, y, xlab, ylab, main) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  if (!length(x)) { plot.new(); title(main = paste(main, "(no valid samples)")); return(invisible(NULL)) }
  bounds <- function(z) {
    r <- range(z); p <- if (diff(r) == 0) max(abs(r), 1) * 0.01 else diff(r) * 0.02
    r + c(-p, p)
  }
  xr <- bounds(x); yr <- bounds(y); bins <- 300L
  ix <- pmin(bins, pmax(1L, floor((x - xr[1]) / diff(xr) * bins) + 1L))
  iy <- pmin(bins, pmax(1L, floor((y - yr[1]) / diff(yr) * bins) + 1L))
  count <- matrix(tabulate(ix + (iy - 1L) * bins, nbins = bins * bins), bins, bins)
  pal <- grDevices::hcl.colors(256, "Inferno")
  color <- matrix(pal[1L + floor(log1p(count) / log1p(max(count)) * 255)], bins, bins)
  color[count == 0] <- "white"
  plot(NA_real_, NA_real_, xlim = xr, ylim = yr, xaxs = "i", yaxs = "i",
       xlab = xlab, ylab = ylab, main = main)
  rasterImage(as.raster(t(color[, bins:1, drop = FALSE])), xr[1], yr[1], xr[2], yr[2], interpolate = FALSE)
  box()
  mtext(sprintf("Log count/bin; peak %s; n = %s", format(max(count), big.mark = ","),
                format(length(x), big.mark = ",")), side = 3, cex = 0.65)
}

plot_inbred <- function(d, title = "Per-sample inbreeding", breaks = 100L) {
  if (length(breaks) != 1 || !is.finite(breaks) || breaks < 2 || breaks != as.integer(breaks))
    stop("breaks must be an integer >= 2")
  old <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.5, 1), oma = c(0, 0, 2, 0))
  on.exit(par(old))
  valid_f <- is.finite(d$F) & d$N_SITES > 0 & is.finite(d$E_HET) & d$E_HET > 0
  f <- d$F[valid_f]
  if (length(f)) {
    hist(f, breaks = breaks, col = "#0072B2", border = "white",
         xlab = "Sample inbreeding coefficient F", main = "F distribution")
    abline(v = 0, lty = 2, col = "grey40")
  } else { plot.new(); title(main = "No valid F estimates") }
  hist(d$N_SITES, breaks = breaks, col = "#666666", border = "white",
       xlab = "Number of informative sites", main = "Site coverage")
  het <- d$N_SITES > 0 & is.finite(d$O_HET) & d$O_HET >= 0 &
    is.finite(d$E_HET) & d$E_HET > 0
  .inbred_density(d$E_HET[het] / d$N_SITES[het], d$O_HET[het] / d$N_SITES[het],
                   "Expected heterozygosity (E_HET / N_SITES)",
                   "Observed heterozygosity (O_HET / N_SITES)", "Observed vs expected")
  if (any(het)) abline(0, 1, lty = 2, col = "#009E73")
  .inbred_density(d$N_SITES[valid_f], f, "Number of informative sites", "F", "F vs site coverage")
  if (length(f)) abline(h = 0, lty = 2, col = "#009E73")
  mtext(sprintf("%s | %s samples | %s valid F estimates", title,
                format(nrow(d), big.mark = ","), format(length(f), big.mark = ",")), outer = TRUE)
  invisible(data.table::data.table(samples = nrow(d), valid_f = length(f),
    excluded_f = nrow(d) - length(f), zero_sites = sum(d$N_SITES == 0),
    median_sites = median(d$N_SITES),
    median_f = if (length(f)) median(f) else NA_real_,
    f_q01 = if (length(f)) unname(quantile(f, 0.01)) else NA_real_,
    f_q99 = if (length(f)) unname(quantile(f, 0.99)) else NA_real_))
}

save_inbred <- function(d, out = "inbred.png", title = "Per-sample inbreeding",
                        breaks = 100L, width = 11, height = 8, res = 180,
                        summary_out = paste0(out, ".summary.tsv")) {
  ext <- tolower(tools::file_ext(out))
  if (!ext %in% c("png", "pdf")) stop("Output must be .png or .pdf")
  if (ext == "png") grDevices::png(out, width = width, height = height, units = "in", res = res)
  else grDevices::pdf(out, width = width, height = height)
  device <- grDevices::dev.cur(); on.exit(grDevices::dev.off(device), add = TRUE)
  summary <- plot_inbred(d, title, breaks)
  if (!is.null(summary_out)) data.table::fwrite(summary, summary_out, sep = "\t")
  invisible(summary)
}

.inbred_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  help <- paste("Usage: Rscript scripts/plot-inbred.R INPUT.inbred [options]",
    "  -o, --out FILE    PNG or PDF (default INPUT.inbred.png)",
    "  --title TEXT      Plot title",
    "  --breaks N        Histogram target bin count (default 100)", sep = "\n")
  if (!length(args) || any(args %in% c("-h", "--help"))) { cat(help, "\n"); return(invisible(NULL)) }
  file <- args[1]; out <- paste0(file, ".png"); title <- "Per-sample inbreeding"; breaks <- 100L; i <- 2L
  while (i <= length(args)) {
    opt <- args[i]
    if (!opt %in% c("-o", "--out", "--title", "--breaks")) stop("Unknown option: ", opt)
    if (i == length(args)) stop("Missing value for ", opt)
    v <- args[i + 1L]
    switch(opt, "-o" = , "--out" = { out <- v }, "--title" = { title <- v },
           "--breaks" = { breaks <- as.numeric(v) })
    i <- i + 2L
  }
  summary <- save_inbred(read_inbred(file), out, title, breaks)
  print(summary); message("Saved ", out)
}
if (sys.nframe() == 0L) .inbred_main()
