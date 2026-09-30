#!/usr/bin/env Rscript
##
## plot-ld-decay.R: LD decay curves from PCAone -R/--print-r2 output.
##
## Reads one or more .ld.gz files (PCAone -R, or plink --ld with
## --ld-window-r2 0), bins the pairs by distance, and plots the mean r2 per bin,
## pooled over chromosomes. Compare standard and ancestry-adjusted LD:
##
##   PCAone -b plink -k 3 -o pcs                                   # the PCs
##   PCAone -b plink -P pcs -R --ld-bp 1000000 -o adj               # adjusted LD
##   PCAone -b plink --ld-stats 1 -R --ld-bp 1000000 -o std         # standard LD
##   Rscript plot-ld-decay.R adj.ld.gz std.ld.gz --labels Adjusted,Standard \
##       -n plink.fam -o ld-decay.png
##
## An input can also be a bin table, written by --save-bins or by the faster
## summarise_ld_r2bin tool (same bins, same format), so a large .ld.gz is
## binned once and re-plotted quickly.
##
## From R:
##   source("plot-ld-decay.R")
##   a <- ld_bins("adj.ld.gz"); s <- ld_bins("std.ld.gz")
##   plot_ld_decay(list(Adjusted = a, Standard = s), nsamples = 400)
##
## Needs R >= 4.0 and the data.table package.

LD_COLS <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
             "#e87ba4", "#008300", "#4a3aa7", "#e34948")
LD_LTY <- c(1, 2, 4, 5, 6, 3, 1, 2)
LD_INK <- "grey20"

## ---- binning ---------------------------------------------------------------

#' Bin edges: [0, min) and then `bins` log-spaced (or, with linear = TRUE,
#' `bins` equal) bins up to max. summarise_ld_r2bin.cpp uses the same edges.
ld_edges <- function(bins = 40L, min = 1000, max = 1e6, linear = FALSE) {
  if (linear) return(seq(0, max, length.out = bins + 1L))
  if (!(min > 0 && max > min)) stop("need 0 < --min < --max")
  c(0, 10^seq(log10(min), log10(max), length.out = bins + 1L))
}

## round a distance up to one significant digit: 999,998 -> 1e6, 2.9e6 -> 3e6
.round_up <- function(x) {
  p <- 10^floor(log10(x))
  ceiling(x / p - 1e-9) * p
}

.fread <- function(file, ...) {
  if (grepl("\\.gz$", file) && nzchar(Sys.which("gzip")))
    data.table::fread(cmd = paste("gzip -dc", shQuote(file)), showProgress = FALSE, ...)
  else data.table::fread(file, showProgress = FALSE, ...)
}

## column names from the first line (file() reads .gz transparently)
.header <- function(file) {
  con <- file(file, "r")
  on.exit(close(con))
  strsplit(trimws(readLines(con, n = 1L, warn = FALSE)), "[ \t]+")[[1]]
}

## TRUE if the file is a bin table rather than pairs
.is_bin_table <- function(file) all(c("chr", "lo", "hi", "n", "mean_r2") %in% .header(file))

#' Bin the pairs of one .ld/.ld.gz file, or read a bin table.
#'
#' Returns one row per chromosome and bin: chr, lo, hi, n, mean_dist, mean_r2,
#' sd_r2. Pairs on different chromosomes (plink --inter-chr) form one row with
#' chr "inter" and lo = hi = NA. Pairs beyond `max` are dropped with a message.
ld_bins <- function(file, bins = 40L, min = 1000, max = NULL, linear = FALSE) {
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("plot-ld-decay.R needs the data.table package: install.packages('data.table')")
  if (!file.exists(file)) stop("cannot find ", file)
  if (.is_bin_table(file)) {
    tab <- as.data.frame(.fread(file, colClasses = list(character = "chr")))
    attr(tab, "edges") <- sort(unique(c(tab$lo, tab$hi)))
    return(tab)
  }
  h <- .header(file)
  need <- c("CHR_A", "BP_A", "BP_B", "R2")
  if (!all(need %in% h))
    stop(file, " is neither an LD pair file (columns ", paste(need, collapse = ", "),
         ") nor a bin table")
  sel <- intersect(c("CHR_A", "BP_A", "CHR_B", "BP_B", "R2"), h)
  d <- .fread(file, select = sel, colClasses = list(character = intersect(c("CHR_A", "CHR_B"), sel)))
  if (!nrow(d)) stop(file, " has no pairs")
  data.table::set(d, j = "dist", value = abs(as.numeric(d$BP_B) - d$BP_A))
  inter <- if ("CHR_B" %in% sel) d$CHR_A != d$CHR_B else FALSE
  if (is.null(max)) max <- .round_up(base::max(if (any(inter)) d$dist[!inter] else d$dist, 1))
  e <- ld_edges(bins, min, max, linear)
  b <- findInterval(d$dist, e, rightmost.closed = TRUE)
  b[b == length(e)] <- -1L   # beyond max
  if (any(inter)) b[inter] <- 0L
  data.table::set(d, j = "bin", value = b)
  rm(b, inter)
  ## aggregate first, then tidy the (few) groups
  g <- d[, list(n = .N, s = sum(dist), r = sum(R2), rr = sum(R2 * R2)), by = list(chr = CHR_A, bin)]
  far <- g$bin == -1L
  if (any(far)) message(sprintf("%s: %s pairs beyond %s bp ignored; raise --max",
                                basename(file), format(sum(g$n[far]), big.mark = ","), format(max, scientific = FALSE)))
  g <- g[!far]
  g[, chr := ifelse(bin == 0L, "inter", sub("^chr", "", chr, ignore.case = TRUE))]
  g <- g[, list(n = sum(n), s = sum(s), r = sum(r), rr = sum(rr)), by = list(chr, bin)]
  tab <- data.frame(chr = g$chr,
                    lo = ifelse(g$bin > 0L, e[pmax(g$bin, 1L)], NA),
                    hi = ifelse(g$bin > 0L, e[pmax(g$bin, 1L) + 1L], NA),
                    n = g$n,
                    mean_dist = g$s / g$n,
                    mean_r2 = g$r / g$n,
                    sd_r2 = sqrt(pmax(0, (g$rr - g$r^2 / g$n) / pmax(g$n - 1, 1))),
                    stringsAsFactors = FALSE)
  tab <- tab[order(tab$chr == "inter", suppressWarnings(as.numeric(tab$chr)), tab$chr, tab$lo), ]
  rownames(tab) <- NULL
  attr(tab, "edges") <- e
  tab
}

## the same text as summarise_ld_r2bin writes
write_bins <- function(tab, file) {
  f <- function(x, d) ifelse(is.na(x), "NA", sprintf(paste0("%.", d, "g"), x))
  out <- data.frame(chr = tab$chr, lo = f(tab$lo, 15), hi = f(tab$hi, 15), n = sprintf("%.0f", tab$n),
                    mean_dist = f(tab$mean_dist, 7), mean_r2 = f(tab$mean_r2, 7), sd_r2 = f(tab$sd_r2, 7))
  utils::write.table(out, file, sep = "\t", quote = FALSE, row.names = FALSE)
}

#' Pool the chromosomes of a bin table into one curve (weighted by pairs).
#' Returns lo, hi, n, mean_dist, mean_r2, sd_r2, with attribute "inter", the
#' mean r2 of pairs on different chromosomes (NA if none).
pool_bins <- function(tab, chr = NULL) {
  inter <- tab[tab$chr == "inter", ]
  tab <- tab[tab$chr != "inter", ]
  if (length(chr)) {
    tab <- tab[tab$chr %in% sub("^chr", "", chr, ignore.case = TRUE), ]
    if (!nrow(tab)) stop("no bins left after --chr")
  }
  g <- data.table::as.data.table(tab)[, list(
    N = sum(n), sd = sum(n * mean_dist), sr = sum(n * mean_r2),
    ss = sum((n - 1) * sd_r2^2 + n * mean_r2^2)), by = list(lo, hi)]
  g <- g[order(lo)]
  m <- g$sr / g$N
  out <- data.frame(lo = g$lo, hi = g$hi, n = g$N, mean_dist = g$sd / g$N, mean_r2 = m,
                    sd_r2 = sqrt(pmax(0, (g$ss - g$N * m^2) / pmax(g$N - 1, 1))))
  attr(out, "inter") <- if (nrow(inter)) sum(inter$n * inter$mean_r2) / sum(inter$n) else NA_real_
  out
}

## sample-size correction of r2 (as in LD score regression): E[r2_hat] is
## about r2 + (1 - r2) / (N - 2), so r2 = r2_hat - (1 - r2_hat) / (N - 2)
correct_r2 <- function(r2, n) r2 - (1 - r2) / (n - 2)

## N from a number or a .fam/.psam file
.nsamples <- function(x) {
  if (is.null(x)) return(NULL)
  vapply(x, function(v) {
    if (file.exists(v)) {
      l <- readLines(v, warn = FALSE)
      sum(nzchar(l) & !startsWith(l, "#"))
    } else as.numeric(v)
  }, 0)
}

#' Summary of one pooled curve: pairs, r2 in the first bin with pairs, the
#' plateau (mean r2 of the bins reaching into the last 20% of the distance
#' range), and the
#' half-decay distance, where r2 first falls halfway from the first bin to the
#' plateau.
decay_summary <- function(curves) {
  do.call(rbind, lapply(names(curves), function(k) {
    p <- curves[[k]]
    top <- max(p$hi)
    far <- p$hi > 0.8 * top   # bins reaching into the last 20% of distances
    plateau <- sum(p$n[far] * p$mean_r2[far]) / sum(p$n[far])
    r0 <- p$mean_r2[1]
    half <- plateau + (r0 - plateau) / 2
    i <- which(p$mean_r2 <= half)[1]
    hd <- NA_real_
    if (!is.na(i) && i > 1) {   # interpolate on log distance
      x <- log10(p$mean_dist[c(i - 1, i)])
      y <- p$mean_r2[c(i - 1, i)]
      hd <- 10^(x[1] + (half - y[1]) * diff(x) / diff(y))
    }
    data.frame(curve = k, pairs = sum(p$n), first_bin_r2 = signif(r0, 4),
               plateau_r2 = signif(plateau, 4), half_decay_bp = round(hd),
               unlinked_r2 = signif(attr(p, "inter"), 4), stringsAsFactors = FALSE)
  }))
}

## ---- plot ------------------------------------------------------------------

.fmt_bp <- function(x) {
  vapply(x, function(v) {
    if (v == 0) "0"
    else if (v >= 1e6) paste(signif(v / 1e6, 6), "Mb")
    else if (v >= 1e3) paste(signif(v / 1e3, 6), "kb")
    else paste(signif(v, 6), "bp")
  }, "")
}

#' Plot LD decay curves.
#'
#' @param tabs     named list of bin tables (from ld_bins); names label curves
#' @param chr      chromosomes to pool (default all)
#' @param nsamples sample size, one or one per curve: draws the r2 expected
#'                 without LD, 1/(N-1), or with correct = TRUE removes it
#' @param correct  apply correct_r2() (needs nsamples)
#' @param baseline r2 of unlinked pairs, one per curve, drawn right of the plot;
#'                 by default taken from inter-chromosomal pairs in the input
#' @param by_chr   also draw each chromosome's curve, thin
#' @param ylog     log y axis
#' @return the pooled curves (list of data frames), invisibly
plot_ld_decay <- function(tabs, chr = NULL, nsamples = NULL, correct = FALSE, baseline = NULL,
                          by_chr = FALSE, ylog = FALSE, ymax = NULL, title = NULL, cex = 1,
                          legend_pos = "topright") {
  if (is.data.frame(tabs)) tabs <- list(tabs)
  if (is.null(names(tabs))) names(tabs) <- paste("curve", seq_along(tabs))
  k <- length(tabs)
  if (k > length(LD_COLS)) stop("at most ", length(LD_COLS), " curves")
  nsamples <- .nsamples(nsamples)
  if (length(nsamples) == 1L) nsamples <- rep(nsamples, k)
  if (length(nsamples) && length(nsamples) != k) stop("give one sample size, or one per curve")
  if (correct && !length(nsamples)) stop("--correct needs the sample size (-n)")
  if (length(baseline) && length(baseline) != k) stop("give one --baseline per curve")

  if (length(chr)) chr <- sub("^chr", "", chr, ignore.case = TRUE)
  curves <- lapply(tabs, pool_bins, chr = chr)
  per_chr <- if (by_chr) lapply(tabs, function(t) {
    t <- t[t$chr != "inter" & (!length(chr) | t$chr %in% chr), ]
    split(t, factor(t$chr, levels = unique(t$chr)))
  })
  base <- vapply(seq_len(k), function(i) if (length(baseline)) as.numeric(baseline[i]) else attr(curves[[i]], "inter"), 0)
  for (i in seq_len(k)) {
    if (correct) {
      Ni <- nsamples[i]
      curves[[i]]$mean_r2 <- correct_r2(curves[[i]]$mean_r2, Ni)
      if (by_chr) per_chr[[i]] <- lapply(per_chr[[i]], function(t) {
        t$mean_r2 <- correct_r2(t$mean_r2, Ni)
        t
      })
      base[i] <- correct_r2(base[i], Ni)
    }
    attr(curves[[i]], "inter") <- base[i]
  }
  expect <- if (length(nsamples) && !correct) unique(1 / (nsamples - 1)) else NULL

  ## axes
  edges <- sort(unique(unlist(lapply(curves, function(p) c(p$lo, p$hi)))))
  w <- unlist(lapply(curves, function(p) p$hi - p$lo))
  linear <- length(w) > 2L && diff(range(w)) < 1e-6 * max(w)   # equal-width bins
  x_all <- unlist(lapply(curves, `[[`, "mean_dist"))
  xlim <- if (linear) c(0, max(edges)) else c(max(min(x_all), 1), max(edges))
  y_all <- c(unlist(lapply(curves, `[[`, "mean_r2")), if (by_chr) unlist(lapply(per_chr, function(l) lapply(l, `[[`, "mean_r2"))),
             base, expect)
  y_all <- y_all[is.finite(y_all)]
  if (ylog) {
    y_pos <- y_all[y_all > 0]
    ylim <- c(min(y_pos) / 1.2, if (is.null(ymax)) max(y_pos) * 1.2 else ymax)
  } else {
    ylim <- c(min(0, y_all), if (is.null(ymax)) max(y_all) * 1.05 else ymax)
  }
  has_base <- any(is.finite(base))
  op <- par(mar = c(3.4, 4.2, if (is.null(title)) 1 else 2.4, if (has_base) 5.2 else 1),
            mgp = c(2.3, 0.5, 0), tcl = -0.25, las = 1)
  on.exit(par(op))
  plot.new()
  plot.window(xlim, ylim, log = paste0(if (linear) "" else "x", if (ylog) "y" else ""))
  grid_x <- if (linear) pretty(xlim) else 10^(ceiling(log10(xlim[1])):floor(log10(xlim[2])))
  abline(v = grid_x, col = "grey92", lwd = 1)
  if (!is.null(expect)) {
    abline(h = expect, col = "grey55", lty = 3, lwd = 1.2)
    text(grconvertX(0.01, "npc", "user"), expect, "1/(N-1)", adj = c(0, -0.4), cex = 0.75 * cex, col = "grey40")
  }
  for (i in seq_len(k)) {
    if (by_chr) for (t in per_chr[[i]]) lines(t$mean_dist, t$mean_r2, col = adjustcolor(LD_COLS[i], 0.25), lwd = 0.8)
    p <- curves[[i]]
    lines(p$mean_dist, p$mean_r2, col = LD_COLS[i], lty = LD_LTY[i], lwd = 2)
    points(p$mean_dist, p$mean_r2, col = LD_COLS[i], pch = 16, cex = 0.45 * cex)
  }
  ## unlinked pairs: short segments right of the plot
  if (has_base) {
    x0 <- grconvertX(1.04, "npc", "user")
    x1 <- grconvertX(1.16, "npc", "user")
    for (i in which(is.finite(base))) segments(x0, base[i], x1, base[i], col = LD_COLS[i], lty = LD_LTY[i], lwd = 2, xpd = NA)
    text(grconvertX(1.10, "npc", "user"), grconvertY(0, "npc", "user"), "unlinked", adj = c(0.5, 1.6),
         cex = 0.75 * cex, col = LD_INK, xpd = NA)
  }
  if (linear) {
    at <- pretty(xlim)
    axis(1, at = at, labels = .fmt_bp(at), col = "grey45", col.axis = LD_INK, cex.axis = 0.85 * cex)
  } else {
    axis(1, at = grid_x, labels = .fmt_bp(grid_x), col = "grey45", col.axis = LD_INK, cex.axis = 0.85 * cex)
    minor <- as.vector(outer(2:9, 10^(floor(log10(xlim[1])):floor(log10(xlim[2])))))
    axis(1, at = minor[minor >= xlim[1] & minor <= xlim[2]], labels = FALSE, tcl = -0.15, col = "grey45")
  }
  axis(2, col = "grey45", col.axis = LD_INK, cex.axis = 0.85 * cex)
  box(bty = "l", col = "grey45")
  title(xlab = "Distance between SNPs", col.lab = LD_INK)
  title(ylab = if (correct) expression("mean " * r^2 * ", corrected for sample size") else expression("mean " * r^2),
        col.lab = LD_INK)
  if (k > 1L)
    legend(legend_pos, legend = names(tabs), col = LD_COLS[seq_len(k)], lty = LD_LTY[seq_len(k)], lwd = 2,
           bty = "n", cex = 0.9 * cex, text.col = LD_INK, seg.len = 2.5)
  if (!is.null(title)) mtext(title, side = 3, line = 0.8, adj = 0, font = 2, cex = 1.05 * cex, col = LD_INK)
  invisible(curves)
}

## ---- command line ----------------------------------------------------------

LD_USAGE <- "Plot LD decay curves from PCAone -R/--print-r2 output.

Usage: Rscript plot-ld-decay.R [options] INPUT [INPUT ...]

INPUT is an .ld.gz/.ld file (PCAone -R, or plink --ld --ld-window-r2 0), or a
bin table written by --save-bins or summarise_ld_r2bin. Each input is a curve.

Binning (ignored for bin tables, which keep their bins)
      --bins N        log-spaced bins between --min and --max [40]
      --min D         lower edge of the first log bin, bp; closer pairs form one bin [1000]
      --max D         upper edge of the last bin, bp; set it to --ld-bp
                      [largest distance, rounded up]
      --linear        equal-width bins from 0 to --max, and a linear x axis
      --save-bins     write each input's bins to <input>.decay.tsv (re-plot from it)

Plot
  -o, --out F         .png, .pdf, .svg, .jpg or .tiff [ld-decay.png]
      --labels L      comma-separated curve names [input file names]
      --chr S         chromosomes to pool, e.g. 1-22 [all]
  -n, --nsamples N    sample size, or a .fam/.psam; one, or one per input.
                      Draws the r2 expected without LD, 1/(N-1)
      --correct       subtract the sampling expectation: r2 - (1 - r2)/(N - 2)
      --baseline V    r2 of unlinked pairs, one per input, drawn right of the
                      plot [mean r2 of inter-chromosomal pairs in the input]
      --by-chr        also draw each chromosome, thin
      --ylog          log y axis
      --ymax Y        upper y limit
      --legend P      legend position: topright, right, bottomleft, ... [topright]
      --title T       plot title
      --width W       inches [6]
      --height H      inches [4.5]
      --res N         png/jpg/tiff resolution [150]
      --cex X         scale text [1]
      --table F       write the plotted (pooled) curves to F (tsv)
  -h, --help

The summary printed per curve: pairs, mean r2 of the first bin, the plateau
(mean r2 over the last 20% of distances), the half-decay distance (where r2
falls halfway from the first bin to the plateau), and the unlinked r2."

.parse_set <- function(s) {
  if (is.null(s)) return(NULL)
  tok <- strsplit(gsub("\\s+", "", s), ",", fixed = TRUE)[[1]]
  unique(unlist(lapply(tok[nzchar(tok)], function(x) {
    ab <- strsplit(x, "-", fixed = TRUE)[[1]]
    if (length(ab) == 2L && all(grepl("^[0-9]+$", ab))) as.character(seq(as.integer(ab[1]), as.integer(ab[2]))) else x
  })))
}

.split <- function(s) if (is.null(s)) NULL else trimws(strsplit(s, ",", fixed = TRUE)[[1]])

.bp <- function(x) {
  if (is.null(x)) return(NULL)
  x <- sub("B$", "", toupper(trimws(x)))
  mult <- c(K = 1e3, M = 1e6)
  u <- substring(x, nchar(x))
  v <- if (u %in% names(mult)) as.numeric(substr(x, 1L, nchar(x) - 1L)) * mult[[u]] else as.numeric(x)
  if (is.na(v)) stop("cannot parse distance '", x, "'", call. = FALSE)
  v
}

parse_ld_args <- function(args) {
  short <- c(o = "out", n = "nsamples", h = "help")
  flags <- c("help", "linear", "save-bins", "correct", "by-chr", "ylog")
  known <- c("out", "labels", "bins", "min", "max", "chr", "nsamples", "baseline", "ymax", "legend",
             "title", "width", "height", "res", "cex", "table", flags)
  opt <- list(inputs = character())
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    val <- NULL
    if (grepl("^--[^=]+=", a)) {
      key <- sub("^--([^=]+)=.*$", "\\1", a)
      val <- sub("^--[^=]+=", "", a)
    } else if (grepl("^--.", a)) {
      key <- substring(a, 3L)
    } else if (grepl("^-[A-Za-z]$", a)) {
      if (!substring(a, 2L) %in% names(short)) stop("unknown option ", a, " (see -h)", call. = FALSE)
      key <- short[[substring(a, 2L)]]
    } else {
      opt$inputs <- c(opt$inputs, a)
      i <- i + 1L
      next
    }
    if (!key %in% known) stop("unknown option --", key, " (see -h)", call. = FALSE)
    if (key %in% flags) {
      opt[[key]] <- TRUE
    } else {
      if (is.null(val)) {
        i <- i + 1L
        if (i > length(args)) stop("--", key, " needs a value", call. = FALSE)
        val <- args[i]
      }
      opt[[key]] <- val
    }
    i <- i + 1L
  }
  opt
}

.open_device <- function(file, width, height, res = 150) {
  ext <- tolower(sub(".*\\.", "", basename(file)))
  switch(ext,
    png = grDevices::png(file, width = width, height = height, units = "in", res = res),
    jpg = , jpeg = grDevices::jpeg(file, width = width, height = height, units = "in", res = res, quality = 95),
    tif = , tiff = grDevices::tiff(file, width = width, height = height, units = "in", res = res, compression = "lzw"),
    pdf = grDevices::pdf(file, width = width, height = height, useDingbats = FALSE),
    svg = grDevices::svg(file, width = width, height = height),
    stop("unknown output format '.", ext, "'; use .png, .pdf, .svg, .jpg or .tiff", call. = FALSE))
}

ld_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  o <- parse_ld_args(args)
  if (isTRUE(o$help) || !length(args)) {
    cat(LD_USAGE, "\n")
    return(invisible())
  }
  if (!length(o$inputs)) stop("give at least one .ld.gz or bin table (see -h)", call. = FALSE)
  labels <- .split(o$labels)
  if (is.null(labels)) labels <- sub("\\.(ld(\\.gz)?|decay\\.tsv|tsv|txt)$", "", basename(o$inputs))
  if (length(labels) != length(o$inputs)) stop("give one label per input", call. = FALSE)
  num <- function(x, def) if (is.null(x)) def else as.numeric(x)

  bins <- as.integer(num(o$bins, 40))
  mn <- num(.bp(o$min), 1000)
  mx <- .bp(o$max)
  linear <- isTRUE(o$linear)
  ## one set of edges for all inputs: without --max, the first .ld.gz sets it
  tabs <- list()
  for (i in seq_along(o$inputs)) {
    f <- o$inputs[i]
    t0 <- Sys.time()
    tabs[[labels[i]]] <- ld_bins(f, bins = bins, min = mn, max = mx, linear = linear)
    if (is.null(mx) && !.is_bin_table(f)) mx <- max(attr(tabs[[labels[i]]], "edges"))
    message(sprintf("%s: %s pairs in %d bins (%.1f s)", f, format(sum(tabs[[i]]$n), big.mark = ","),
                    length(unique(tabs[[i]]$lo[!is.na(tabs[[i]]$lo)])), as.numeric(Sys.time() - t0, units = "secs")))
    if (isTRUE(o$`save-bins`) && !.is_bin_table(f)) {
      fb <- paste0(sub("\\.ld(\\.gz)?$", "", f), ".decay.tsv")
      write_bins(tabs[[labels[i]]], fb)
      message("wrote ", fb)
    }
  }
  out <- if (is.null(o$out)) "ld-decay.png" else o$out
  .open_device(out, num(o$width, 6), num(o$height, 4.5), num(o$res, 150))
  curves <- tryCatch(
    plot_ld_decay(tabs, chr = .parse_set(o$chr), nsamples = .split(o$nsamples), correct = isTRUE(o$correct),
                  baseline = .split(o$baseline), by_chr = isTRUE(o$`by-chr`), ylog = isTRUE(o$ylog),
                  ymax = if (is.null(o$ymax)) NULL else as.numeric(o$ymax), title = o$title,
                  cex = num(o$cex, 1), legend_pos = if (is.null(o$legend)) "topright" else o$legend),
    finally = grDevices::dev.off())
  message("wrote ", out)
  print(decay_summary(curves), row.names = FALSE)
  if (!is.null(o$table)) {
    tt <- do.call(rbind, lapply(names(curves), function(k) cbind(curve = k, curves[[k]])))
    for (v in c("mean_dist", "mean_r2", "sd_r2")) tt[[v]] <- signif(tt[[v]], 7)
    utils::write.table(tt, o$table, sep = "\t", quote = FALSE, row.names = FALSE)
    message("wrote ", o$table)
  }
  invisible(curves)
}

if (sys.nframe() == 0L && !interactive()) {
  tryCatch(ld_main(), error = function(e) {
    message("error: ", conditionMessage(e))
    quit(status = 1)
  })
}
