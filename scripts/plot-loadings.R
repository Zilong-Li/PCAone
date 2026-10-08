#!/usr/bin/env Rscript
##
## plot-loadings.R: plot PCAone SNP loadings along the genome.
##
## A PCAone run writes <prefix>.loadings (one row per variant, one column
## per PC) and <prefix>.mbim (chr, id, cM, bp, A1, A2, freq). For each PC the
## script weights each PC by its singular value from <prefix>.sigvals and
## plots the max |weighted loading| in consecutive bins of variants. Bins never
## span two chromosomes, so millions of variants plot quickly and no peak is
## lost. Each PC's top variant is marked and printed as a table.
##
## Command line (see -h for all options):
##   Rscript plot-loadings.R -p pcaone
##   Rscript plot-loadings.R -p pcaone --mode panel --pcs 1-8
##   Rscript plot-loadings.R -p pcaone --highlight 1,2 \
##       --groups "Population structure:1,2,4;HLA:6,12"
##   Rscript plot-loadings.R -p pcaone --region 6:25M-35M --mode panel --pcs 1-8
##
## The categories of the original 40-PC figure:
##   --groups "Population structure:1,2,4;Centromere:3,5,8-11,13,14,19,20,22,29,32,34,35,37;
##             HLA:6,12,15,16,27,30,38;Other structure:7,17,18,23-26,28,31,33,36,39,40;Inversion:21"
##   (one argument, without the line break)
##
## From R:
##   source("plot-loadings.R")
##   d <- read_loadings("pcaone.loadings", "pcaone.mbim", pcs = 1:10)
##   d <- subset_variants(d, d$snp %in% c("rs123", "rs456"))
##   plot_loadings(d, highlight = 1:2, groups = list(HLA = c(6, 12)))
##   loading_peaks(d)
##
## Needs R >= 4.0 and the data.table package.

## Categorical colours, assigned in this order: groups first, then
## highlighted PCs.
PL_COLS <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
             "#e87ba4", "#008300", "#4a3aa7", "#e34948")
PL_GREY <- "grey62"   # lines of PCs that are neither highlighted nor grouped
PL_INK  <- "grey20"   # text and axes

## ---- parsers ---------------------------------------------------------------

## "1-3,X" -> c("1", "2", "3", "X"). Numeric ranges are expanded.
parse_set <- function(s) {
  if (is.null(s) || !length(s)) return(NULL)
  if (!is.character(s)) return(as.character(s))
  tok <- strsplit(gsub("\\s+", "", paste(s, collapse = ",")), ",", fixed = TRUE)[[1]]
  tok <- tok[nzchar(tok)]
  out <- unlist(lapply(tok, function(x) {
    ab <- strsplit(x, "-", fixed = TRUE)[[1]]
    if (length(ab) == 2L && all(grepl("^[0-9]+$", ab))) as.character(seq(as.integer(ab[1]), as.integer(ab[2])))
    else x
  }))
  unique(out)
}

parse_int_set <- function(s) {
  if (is.null(s)) return(NULL)
  v <- suppressWarnings(as.integer(parse_set(s)))
  if (anyNA(v)) stop("cannot parse '", paste(s, collapse = ","), "' as a list of integers like 1-5,8")
  v
}

## "25M", "250k", "2.5e7", "25Mb" -> bp
parse_bp <- function(x) {
  x <- sub("B$", "", toupper(trimws(x)))
  mult <- c(K = 1e3, M = 1e6, G = 1e9)
  u <- substring(x, nchar(x))
  v <- if (u %in% names(mult)) as.numeric(substr(x, 1L, nchar(x) - 1L)) * mult[[u]] else as.numeric(x)
  if (is.na(v)) stop("cannot parse position '", x, "'")
  v
}

## "6:25M-35M", "chr6:25000000-35000000" or just "6"
parse_region <- function(s) {
  if (is.null(s)) return(NULL)
  if (is.list(s)) return(s)
  p <- strsplit(sub("^chr", "", s, ignore.case = TRUE), ":", fixed = TRUE)[[1]]
  if (length(p) == 1L) return(list(chr = p, from = -Inf, to = Inf))
  r <- strsplit(p[2], "-", fixed = TRUE)[[1]]
  if (length(p) != 2L || length(r) != 2L) stop("region must look like 6:25M-35M, got '", s, "'")
  list(chr = p[1], from = parse_bp(r[1]), to = parse_bp(r[2]))
}

## PC groups, as a named list (label -> PCs), a file, or an inline spec.
## Inline: "Population structure:1,2,4;HLA:6,12,15-16".
## File: one "label: PCs" per line, or two columns "<PC> <label>".
parse_groups <- function(spec) {
  if (is.null(spec)) return(NULL)
  if (is.list(spec)) return(lapply(spec, as.integer))
  lines <- if (length(spec) == 1L && file.exists(spec)) readLines(spec, warn = FALSE)
           else strsplit(paste(spec, collapse = ";"), ";", fixed = TRUE)[[1]]
  lines <- trimws(sub("#.*", "", lines))
  lines <- lines[nzchar(lines)]
  g <- list()
  for (ln in lines) {
    if (grepl(":", ln, fixed = TRUE)) {
      lab <- trimws(sub(":.*", "", ln))
      pcs <- parse_int_set(sub("^[^:]*:", "", ln))
    } else {
      lab <- trimws(sub("^\\S+\\s+", "", ln))
      pcs <- parse_int_set(sub("\\s.*", "", ln))
    }
    if (!nzchar(lab)) stop("group '", ln, "' has no name")
    g[[lab]] <- unique(c(g[[lab]], pcs))
  }
  dup <- unlist(g)[duplicated(unlist(g))]
  if (length(dup)) warning("PC ", paste(unique(dup), collapse = ","), " is in more than one group; the first one is used")
  g
}

## group label of each PC, NA if none
pc_group <- function(pcs, groups) {
  if (!length(groups)) return(rep(NA_character_, length(pcs)))
  lab <- rep(names(groups), lengths(groups))
  lab[match(pcs, unlist(groups))]
}

## ---- data ------------------------------------------------------------------

#' Read a PCAone .loadings file and the matching .mbim/.bim.
#'
#' @param loadings .loadings file (no header, one column per PC)
#' @param bim      .mbim or .bim with the same rows; NULL plots by variant index
#' @param pcs      PCs to read (default all)
#' @param chr      keep only these chromosomes, e.g. c(1:22, "X")
#' @param region   keep one region, "6:25M-35M" or list(chr, from, to)
#' @param sigvals  .sigvals file (default: matching .loadings prefix)
#' @param weighted multiply each PC by its singular value (default TRUE)
read_loadings <- function(loadings, bim = NULL, pcs = NULL, chr = NULL, region = NULL,
                          sigvals = NULL, weighted = TRUE) {
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("plot-loadings.R needs the data.table package: install.packages('data.table')")
  nc <- ncol(data.table::fread(loadings, nrows = 1L, header = FALSE))
  pcs <- if (is.null(pcs)) seq_len(nc) else sort(unique(as.integer(pcs)))
  bad <- pcs[pcs < 1L | pcs > nc]
  if (length(bad)) stop(sprintf("%s has %d PCs; cannot plot PC %s", loadings, nc, paste(bad, collapse = ",")))
  V <- data.table::fread(loadings, header = FALSE, select = pcs, showProgress = FALSE)
  data.table::setnames(V, paste0("PC", pcs))
  n <- nrow(V)
  weights <- rep(1, length(pcs))
  if (weighted) {
    if (is.null(sigvals)) sigvals <- paste0(sub("\\.loadings$", "", loadings), ".sigvals")
    if (!file.exists(sigvals))
      stop("cannot find ", sigvals, "; give --sigvals or use --no-weights (weighted = FALSE in R)")
    ## PCAone's #nsamples,nsnps[,key=value] header is a comment to scan().
    values <- scan(sigvals, what = double(), comment.char = "#", quiet = TRUE)
    if (length(values) != nc)
      stop(sprintf("%s has %d singular values but %s has %d PCs", sigvals, length(values), loadings, nc))
    if (any(!is.finite(values) | values < 0))
      stop("singular values must be finite and non-negative in ", sigvals)
    weights <- values[pcs]
    for (j in seq_along(pcs)) data.table::set(V, j = j, value = V[[j]] * weights[j])
  }

  chrv <- snp <- bp <- NULL
  if (!is.null(bim)) {
    map <- data.table::fread(bim, header = FALSE, select = c(1L, 2L, 4L), showProgress = FALSE)
    if (nrow(map) != n)
      stop(sprintf("%s has %d rows but %s has %d; they must describe the same variants",
                   bim, nrow(map), loadings, n))
    chrv <- sub("^chr", "", as.character(map[[1]]), ignore.case = TRUE)
    snp <- as.character(map[[2]])
    bp <- as.numeric(map[[3]])
    rm(map)
  } else {
    chrv <- rep("", n)
  }

  region <- parse_region(region)
  chr <- parse_set(chr)
  if ((length(chr) || !is.null(region)) && is.null(bim)) stop("--chr/--region need the .mbim/.bim file")
  keep <- NULL
  if (length(chr)) keep <- chrv %in% sub("^chr", "", chr, ignore.case = TRUE)
  if (!is.null(region)) {
    k <- chrv == region$chr & bp >= region$from & bp <= region$to
    keep <- if (is.null(keep)) k else keep & k
  }
  row <- if (is.null(keep)) seq_len(n) else which(keep)
  if (!length(row)) stop("no variants left after --chr/--region")
  if (length(row) < n) {
    V <- V[row]
    chrv <- chrv[row]
    snp <- snp[row]
    bp <- bp[row]
  }
  structure(list(V = V, chr = chrv, snp = snp, bp = bp, row = row, pcs = pcs, file = loadings,
                 weighted = weighted, weights = weights, sigvals = if (weighted) sigvals else NULL),
            class = "pcaone_loadings")
}

#' Subset variants while keeping loadings and variant metadata aligned.
#'
#' @param d        A pcaone_loadings object returned by read_loadings().
#' @param keep     A logical mask of length nrow(d$V), or integer row indices.
#' @return A pcaone_loadings object retaining the original file row numbers.
#' @export
subset_variants <- function(d, keep) {
  if (is.logical(keep)) {
    stopifnot(length(keep) == nrow(d$V), !anyNA(keep))
    keep <- which(keep)
  }
  stopifnot(
    is.numeric(keep), length(keep) > 0L, !anyNA(keep),
    all(keep == floor(keep)),
    all(keep >= 1L & keep <= nrow(d$V))
  )

  d$V <- d$V[keep, ]
  for (field in c("chr", "snp", "bp", "row")) {
    if (!is.null(d[[field]])) d[[field]] <- d[[field]][keep]
  }
  d
}

#' Max |loading| per bin of `window` consecutive variants (per chromosome).
#' `window = NULL` picks it so that there are about `target` bins.
bin_loadings <- function(d, window = NULL, target = 4000L, xaxis = c("index", "bp")) {
  xaxis <- match.arg(xaxis)
  if (xaxis == "bp" && is.null(d$bp)) stop("an x axis in bp needs the .mbim/.bim file")
  n <- length(d$chr)
  if (is.null(window) || is.na(window)) window <- max(1L, ceiling(n / target))
  window <- as.integer(window)

  r <- rle(d$chr)
  if (anyDuplicated(r$values))
    warning("variants are not sorted by chromosome; a split chromosome is labelled once per block")
  len <- r$lengths
  nrun <- length(len)
  run <- rep.int(seq_len(nrun), len)
  nb <- (len - 1L) %/% window + 1L
  bin <- c(0L, cumsum(nb))[run] + (sequence(len) - 1L) %/% window + 1L
  last <- cumsum(len)
  first <- last - len + 1L

  if (xaxis == "bp") {
    rng <- data.table::data.table(run = run, bp = d$bp)[, list(lo = min(bp), hi = max(bp)), by = run]
    span <- rng$hi - rng$lo
    gap <- if (nrun > 1L) 0.004 * sum(span) else 0
    ## lay chromosomes end to end; the first keeps its own coordinates
    off <- c(0, cumsum(span + gap))[seq_len(nrun)] - rng$lo + rng$lo[1]
    x <- d$bp + off[run]
    bounds <- data.frame(chr = r$values, start = rng$lo + off, end = rng$hi + off)
  } else {
    x <- seq_len(n)
    bounds <- data.frame(chr = r$values, start = first - 0.5, end = last + 0.5)
  }
  bounds$mid <- (bounds$start + bounds$end) / 2

  xb <- data.table::data.table(b = bin, x = x)[, list(lo = min(x), hi = max(x)), by = b]
  M <- matrix(NA_real_, nrow(xb), length(d$pcs), dimnames = list(NULL, names(d$V)))
  for (j in seq_along(d$pcs))
    M[, j] <- data.table::data.table(b = bin, v = abs(d$V[[j]]))[, list(m = max(v)), by = b]$m
  list(x = (xb$lo + xb$hi) / 2, run = run[!duplicated(bin)],
       M = M, bounds = bounds, window = window, xaxis = xaxis, pcs = d$pcs, weighted = isTRUE(d$weighted))
}

#' The top variant of each PC, with the share of the PC's sum of squared
#' loadings that falls on that variant's chromosome (high = one region drives
#' the PC, e.g. an inversion or the HLA).
loading_peaks <- function(d, groups = NULL) {
  groups <- parse_groups(groups)
  multi <- length(unique(d$chr)) > 1L
  res <- lapply(seq_along(d$pcs), function(j) {
    v <- d$V[[j]]
    i <- which.max(abs(v))
    ss <- v * v
    data.frame(PC = d$pcs[j],
               chr = d$chr[i],
               snp = if (is.null(d$snp)) NA_character_ else d$snp[i],
               bp = if (is.null(d$bp)) NA_real_ else d$bp[i],
               row = d$row[i],
               loading = signif(v[i], 4),
               chr_share = if (multi) round(sum(ss[d$chr == d$chr[i]]) / sum(ss), 3) else NA_real_,
               stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, res)
  res$group <- pc_group(res$PC, groups)
  res
}

## ---- drawing helpers -------------------------------------------------------

.ink_on <- function(col) {
  rgb <- grDevices::col2rgb(col) / 255
  ifelse(colSums(rgb * c(0.299, 0.587, 0.114)) > 0.6, "black", "white")
}

.lighter <- function(col, f = 0.45) {
  rgb <- grDevices::col2rgb(col) / 255
  grDevices::rgb(t(rgb + (1 - rgb) * f))
}

.peak_label <- function(pk) {
  if (is.na(pk$bp)) return(sprintf("variant %s", format(pk$row, big.mark = ",")))
  sprintf("%s%s:%.2f Mb", if (grepl("^[0-9XYM]", pk$chr)) "chr" else "", pk$chr, pk$bp / 1e6)
}

.shade_chr <- function(b, ylim) {
  bd <- b$bounds
  if (nrow(bd) < 2L) return(invisible())
  even <- seq_len(nrow(bd)) %% 2L == 0L
  rect(bd$start[even], ylim[1], bd$end[even], ylim[2], col = "grey94", border = NA)
}

.x_axis <- function(b, cex, outer = FALSE) {
  bd <- b$bounds
  ax <- list(side = 1, col = "grey45", col.axis = PL_INK, cex.axis = 0.85 * cex)
  line <- if (nrow(bd) > 1L) 1.3 else 1.9   # chromosome names sit closer than ticks
  lab <- function(txt) if (outer) mtext(txt, side = 1, line = line, outer = TRUE, col = PL_INK, cex = 0.9 * cex)
                       else title(xlab = txt, col.lab = PL_INK, line = line)
  if (nrow(bd) == 1L) {
    if (b$xaxis == "bp") {
      at <- pretty(c(bd$start, bd$end))
      at <- at[at >= bd$start & at <= bd$end]
      do.call(axis, c(ax, list(at = at, labels = format(at / 1e6))))
      lab(sprintf("%s%s position (Mb)", if (nzchar(bd$chr)) "chr" else "", bd$chr))
    } else {
      at <- pretty(c(bd$start, bd$end))
      at <- at[at >= bd$start & at <= bd$start + 0.97 * (bd$end - bd$start)]
      do.call(axis, c(ax, list(at = at, labels = format(at, big.mark = ",", scientific = FALSE, trim = TRUE))))
      lab(if (nzchar(bd$chr)) sprintf("Variants on chr%s", bd$chr) else "Variant index")
    }
  } else {
    do.call(axis, c(ax, list(at = bd$mid, labels = bd$chr, tick = FALSE, line = -0.7, gap.axis = 0.3)))
    lab("Chromosome")
  }
}

## loadings are small (~1/sqrt(M)); show them as multiples of 10^e
.yscale <- function(M) {
  e <- floor(log10(max(M, na.rm = TRUE)))
  if (!is.finite(e) || e > -2) e <- 0
  e
}

.ylab <- function(b, e) {
  w <- format(b$window, big.mark = ",")
  quantity <- if (isTRUE(b$weighted)) "|loading x singular value|" else "|loading|"
  what <- if (b$window > 1L) paste0("max ", quantity, " per ", w, " SNPs") else quantity
  if (e == 0) what else bquote(.(what) ~ "(" * symbol("\264") * 10^.(e) * ")")
}

## fewest legend rows that fit `avail` inches (legend fills column by column)
.legend_layout <- function(lab, cex, avail) {
  n <- length(lab)
  if (!n) return(list(nrow = 0L, ncol = 0L))
  w <- strwidth(lab, units = "inches", cex = cex)
  pad <- strwidth("MMMM", units = "inches", cex = cex)
  for (nr in seq_len(n)) {
    nc <- ceiling(n / nr)
    colw <- vapply(split(w, rep(seq_len(nc), each = nr, length.out = n)), max, 0)
    if (sum(colw + pad) <= avail) break
  }
  list(nrow = ceiling(n / nc), ncol = nc, width = colw + pad - strwidth("MM", units = "inches", cex = cex))
}

## ---- main plot -------------------------------------------------------------

#' Plot loadings.
#'
#' @param d         from read_loadings()
#' @param mode      "overlay": all PCs in one panel, with a numbered marker at
#'                  each PC's peak. "panel": one row per PC, chromosomes in
#'                  alternating shades, peak position labelled.
#' @param highlight PCs drawn in colour in overlay mode (others are grey)
#' @param groups    PC categories (see parse_groups); colour the peak markers
#' @param window    variants per bin; NULL gives ~`target` bins
#' @param xaxis     "index" (variant order) or "bp"; NULL: bp for one chromosome
#' @param peaks     mark each PC's peak
#' @param fixed_y   panel mode: the same y axis for every PC
#' @return the loading_peaks() table, invisibly
plot_loadings <- function(d, mode = c("overlay", "panel"), highlight = NULL, groups = NULL,
                          window = NULL, target = 4000L, xaxis = NULL, peaks = TRUE,
                          fixed_y = FALSE, title = NULL, cex = 1) {
  mode <- match.arg(mode)
  groups <- parse_groups(groups)
  if (is.null(xaxis)) xaxis <- if (!is.null(d$bp) && length(unique(d$chr)) == 1L) "bp" else "index"
  b <- bin_loadings(d, window = window, target = target, xaxis = xaxis)
  b$e <- .yscale(b$M)
  b$M <- b$M / 10^b$e
  pk <- loading_peaks(d, groups)

  highlight <- parse_int_set(highlight)
  if (mode == "panel") highlight <- NULL
  miss <- setdiff(highlight, d$pcs)
  if (length(miss)) warning("PC ", paste(miss, collapse = ","), " not plotted; not highlighted")
  highlight <- intersect(highlight, d$pcs)

  ## colours: one slot per group (in the order given, whether plotted or not),
  ## then one per highlighted PC
  nslot <- length(groups) + length(highlight)
  if (nslot > length(PL_COLS))
    stop(sprintf("at most %d groups + highlighted PCs can get their own colour (got %d)", length(PL_COLS), nslot))
  gcol <- stats::setNames(PL_COLS[seq_along(groups)], names(groups))
  hcol <- PL_COLS[length(groups) + seq_along(highlight)]

  if (mode == "overlay") .plot_overlay(b, pk, gcol, highlight, hcol, peaks, title, cex)
  else .plot_panels(b, pk, gcol, peaks, fixed_y, title, cex)
  invisible(pk)
}

.plot_overlay <- function(b, pk, gcol, hl, hcol, peaks, title, cex) {
  M <- b$M
  x <- b$x
  pcs <- b$pcs
  grp <- pk$group
  used <- names(gcol) %in% grp

  ## legend entries: groups present, "Other" peaks, highlighted PCs, other PCs
  leg <- list(lab = character(), pch = integer(), bg = character(), col = character(), lty = integer(), lwd = numeric())
  add <- function(lab, pch = NA, bg = NA, col, lty = NA, lwd = NA) {
    leg$lab <<- c(leg$lab, lab); leg$pch <<- c(leg$pch, pch); leg$bg <<- c(leg$bg, bg)
    leg$col <<- c(leg$col, col); leg$lty <<- c(leg$lty, lty); leg$lwd <<- c(leg$lwd, lwd)
  }
  if (peaks && any(used)) {
    for (g in names(gcol)[used]) add(g, pch = 21, bg = gcol[[g]], col = "white")
    if (anyNA(grp)) add("Ungrouped", pch = 21, bg = "grey80", col = "white")
  }
  for (k in seq_along(hl)) add(paste0("PC", hl[k]), col = hcol[k], lty = 1, lwd = 2)
  if (length(hl) && length(hl) < length(pcs)) add("Other PCs", col = PL_GREY, lty = 1, lwd = 1.5)

  ## margins: legend rows + title above the plot
  lay <- .legend_layout(leg$lab, 0.9 * cex, par("din")[1] - 5.4 * par("csi"))
  nrow_leg <- lay$nrow
  top <- 0.8 + 1.3 * nrow_leg + if (is.null(title)) 0 else 1.6
  op <- par(mar = c(3.2, 4.4, top, 1), mgp = c(2.4, 0.5, 0), tcl = -0.25, las = 1)
  on.exit(par(op))

  ymax <- max(M, na.rm = TRUE) * 1.08
  xlim <- range(b$bounds$start, b$bounds$end)
  plot.new()
  plot.window(xlim, c(0, ymax), xaxs = "i", yaxs = "i")
  .shade_chr(b, c(0, ymax))
  for (j in which(!pcs %in% hl)) lines(x, M[, j], col = PL_GREY, lwd = 0.8)
  for (k in seq_along(hl)) lines(x, M[, match(hl[k], pcs)], col = hcol[k], lwd = 1.2)
  axis(2, col = "grey45", col.axis = PL_INK, cex.axis = 0.85 * cex)
  box(bty = "l", col = "grey45")
  title(ylab = .ylab(b, b$e), col.lab = PL_INK)
  .x_axis(b, cex)

  if (peaks) {
    w <- apply(M, 2, which.max)
    y <- M[cbind(w, seq_along(w))]
    fill <- ifelse(is.na(grp), "grey80", gcol[grp])
    ishl <- pcs %in% hl
    fill[ishl & is.na(grp)] <- hcol[match(pcs[ishl & is.na(grp)], hl)]
    ## draw the tallest first so lower markers stay readable where they overlap
    o <- order(-y)
    points(x[w][o], y[o], pch = 21, bg = fill[o], col = "white", cex = 2.3 * cex, lwd = 1, xpd = NA)
    text(x[w][o], y[o], pcs[o], col = .ink_on(fill[o]), cex = 0.62 * cex, font = 2, xpd = NA)
  }

  if (length(leg$lab)) {
    usr <- par("usr")
    legend(usr[1], usr[4] + strheight("M", cex = 0.5), xjust = 0, yjust = 0, legend = leg$lab,
           pch = leg$pch, pt.bg = leg$bg, col = leg$col, lty = leg$lty, lwd = leg$lwd,
           ncol = lay$ncol, text.width = lay$width * diff(usr[1:2]) / par("pin")[1], bty = "n", xpd = NA, cex = 0.9 * cex,
           pt.cex = 1.8, seg.len = 1.6, x.intersp = 0.6, text.col = PL_INK)
  }
  if (!is.null(title)) mtext(title, side = 3, line = top - 1.3, adj = 0, font = 2, cex = 1.1 * cex, col = PL_INK)
}

.plot_panels <- function(b, pk, gcol, peaks, fixed_y, title, cex) {
  M <- b$M
  x <- b$x
  np <- ncol(M)
  ## group names go in the right margin, one per panel
  glab <- unique(pk$group[!is.na(pk$group)])
  right <- if (length(glab)) 2.2 + max(strwidth(glab, units = "inches", cex = 0.8 * cex)) / par("csi") else 1
  op <- par(mfrow = c(np, 1), mar = c(0.3, 4.4, 0.3, right),
            oma = c(3.2, 1.4, if (is.null(title)) 0.6 else 2.2, 0),
            mgp = c(2.4, 0.5, 0), tcl = -0.25, las = 1)
  on.exit(par(op))
  xlim <- range(b$bounds$start, b$bounds$end)
  yall <- max(M, na.rm = TRUE) * 1.15
  seg <- split(seq_along(x), b$run)   # bins of each chromosome, drawn as one line
  for (j in seq_len(np)) {
    y <- M[, j]
    g <- pk$group[j]
    base <- if (is.na(g)) "grey30" else gcol[[g]]
    cols <- c(base, .lighter(base))
    ylim <- c(0, if (fixed_y) yall else max(y, na.rm = TRUE) * 1.15)
    plot.new()
    plot.window(xlim, ylim, xaxs = "i", yaxs = "i")
    for (r in seq_along(seg)) {
      i <- seg[[r]]
      if (length(i) == 1L) points(x[i], y[i], pch = 16, cex = 0.3, col = cols[(r - 1L) %% 2L + 1L])
      else lines(x[i], y[i], col = cols[(r - 1L) %% 2L + 1L], lwd = 0.9)
    }
    at <- pretty(ylim, n = 3)
    axis(2, at = at[at < ylim[2]], col = "grey45", col.axis = PL_INK, cex.axis = 0.75 * cex)
    box(bty = "l", col = "grey45")
    mtext(paste0("PC", b$pcs[j]), side = 2, line = 3.2, las = 0, font = 2, cex = 0.75 * cex, col = PL_INK)
    if (!is.na(g)) {
      usr <- par("usr")
      legend(usr[2], mean(usr[3:4]), xjust = 0, yjust = 0.5, legend = g, pch = 15, col = base,
             bty = "n", cex = 0.8 * cex, pt.cex = 1.4, text.col = PL_INK, x.intersp = 0.6, xpd = NA)
    }
    if (peaks) {
      w <- which.max(y)
      points(x[w], y[w], pch = 21, bg = base, col = "white", cex = 1.3 * cex, lwd = 1)
      onright <- x[w] > xlim[1] + 0.8 * diff(xlim)
      text(x[w], y[w], .peak_label(pk[j, ]), pos = if (onright) 2 else 4, offset = 0.5,
           cex = 0.75 * cex, col = PL_INK, xpd = NA)
    }
    if (j == np) {
      op2 <- par(xpd = NA)
      .x_axis(b, cex, outer = TRUE)
      par(op2)
    }
  }
  mtext(.ylab(b, b$e), side = 2, line = 0.2, outer = TRUE, las = 0, col = PL_INK, cex = 0.8 * cex)
  if (!is.null(title)) mtext(title, side = 3, line = 0.6, outer = TRUE, adj = 0, font = 2, cex = 1.1 * cex, col = PL_INK)
}

## ---- command line ----------------------------------------------------------

PL_USAGE <- "Plot PCAone SNP loadings (.loadings) along the genome.

Usage: Rscript plot-loadings.R -p PREFIX [options]
       Rscript plot-loadings.R -l FILE.loadings [-b FILE.mbim] [options]

Input
  -p, --prefix P      reads P.loadings, P.sigvals and P.mbim (or P.bim)
  -l, --loadings F    loadings file (overrides --prefix)
      --sigvals F     singular values [matching loadings prefix.sigvals]
      --no-weights    plot raw loadings instead of loading x singular value
  -b, --bim F         .mbim/.bim with the same variants; without it the x axis
                      is the variant index
      --pcs S         PCs to plot, e.g. 1-10 or 1,3,5 [all]
      --chr S         chromosomes to keep, e.g. 1-22 or 6
      --region R      one region, e.g. 6:25M-35M

Plot
  -o, --out F         .png, .pdf, .svg, .jpg or .tiff [PREFIX.loadings.png]
      --mode M        overlay: all PCs in one panel; panel: one row per PC [overlay]
      --highlight S   overlay: draw these PCs in colour, e.g. 1,2
      --groups G      colour PCs by category, inline \"Pop structure:1,2,4;HLA:6,12\"
                      or a file with lines \"label: PCs\" or \"PC label\"
      --window N      variants per bin (max |loading| per bin) [about 4000 bins]
      --xaxis X       index or bp [bp for one chromosome, else index]
      --fixed-y       panel: same y axis for all PCs
      --no-peaks      do not mark each PC's top variant
      --title T       plot title
      --width W       inches [11]
      --height H      inches [overlay 4.5; panel 0.9 + 1 per PC]
      --res N         png/jpg/tiff resolution [150]
      --cex X         scale all text and markers [1]
      --peaks-out F   also write the peak table to F (tsv)
  -h, --help

The peak table (printed) gives each PC's top variant and chr_share, the share
of the PC's sum of squared loadings on that chromosome: high values mean one
region (an inversion, the HLA, a centromere) drives the PC."

parse_args <- function(args) {
  short <- c(p = "prefix", l = "loadings", b = "bim", o = "out", h = "help")
  flags <- c("help", "no-peaks", "fixed-y", "no-weights")
  known <- c("prefix", "loadings", "sigvals", "bim", "out", "pcs", "chr", "region", "mode", "highlight",
             "groups", "window", "xaxis", "title", "width", "height", "res", "cex", "peaks-out", flags)
  opt <- list()
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    val <- NULL
    if (grepl("^--[^=]+=", a)) {
      key <- sub("^--([^=]+)=.*$", "\\1", a)
      val <- sub("^--[^=]+=", "", a)
    } else if (grepl("^--.", a)) {
      key <- substring(a, 3L)
    } else if (grepl("^-[A-Za-z]$", a) && substring(a, 2L) %in% names(short)) {
      key <- short[[substring(a, 2L)]]
    } else stop("unexpected argument '", a, "' (see -h)", call. = FALSE)
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

open_device <- function(file, width, height, res = 150) {
  ext <- tolower(sub(".*\\.", "", basename(file)))
  switch(ext,
    png = grDevices::png(file, width = width, height = height, units = "in", res = res),
    jpg = , jpeg = grDevices::jpeg(file, width = width, height = height, units = "in", res = res, quality = 95),
    tif = , tiff = grDevices::tiff(file, width = width, height = height, units = "in", res = res, compression = "lzw"),
    pdf = grDevices::pdf(file, width = width, height = height, useDingbats = FALSE),
    svg = grDevices::svg(file, width = width, height = height),
    stop("unknown output format '.", ext, "'; use .png, .pdf, .svg, .jpg or .tiff", call. = FALSE))
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  o <- parse_args(args)
  if (isTRUE(o$help) || !length(args)) {
    cat(PL_USAGE, "\n")
    return(invisible())
  }
  loadings <- o$loadings
  if (is.null(loadings)) {
    if (is.null(o$prefix)) stop("give -p PREFIX or -l FILE.loadings", call. = FALSE)
    loadings <- paste0(o$prefix, ".loadings")
  }
  if (!file.exists(loadings)) stop("cannot find ", loadings, call. = FALSE)
  base <- sub("\\.loadings$", "", loadings)
  bim <- o$bim
  if (is.null(bim)) {
    cand <- paste0(c(o$prefix, base), rep(c(".mbim", ".bim"), each = length(c(o$prefix, base))))
    bim <- cand[file.exists(cand)][1]
    if (is.na(bim)) {
      bim <- NULL
      message("no .mbim/.bim found next to ", loadings, "; plotting by variant index")
    }
  }
  mode <- if (is.null(o$mode)) "overlay" else match.arg(o$mode, c("overlay", "panel"))
  out <- if (is.null(o$out)) paste0(base, ".loadings.png") else o$out
  num <- function(x, def) if (is.null(x)) def else as.numeric(x)

  d <- read_loadings(loadings, bim, pcs = parse_int_set(o$pcs), chr = o$chr, region = o$region,
                     sigvals = o$sigvals, weighted = !isTRUE(o$`no-weights`))
  message(sprintf("read %s variants x %d PCs from %s", format(length(d$row), big.mark = ","),
                  length(d$pcs), loadings))
  width <- num(o$width, 11)
  height <- num(o$height, if (mode == "panel") 0.9 + 1.0 * length(d$pcs) else 4.5)
  open_device(out, width, height, num(o$res, 150))
  pk <- tryCatch(
    plot_loadings(d, mode = mode, highlight = o$highlight, groups = o$groups,
                  window = if (is.null(o$window)) NULL else as.integer(o$window),
                  xaxis = o$xaxis, peaks = !isTRUE(o$`no-peaks`), fixed_y = isTRUE(o$`fixed-y`),
                  title = o$title, cex = num(o$cex, 1)),
    finally = grDevices::dev.off())
  message("wrote ", out)
  print(pk, row.names = FALSE)
  if (!is.null(o$`peaks-out`)) {
    utils::write.table(pk, o$`peaks-out`, sep = "\t", quote = FALSE, row.names = FALSE)
    message("wrote ", o$`peaks-out`)
  }
  invisible(pk)
}

if (sys.nframe() == 0L && !interactive()) {
  tryCatch(main(), error = function(e) {
    message("error: ", conditionMessage(e))
    quit(status = 1)
  })
}
