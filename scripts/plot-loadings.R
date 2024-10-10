library(data.table)
library(caTools)

mycols <- c("#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#BEBADA", "#FFED6F", "#8DD3C7")
palette(mycols)

df <- fread("pcaone.loadings", h = F, data.table = F)
npc <- ncol(df)
bim <- fread("plink.bim", h = F, data.table = F)
bim <- bim[, c(1, 2, 4)]
stopifnot(nrow(bim) == nrow(df))
df <- cbind(bim, df)
colnames(df) <- c("chr", "snp", "bp", paste0("pc", 1:npc))
chrtable <- table(df$chr)

window <- 1000
mat <- apply(df[, -(1:3)], 2, function(x) abs(caTools::runmax(x, k = window)))
groups <- list(
  "Population structure" = c(1, 2, 4),
  "Centromere" = c(3, 5, 8, 9, 10, 11, 13, 14, 19, 20, 22, 29, 32, 34, 35, 37),
  "HLA" = c(6, 12, 15, 16, 27, 30, 38),
  "Other structure" = c(7, 17, 18, 23:26, 28, 31, 33, 36, 39, 40),
  "Inversion" = c(21)
)

# anno : a vecotor of pc to be colored
plotLegend <- function(anno = c(1, 2)) {
  par(mar = c(0, 0, 0, 0), mgp = c(2, 0.6, 0))
  plot(0, 0, col = "transparent", axes = F, xlab = "", ylab = "")
  cats <- names(groups)
  letext <- c(cats, paste0("PC", anno))
  legend("bottomleft",
    legend = letext, col = 1:length(letext), pch = c(rep(16, length(cats)), rep(NA, 3)),
    lty = c(rep(NA, length(cats)), rep(1, 3)), adj = 0, lwd = 4, pt.cex = 3.0, cex = 1.2, text.font = 2,
    horiz = T, xjust = 0, yjust = 0, xpd = T, inset = c(0, 0), bty = "n", x.intersp = 0.1
  )
}


plotLoadings <- function(mat, chrtable, groups, anno = c(1, 2), title = NULL) {
  if (!is.list(groups)) stop("groups must be a list")
  cats <- names(groups)
  ymax <- c(0, max(mat)) * 1.05
  nchr <- length(chrtable)
  par(mar = c(3.5, 3.5, 4, 1), mgp = c(2, 0.6, 0))

  plot(1:nrow(mat), mat[, 1],
    type = "l", col = "gray60", lwd = 2, cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5,
    ylim = ymax, axes = F, xlab = "Chromosomes", ylab = "SNP loadings", main = title
  )
  axis(2, cex.lab = 1.5, cex.axis = 1.5)
  ## text(cumsum(as.vector(chrtable)) - chrtable/2, rep(-0.001,22), 1:22, pch = "|", xpd = T, ...)
  text(cumsum(as.vector(chrtable)) - chrtable / 2, rep(-ymax[2] / 30, nchr), 1:nchr, pch = "|", xpd = T, cex = 1.1)
  apply(mat, 2, function(x) lines(1:nrow(mat), x, lwd = 2, col = "gray60"))
  if (is.vector(anno)) {
    sapply(1:length(anno), function(i) lines(1:nrow(mat), mat[, anno[i]], lwd = 2, col = i + length(cats)))
  }
  abline(v = cumsum(as.vector(chrtable)), col = "white", lwd = 3)
  w <- apply(mat, 2, which.max)
  for (i in 1:ncol(mat)) {
    c <- which(cats == names(which(sapply(groups, function(x) i %in% x))))
    points(w[i], mat[w[i], i], col = c, pch = 16, cex = 3)
    text(w[i], mat[w[i], i], i)
  }
}

plotLoadings(mat, chrtable, groups)
