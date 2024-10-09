
library(data.table)

## sample size correction
bs <- function(x,n) x-(1-x)/(n-2)

parse_ld_bins <- function(fn, n=NULL){
  nbin <- 500
  all <- read.table(fn, head=T)
  chrs <- unique(all$chr)
  num <- sapply(chrs, function(c) all[all$chr==c, "num"][1:nbin])

  weight <- num / rowSums(num)

  r2chr <- sapply(chrs, function(c) all[all$chr==c,"avg_R2"][1:nbin])
  r2all <- rowSums(r2chr*weight)

  dist <- all[all$chr==1,"distance"][1:nbin]
  r2 <- r2all
  if(!is.null(n)) r2 <- bs(r2all,n)
  dat <- cbind(dist=dist, r2=r2)
  as.data.frame(dat)
}

plot_ld_curve <- function(dat, ...){
  ymax <- min(c(max(dat$r2), 0.4))
  ymax <- 0.3
  bins <- c("5Kb"=5e3, "10Kb"=1e4, "50Kb"=5e4, "100Kb"=1e5, "500Kb"=5e5, "1Mb"=1e6, "5Mb"=5e6)
  plot(dat$dist, dat$r2, ylim = c(0, ymax), type = "l", log = "x", axes = F, xlab = "", ylab="", ...)
  axis(1, at=bins, labels = NA)
  mtext(names(bins),side=1,line=1:2*0.9,at=bins)
  ylabr2 <- expression(paste(mean~tilde(r)^2))
  mtext(ylabr2,2,line=2,cex=2)
  mtext("Distance",side=1,line=4,cex=2)
  axis(2,at=0:(as.integer(ymax*10))/10)
}

## https://logfc.wordpress.com/2017/03/15/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/
addlabel <- function(label, adj=0, padj=0, font=1){
  par(xpd=NA)
  di <- dev.size("in")
  x <- grconvertX(c(0, di[1]), from="in", to="user")
  y <- grconvertY(c(0, di[2]), from="in", to="user")

  fig <- par("fig")
  x <- x[1] + (x[2] - x[1]) * fig[1:2]
  y <- y[1] + (y[2] - y[1]) * fig[3:4]

  ax <- x[1] + strwidth(label, cex=3) / 2
  ay <- y[2] - strheight(label, cex=3) / 2
  text(ax-adj, ay-padj, label, cex=3, font=font)
}


