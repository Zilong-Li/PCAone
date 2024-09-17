#!/usr/bin/env Rscript

## this script calculates the adjusted LD (r2).
## the input is the binary file with suffix .residuals outputted by PCAone

args <- commandArgs(trailingOnly = T)

library(Rcpp)
library(RcppEigen)

sourceCpp("misc/ld-scores.cpp")

cppFunction('List ld_curve_cpp(Eigen::MatrixXd m, IntegerVector pos, IntegerVector bins){
int nsnps = m.cols(), nb = bins.size(), winsize = max(bins);
NumericVector r_sum(nb), r2_sum(nb);
IntegerVector counts(nb);
const int df = m.rows() - 1; // N-1
Eigen::VectorXd sds = (m.array().square().colwise().sum() / df).sqrt();
int i, j, b, dist;
for(i = 0; i < nsnps; i++) {
  b = 0;
  for(j = i + 1; j < nsnps; j++) {
    dist = pos(j) - pos(i);
    if(dist > winsize) break;
    if(dist > bins[b]) b = b + 1;
    double r  = (m.col(j).array() * m.col(i).array() / (sds(j) * sds(i))).sum() / df;
    r_sum[b] += r;
    r2_sum[b] += r * r;
    counts[b] += 1;
  }
}
return List::create(Named("bin") = bins, Named("counts") = counts,
                    Named("r_sum") = r_sum, Named("r2_sum") = r2_sum);
}', depends = c("RcppEigen"), rebuild = TRUE)

resfile <- "pcaone.residuals"

bimfile <- "example/plink.bim"

bim <- data.table::fread(bimfile, h =F, data.table = F)

fed(bim)

chr1 <- bim[bim[,1]==1,]

con <- file(resfile, "rb")
dims <- readBin(con, integer(), n = 2)
nsnps <- dims[1]
nsamples <- dims[2]

snps <- nsnps / 10
snps <- nrow(chr1)
m <- matrix(readBin(con, numeric(), n = nsamples * snps, size = 4), ncol = snps)

pos <- chr1[,4]

bins <- c(1000, 5000, 10000, 50000, 100000, 500000, 1e6, 5e6)

system.time(d <- ld_curve_cpp(m, pos, bins))

str(d)
plot(o[,2]/o[,1])

dev.off()
gc()


close(con)

bytes_per_snp <- nsamples * 4

a <- readBin(con, numeric(), n = nsamples, size = 4)
b <- readBin(con, numeric(), n = nsamples, size = 4)
o <- cor(a,b, method = "pearson")**2


n <- 100 

for(i in 1:n){
  b <- readBin(con, numeric(), n = nsamples, size = 4)
  o <- c(o, cor(a,b, method = "pearson")**2)
}

cummean <- function(x) cumsum(x)/1:length(x)

plot(o)

dev.off()

