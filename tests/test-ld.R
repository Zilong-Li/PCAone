library(genio)
library(Rcpp)

sourceCpp("tests/test-cor.cpp")

## n <- 10000
## x <- runif(n)
## y <- runif(n)
## cor(x, y)
## m <- as.matrix(cbind(x, y))
## rcppeigen_cor(m)

sourceCpp("tests/test-cor.cpp")
args <- commandArgs(trailingOnly = TRUE)

bfile <- args[1]
# read annotation tables
bim <- read_bim(paste0(bfile, ".bim"))
fam <- read_bim(paste0(bfile, ".fam"))
# read an existing Plink *.bim file
# pass locus and individual IDs as vectors, setting data dimensions too
G <- read_bed(paste0(bfile, ".bed"), bim$id, fam$id)
G <- sweep(G, 1, rowMeans(G), "-")
keep <- prune(t(G), bim$pos, 1000000, 0.2)
out <- bim[as.logical(keep),]
write.table(out, paste0(bfile, ".ld.prune.in"), quote=F, row.name=F, col.names=F, sep="\t")
saveRDS(list(keep, bim), paste0(bfile, ".rds"))

bin <- "pcaone.ld.chr1"
con <- file(bin, "rb")
n <- readBin(con, integer(), n=1)

m <- array(0, dim=c(n, n))
for(i in 1:n){
  c = n - i + 1
  m[i,i:n] <- readBin(con, numeric(), n=c, size = 8)
}

S <- file.size(bin) # file size
S <- S - 4
L <- n * (n + 1.0) / 2.0
U <- S / L
m[upper.tri(m, 1)]

