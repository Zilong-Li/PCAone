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

