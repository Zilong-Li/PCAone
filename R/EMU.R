#!/usr/bin/env Rscript

library(RSpectra)
library(Rcpp)
library(genio)
sourceCpp("./mat.cpp")
## dyn.load('mat.so')

# control the threads used by high-performance parallelized matrix library
library(RhpcBLASctl)
blas_set_num_threads(10)
######### function #########

read_bed <- function(path) {
	mat <- read_plink(path)$X
}

# G is m*n
estimated_F <- function(G) {
	f <- apply(G, 1, function(v){ sum(v, na.rm=T) / (2 * sum(!is.na(v))) })
}

init_E <- function(G, f) {
	E <- G - 2 * f
	E[is.na(G)] <- 0
	E
}

rsem <- function(V1, V2) {
	m = dim(V1)[1]
	k = dim(V1)[2]
	res = sqrt(sum((V2 - V1)^2) / (k * m))
	res
}

fa = function(x, args) {
    ## as.numeric(args %*% x)
    ## .Call(eigenMapMatMult, args ,x)
    eigenMapMatMult(args ,x)
}

fatr = function(x, args) {
    ## as.numeric(crossprod(args, x))
    eigenCrossProd(args, x)
}

# EMU for diploid coding
EMU <- function(filepath, k, maxiter=100, tole=5e-7) {
	# diploid mode
	G <- read_bed(filepath)
	m <- dim(G)[1]
	n <- dim(G)[2]
	f <- estimated_F(G)
	stopifnot(m == length(f))
	print(sprintf(fmt = "extracting data done and begin to inital E -- %s", Sys.time()))

	# init E
	E <- init_E(G, f)
	print(sprintf("inital E done -- %s", Sys.time()))

	# do SVD
	## s <- svds(E, k)
	print(sprintf("begin to do svds -- %s", Sys.time()))
    s <- svds(fa, k, Atrans = fatr, args = E, dim = dim(E))
	u1 <- s$u
	d1 <- s$d
	v1 <- s$v
    print(d1)
	# here G is turned into {0,1}, in which 0 means non-missing genotype and 1 means missing genotype
	G <- is.na(G) + 0
	eigenUpdateE(E, G, u1, v1, d1, f)
	# do interation
	print(sprintf("begin to do iteration -- %s", Sys.time()))
	for (i in 1:maxiter) {
		print(sprintf("activate the interation (%s)-- %s", i, Sys.time()))
		s <- svds(E, k)
		u2 <- s$u
		d2 <- s$d
		v2 <- s$v
		print(sprintf("SVD interation (%s)-- %s", i, Sys.time()))
		eigenUpdateE(E, G, u2, v2, d2, f)
		print(sprintf("eigenUpdateE interation (%s)-- %s", i, Sys.time()))
		diff <- rsem(v1, v2)
		print(sprintf("Individual allele frequencies estimated (%s), RMSE=%s. -- %s", i, diff, Sys.time()))
		if (diff < tole) {
			print("come to converged!")
			break
		}
		v1 <- v2
		print(sprintf("done the interation (%s)-- %s", i, Sys.time()))
	}

	# standardize matrix
	print(sprintf("standardize the matrix -- %s", Sys.time()))
	E <- sweep(E, 1, sqrt(f*(1-f)), FUN="/")

	# final SVD on standardized matrix
	s <- svds(E, k)
	print(sprintf("final SVD done -- %s", Sys.time()))
	s
}


######## main #########

args <- commandArgs(T)
bedfile <- args[1]
k <- as.numeric(args[2])
outpre <- args[3]

s <- EMU(bedfile, k)
m_snps  <- dim(s$u)[1]
eigval <- (s$d)^2 / m_snps
eigvec <- s$v

outval <- paste0(outpre, ".eigval")
outvec <- paste0(outpre, ".eigvec")
write.table(as.matrix(eigval), outval, sep="\t", row.names=F, col.names=F)
write.table(eigvec, outvec, sep="\t", row.names=F, col.names=F)

outplot <- paste0(outpre, ".png")
png(outplot)
plot(eigvec[,1:2], main="EMU-R", xlab="PC1", ylab="PC2")
dev.off()

