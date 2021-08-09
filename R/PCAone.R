## R wrapper for PCAone
PCAone <- function (
    bfile = NULL,
    bgen  = NULL,
    beagle = NULL,
    k = 10,
    threads = 1,
    blocksize = 0,
    powers = 3,
    band =  4,
    arnoldi = FALSE,
    halko = TRUE,
    fancy = FALSE,
    emu = FALSE,
    pcangsd = FALSE,
    out = NULL,
    verbose = FALSE,
    test = FALSE,
    program = "PCAone"
    ) {

    infile <- NULL
    if (!is.null(bfile)) infile <- paste("--bfile", bfile)
    if (!is.null(bgen)) infile <- paste("--bgen", bgen)
    if (!is.null(beagle)) infile <- paste("--beagle", beagle)
    stopifnot(!is.null(infile))

    if (is.null(out)) out = tempfile()
    callCmds <- paste(infile, "-k", k, "-n", threads, "-o", out)
    if (blocksize > 0) callCmds <- paste(callCmds, "-b", blocksize)
    if (arnoldi) {
        callCmds <- paste(callCmds, "-A")
    } else if (halko) {
        callCmds <- paste(callCmds, "-p", powers)
        if (fancy) callCmds <- paste(callCmds, "--fancy", "--band", band)
    }
    if (emu) {
        callCmds <- paste(callCmds, "--emu")
    } else if (pcangsd) {
        callCmds <- paste(callCmds, "--pcangsd")
    }
    if (verbose) callCmds <- paste(callCmds, "-v")
    if (test) callCmds <- paste(callCmds, "--test")

    print(paste("Running:", program, callCmds))
    system2(program, callCmds, stdout = "", stderr = "")
}

norm_vec <- function(x) {
    sqrt(sum(x^2))
}

mev <- function(X, Y){
    o = mean(apply(Y, 2, function(x) { norm_vec(t(X) %*% x) }))
    o
}


test_halko_batch <- function() {
    bfile <- "test.n10k.m10k"
    out_arnoldi <- "out.A"
    out_halko <- "out.H"
    k <-  10
    threads <- 6

    system.time(PCAone(bfile = bfile, k = k, threads = threads, out = out_arnoldi, arnoldi = T))
    evec_A <- "out.A.eigvecs"
    A <- read.table(evec_A)

    system.time(PCAone(bfile = bfile, k = k, threads = threads, out = out_halko, powers = 40))
    evec_H <- "out.H.eigvecs"
    H <- read.table(evec_H)
    mev(A[1:2], H[1:2])
}


test_halko_block <- function() {
    bfile <- "test.n10k.m10k"
    out_arnoldi <- "out.A"
    out_halko <- "out.H"
    k <-  10
    threads <- 6
    powers <- 40

    system.time(PCAone(bfile = bfile, k = k, threads = threads, out = out_arnoldi, arnoldi = T))
    evec_A <- "out.A.eigvecs"
    A <- read.table(evec_A)

    system.time(PCAone(bfile = bfile, k = k, threads = threads, out = out_halko, powers = powers))
    evec_H <- "out.H.eigvecs"
    H1 <- read.table(evec_H)
    print(mev(A[1:5], H1[1:5]))

    system.time(PCAone(bfile = bfile, k = k, threads = threads, out = out_halko, powers = powers, blocksize = 5000))
    evec_H <- "out.H.eigvecs"
    H2 <- read.table(evec_H)
    print(mev(H1[1:5], H2[1:5]))

    system.time(PCAone(bfile = bfile, k = k, threads = threads, out = out_halko, powers = powers, fancy = TRUE, blocksize = 500))
    evec_H <- "out.H.eigvecs"
    H3 <- read.table(evec_H)
    print(mev(A[1:5], H3[1:5]))
}

