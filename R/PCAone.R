## R wrapper for PCAone
PCAone <- function (
    program = "PCAone"
    bfile = NULL,
    bgen  = NULL,
    beagle = NULL,
    k = 10,
    threads = 1,
    memory = 0,
    maxp = 20,
    arnoldi = FALSE,
    halko = TRUE,
    fast = FALSE,
    emu = FALSE,
    pcangsd = FALSE,
    out = NULL,
    verbose = FALSE
    ) {

    infile <- NULL
    if (!is.null(bfile)) infile <- paste("--bfile", bfile)
    if (!is.null(bgen)) infile <- paste("--bgen", bgen)
    if (!is.null(beagle)) infile <- paste("--beagle", beagle)
    stopifnot(!is.null(infile))

    if (is.null(out)) out = tempfile()
    callCmds <- paste(infile, "-k", k, "-n", threads, "-o", out)
    if (memory > 0) callCmds <- paste(callCmds, "-m", memory)
    if (arnoldi) {
        callCmds <- paste(callCmds, "-a")
    } else if (halko) {
        callCmds <- paste(callCmds, "--maxp", maxp)
        if (fancy) callCmds <- paste(callCmds, "-f")
    }
    if (emu) {
        callCmds <- paste(callCmds, "--emu")
    } else if (pcangsd) {
        callCmds <- paste(callCmds, "--pcangsd")
    }
    if (verbose) callCmds <- paste(callCmds, "-v")

    print(paste("Calling:", program, callCmds))
    system2(program, callCmds, stdout = "", stderr = "")
}

norm_vec <- function(x) {
    sqrt(sum(x^2))
}

mev <- function(X, Y){
    o = mean(apply(Y, 2, function(x) { norm_vec(t(X) %*% x) }))
    o
}

