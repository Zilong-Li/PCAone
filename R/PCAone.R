## R wrapper for PCAone
library(data.table)
PCAone <- function (
                    program = "PCAone",
                    k = 10,
                    bfile = NULL,
                    bgen  = NULL,
                    beagle = NULL,
                    out = tempfile(),
                    arnoldi = FALSE,
                    threads = 1,
                    fast = FALSE,
                    emu = FALSE,
                    pcangsd = FALSE,
                    memory = 0,
                    maxp = 20,
                    bands = 128,
                    verbose = FALSE,
                    tol_em = 1e-5,
                    tol_maf =  1e-4,
                    tol_halko = 1e-4,
                    maxiter = 100,
                    imaxiter = 1000,
                    itol = 1e-6,
                    ncv = 20
                    ) {
  infile <- NULL
  if (!is.null(bfile)) infile <- paste("--bfile", bfile)
  if (!is.null(bgen)) infile <- paste("--bgen", bgen)
  if (!is.null(beagle)) infile <- paste("--beagle", beagle)
  stopifnot(!is.null(infile))

  callCmds <- paste(infile, "-k", k, "-n", threads, "-o", out)
  if (memory > 0) callCmds <- paste(callCmds, "-m", memory)
  if (arnoldi) {
    callCmds <- paste(callCmds, "-a")
  } else {
    callCmds <- paste(callCmds, "--maxp", maxp)
    if (fast) callCmds <- paste(callCmds, "-f")
  }
  if (emu) {
    callCmds <- paste(callCmds, "--emu")
  } else if (pcangsd) {
    callCmds <- paste(callCmds, "--pcangsd")
  }
  if (verbose) callCmds <- paste(callCmds, "-v")
  if (tol_em != 1e-5) callCmds <- paste(callCmds, "--tol-em", tol_em)
  if (tol_maf != 1e-4) callCmds <- paste(callCmds, "--tol-maf", tol_maf)
  if (tol_halko != 1e-4) callCmds <- paste(callCmds, "--tol-halko", tol_halko)
  if (itol != 1e-6) callCmds <- paste(callCmds, "--itol", itol)
  if (maxiter != 100) callCmds <- paste(callCmds, "--maxiter", maxiter)
  if (imaxiter != 1000) callCmds <- paste(callCmds, "--imaxiter", imaxiter)
  if (ncv != 20) callCmds <- paste(callCmds, "--ncv", ncv)

  print(paste("Calling:", program, callCmds))
  ## system2(program, callCmds, stdout = paste0(out, ".o"), stderr = paste0(out, ".e"))
  system2(program, callCmds, stdout = paste0(out, ".log"))
  U <- fread(paste0(out, ".eigvecs"))
  S <- fread(paste0(out, ".eigvals"))
  res <-  list(u = as.matrix(U), d = as.vector(S$V1))
}
