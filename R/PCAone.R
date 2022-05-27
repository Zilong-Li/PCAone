## R wrapper for PCAone
library(data.table)
PCAone <- function(program = "PCAone",
                   k = 10,
                   bfile = NULL,
                   bgen = NULL,
                   beagle = NULL,
                   csv = NULL,
                   arnoldi = TRUE,
                   fast = FALSE,
                   halko = FALSE,
                   out = tempfile(),
                   cpmed = FALSE,
                   printv = FALSE,
                   shuffle = FALSE,
                   emu = FALSE,
                   pcangsd = FALSE,
                   threads = 1,
                   memory = 0,
                   maxp = 20,
                   oversamples = 10,
                   bands = 64,
                   tol_em = 1e-4,
                   tol_maf = 1e-4,
                   tol_rsvd = 1e-4,
                   maxiter = 100,
                   imaxiter = 1000,
                   itol = 1e-6,
                   ncv = 20,
                   tmp = NULL,
                   verbose = FALSE) {
  if (file.exists(paste0(out, ".eigvecs")) && file.exists(paste0(out, ".eigvals"))) {
    if (printv) {
      res <- list(u = as.matrix(fread(paste0(out, ".eigvecs"))), d = as.vector(fread(paste0(out, ".eigvals"))[,1]), v = as.matrix(fread(paste0(out, ".loadings"))))
    } else {
      res <- list(u = as.matrix(fread(paste0(out, ".eigvecs"))), d = as.vector(fread(paste0(out, ".eigvals"))[,1]))
    }

    return(res)
  }
  infile <- NULL
  if (!is.null(bfile)) infile <- paste("--bfile", bfile)
  if (!is.null(bgen)) infile <- paste("--bgen", bgen)
  if (!is.null(beagle)) infile <- paste("--beagle", beagle)
  if (!is.null(csv)) infile <- paste("--csv", csv)
  stopifnot(!is.null(infile))

  callCmds <- paste(
    infile, "-k", k, "-n", threads, "-o", out
  )
  if (arnoldi) {
    callCmds <- paste(callCmds, "-a")
  } else if (fast) {
    callCmds <- paste(callCmds, "-f")
  } else if (halko) {
    callCmds <- paste(callCmds, "-h")
  }
  if (emu) {
    callCmds <- paste(callCmds, "--emu")
  } else if (pcangsd) {
    callCmds <- paste(callCmds, "--pcangsd")
  }
  if (!is.null(tmp)) callCmds <- paste(callCmds, "--tmp", tmp)
  if (memory > 0) callCmds <- paste(callCmds, "-m", memory)
  if (ncv != 20) callCmds <- paste(callCmds, "--ncv", ncv)
  if (maxp != 20) callCmds <- paste(callCmds, "--maxp", maxp)
  if (bands != 64) callCmds <- paste(callCmds, "--bands", bands)
  if (tol_rsvd != 1e-4) callCmds <- paste(callCmds, "--tol-rsvd", tol_rsvd)
  if (tol_em != 1e-4) callCmds <- paste(callCmds, "--tol-em", tol_em)
  if (tol_maf != 1e-4) callCmds <- paste(callCmds, "--tol-maf", tol_maf)
  if (itol != 1e-6) callCmds <- paste(callCmds, "--itol", itol)
  if (maxiter != 100) callCmds <- paste(callCmds, "--maxiter", maxiter)
  if (imaxiter != 100) callCmds <- paste(callCmds, "--imaxiter", imaxiter)
  if (oversamples != 10) callCmds <- paste(callCmds, "--oversamples", oversamples)
  if (printv) callCmds <- paste(callCmds, "--printv")
  if (shuffle) callCmds <- paste(callCmds, "--shuffle")
  if (verbose) callCmds <- paste(callCmds, "-v")

  print(paste("Calling:", program, callCmds))
  ## system2(program, callCmds, stdout = paste0(out, ".o"), stderr = paste0(out, ".e"))
  system2(program, callCmds, stdout = paste0(out, ".llog"))
  if (printv) {
    res <- list(u = as.matrix(fread(paste0(out, ".eigvecs"))), d = as.vector(fread(paste0(out, ".eigvals"))[,1]), v = as.matrix(fread(paste0(out, ".loadings"))))
  } else {
    res <- list(u = as.matrix(fread(paste0(out, ".eigvecs"))), d = as.vector(fread(paste0(out, ".eigvals"))[,1]))
  }

  return(res)
}
