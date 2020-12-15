
#' @describeIn QUBIC Performs a QUalitative BIClustering.
#'
#' @usage qubic(i, R = FALSE, F = FALSE, d = FALSE, f = 0.85, k = 13, c = 0.90, o = 5000)
#'
#' @importFrom Rcpp evalCpp


qubic <- function(i = NULL, N = FALSE, R = FALSE, Fa = FALSE, d = FALSE, D = FALSE, C = FALSE, n = FALSE, q = 0.05, f = 0.85, k = 13, c = 0.9, o = 100) {
    vec <- c("./qubic", "-i", i)
    if (N) 
        vec <- c(vec, "-N")
    if (R) 
        vec <- c(vec, "-R")
    if (Fa) 
        vec <- c(vec, "-F")  # only do discretization
    if (d) 
        vec <- c(vec, "-d")  # input binrary e.g. chars
    if (D) 
        vec <- c(vec, "-D")
    if (n) 
        vec <- c(vec, "-n")
    if (C) 
        vec <- c(vec, "-C")
    vec <- c(vec, "-q", as.character(q))
    vec <- c(vec, "-f", as.character(f))
    vec <- c(vec, "-k", as.character(k))
    vec <- c(vec, "-c", as.character(c))
    vec <- c(vec, "-o", as.character(o))
    
    unloadNamespace("IRISFGM")
    ret <- .main(vec)
    if (ret == 42) 
        return(IRISFGM::qubic(paste0(i, ".chars"), d = TRUE))
    # return (ret)
}


.onUnload <- function(libpath) {
    library.dynam.unload("IRISFGM", libpath)
}






