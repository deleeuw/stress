dyn.load ("amalgm.so")

amalgm <- function (x, w = rep (1, length (x))) {
    n <- length (x)
    a <- rep (0, n)
    b <- rep (0, n)
    y <- rep (0, n)
    lf <- .Fortran ("AMALGM", n = as.integer (n), x = as.double (x), w = as.double (w), a = as.double (a), b = as.double (b), y = as.double (y), tol = as.double(1e-15), ifault = as.integer(0))
    return (lf$y)
}
