dyn.load ("splinebasis.so")

splineBasis <- function (x, degree, innerknots, lowknot = min(x), highknot = max(x), type = 1) {
    knots <- unique (sort (innerknots))
    knots <- c(rep(lowknot, degree + 1), knots, rep (highknot, degree + 1))
    n <- length (x)
    m <- length (knots)
    nf <- m - (degree + 1)
    basis <- rep (0,  n * nf)
    res <- .C("splinebasis", d = as.integer(degree),
        n = as.integer(n), m = as.integer (m), x = as.double (x), knots = as.double (knots), as.integer (type), basis = as.double (basis))
    basis <- matrix (res$basis, n, nf)
#    basis <- basis[, which(colSums(basis) > 0), drop = FALSE]
    return (basis)
}

knotsQ <- function (x, n = 5) {
  do <- function (z, n) {
    y <- quantile (z, probs = seq(0, 1, length = max (2, n)))
    return (y[-c(1, length(y))])
  }
  if (ncol(x) > 0)
    n <- rep (n, ncol (x))
  if (ncol (x) > 0)
    lapply (1:ncol(x), function (i) do (x[,i], n[i]))
}

knotsR <- function (x, n = rep (5, ncol (x))) {
  do <- function (i) {
    y <- seq (min(x[, i]), max(x[, i]), length = max (2, n[i]))
    return (y[-c(1, length(y))])
  }
  lapply (1:ncol(x),  function (i)
    do (i))
}

knotsE <- function (x = NULL) {
  lapply (1:max(1, ncol(x)), function (i)
    numeric(0))
}

knotsD <- function (x) {
  do <- function (i) {
    y <- sort (unique (x[, i]))
    return (y[-c(1, length(y))])
  }
  lapply (1:ncol(x),  function (i)
    do (i))
}

stepPlotter <- function (x, y, knots, xlab) {
  y <- as.matrix (y)
  plot (x,
        y[, 1],
        type = "n",
        xlab = xlab,
        ylab = "Transform")
  nknots <- length (knots)
  knots <- c(min(x) - 1, knots, max(x) + 1)
  for (i in 1:(nknots + 1)) {
    ind <- which ((x >= knots [i]) & (x < knots[i + 1]))
    lev <- median (y [ind, 1])
    lines (rbind (c(knots[i], lev), c (knots[i + 1], lev)), col = "RED", lwd = 3)
    if (ncol (y) == 2) {
      lev <- median (y [ind, 2])
      lines (rbind (c(knots[i], lev), c (knots[i + 1], lev)), col = "BLUE", lwd = 3)
    }
  }
}

