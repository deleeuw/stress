dyn.load("sinhaSpline.so")

sinhaBasis <-
  function (x, degree, innerknots, lowknot = min(x,innerknots), highknot = max(x,innerknots)) {
    innerknots <- unique (sort (innerknots))
    knots <-
      c(rep(lowknot, degree + 1), innerknots, rep(highknot, degree + 1))
    n <- length (x)
    m <- length (innerknots) + 2 * (degree + 1)
    nf <- length (innerknots) + degree + 1
    basis <- rep (0,  n * nf)
    res <- .C(
      "sinhaBasis", d = as.integer(degree),
      n = as.integer(n), m = as.integer (m), x = as.double (x), knots = as.double (knots), basis = as.double(basis)
    )
    basis <- matrix (res$basis, n, nf)
    basis <- basis[,which(colSums(basis) > 0)]
    return (basis)
  }

