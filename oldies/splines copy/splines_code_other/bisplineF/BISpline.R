

dyn.load ("BISpline.so")

BISpline <-
  function (x,
            degree,
            innerknots,
            lowknot = min(x),
            highknot = max(x),
            type = 1) {
    knots <- unique (sort (innerknots))
    order <- degree + 1
    knots <- c(rep(lowknot, order), knots, rep (highknot, order))
    n <- length (x)
    m <- length (knots)
    nf <- m -  order
    basis <- matrix (0, n, nf)
    iasis <- matrix (0, n, nf)
    h <-
      .Fortran(
        "BSPLINE",
        as.integer (n),
        as.double (x) ,
        as.integer (m),
        as.double (knots),
        as.integer (order),
        as.integer (type),
        basis = as.double (basis),
        iasis = as.double (iasis)
      )
    basis <- matrix (h$basis, n, nf)
    iasis <- matrix (h$iasis, n, nf)
    basis <- basis[, which (colSums (basis) > 0), drop = FALSE]
    iasis <- iasis[, which (colSums (iasis) > 0), drop = FALSE]
    return (list (basis = basis, iasis = order * iasis))
  }
