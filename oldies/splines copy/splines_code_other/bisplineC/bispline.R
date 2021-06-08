dyn.load ("bispline.so")

bispline <-
  function (x,
            degree,
            innerknots,
            lowknot = min(x),
            highknot = max(x),
            norm = 1) {
    innerknots <- unique (sort (innerknots))
    knots <-
      c(rep(lowknot, degree + 1), innerknots, rep(highknot, degree + 1))
    n <- length (knots)
    v <- rep (0, n)
    vint <- rep (0, n)
    h <-
      .C(
        "bisplineC",
        as.integer (n),
        as.double (knots),
        as.integer (degree + 1),
        as.integer (norm),
        as.double (x),
        jint = as.integer (0),
        v = as.double (v),
        vint = as.double (vint),
        ifail = as.integer (0)
      )
    if (h$ifail == 1)
      stop("degree too high for number of knots")
    if (h$ifail == 2)
      stop("degree negative")
    if (h$ifail == 3)
      stop("knots not sorted")
    if (h$ifail == 4)
      stop("too many knots coincide")
    if (h$ifail == 5)
      stop("norm not one or two")
    return (list (
      jint = h$jint,
      v = h$v,
      vint = h$vint
    ))
  }

bisplineBasis <-
  function (x,
            degree,
            innerknots,
            lowknot = min(x, innerknots),
            highknot = max(x, innerknots),
            norm = 2) {
    innerknots <- unique (sort (innerknots))
    knots <-
      c(rep(lowknot, degree + 1), innerknots, rep(highknot, degree + 1))
    n <- length (knots)
    m <- length (x)
    bbasis <- matrix (0, m, n)
    ibasis <- matrix (0, m, n)
    h <- .C (
      "bisplineBasisC",
      as.integer (n),
      as.double (knots),
      as.integer (degree + 1),
      as.integer (m),
      as.double (x),
      as.integer (norm),
      bbasis = as.double (bbasis),
      ibasis = as.double (ibasis)
    )
    bbasis <- matrix (h$bbasis, m, n)
    bbasis <- bbasis[, which (colSums (bbasis) > 0), drop = FALSE]
    ibasis <- matrix (h$ibasis, m, n)
    ibasis <- ibasis[, which (colSums (ibasis) > 0), drop = FALSE]
    ibasis <- iComplete (x, ibasis, norm)
    return (list (bbasis = bbasis, ibasis = ibasis))
  }

iComplete <- function (x, h, norm) {
  m <- ncol (h)
  for (j in 1:m) {
    wj <- which.max (h[, j])
    mj <- max (h[, j])
    h[which (x >= x[wj]), j] <- mj
  }
  return (h)
}
