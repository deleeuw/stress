dyn.load("deboorU.so")

checkIncreasing <- function (innerknots, lowend, highend) {
  h <- .C(
    "checkIncreasing",
    as.double(innerknots),
    as.double (lowend),
    as.double (highend),
    as.integer(length (innerknots)),
    fail = as.integer(0)
  )
  return (h$fail)
}

extendPartition <-
  function (innerknots,
            multiplicities,
            order, 
            lowend,
            highend) {
    ninner <- length (innerknots)
    kk <- sum(multiplicities)
    nextended <- kk + 2 * order
    if (max (multiplicities) > order)
      stop ("multiplicities too high")
    if (min (multiplicities) < 1)
      stop ("multiplicities too low")
    if (checkIncreasing (innerknots, lowend, highend))
      stop ("knot sequence not increasing")
    h <-
      .C(
        "extendPartition",
        as.double (innerknots),
        as.integer (multiplicities),
        as.integer (order),
        as.integer (ninner),
        as.double (lowend),
        as.double (highend),
        extendPartition = as.double (rep (0, nextended))
      )
    return (h)
  }

bisect <- function (x, knots, lowindex = 1, highindex = length (knots)) {
  h <- .C("bisect",
     as.double (x),
     as.double (knots),
     as.integer (lowindex),
     as.integer (highindex),
     index = as.integer (0))
  return (h$index)
}

bsplines <- function (x, knots, order) {
  if ((x > knots[length(knots)]) || (x < knots[1]))
    stop ("argument out of range")
  h <- .C("bsplines",
          as.double (x),
          as.double (knots),
          as.integer (order),
          as.integer (length (knots)),
          index = as.integer (0),
          q = as.double (rep(0, order)))
  return (list (q = h$q, index = h$ind))
}

bsplineBasis <- function (x, innerknots, multiplicities = rep (1, length (innerknots)), order, lowend, highend) {
  n <- length (x)
  k <- sum (multiplicities)
  m <- k + order + 1
  result <- rep (0, n * m)
  h <- .C("bsplineBasis",
          as.double (x),
          as.double (innerknots),
          as.integer (multiplicities),
          as.integer (order),
          as.integer (length(innerknots)),
          as.integer (n),
          as.double (lowend),
          as.double (highend),
          result = as.double (result))
  return (matrix (h$result, n, m))
}

bsplineBasisU <- function (x, knots, order) {
  n <- length (x)
  k <- length (knots)
  m <- k - order
  result <- rep (0, n * m)
  h <- .C("bsplineBasis",
          as.double (x),
          as.double (knots),
          as.integer (order),
          as.integer (k),
          as.integer (n),
          result = as.double (result))
  return (matrix (h$result, n, m))
}
