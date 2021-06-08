ZbSpline <- function(x, knots, k = 1) {
  ZbSplineSingle <- function(x, knots, k = 1) {
    k0 <- knots[k]
    k1 <- knots[k + 1]
    if ((x > k0) && (x <= k1)) {
      return(1)
    }
    return(0)
  }
  return(sapply(x, function(z)
    ZbSplineSingle(z, knots, k)))
}

LbSpline <- function(x, knots, k = 1) {
  LbSplineSingle <- function(x, knots, k = 1) {
    k0 <- knots[k]
    k1 <- knots[k + 1]
    k2 <- knots[k + 2]
    f1 <- function(x)
      (x - k0) / (k1 - k0)
    f2 <- function(x)
      (k2 - x) / (k2 - k1)
    if ((x > k0) && (x <= k1)) {
      return(f1(x))
    }
    if ((x > k1) && (x <= k2)) {
      return(f2(x))
    }
    return(0)
  }
  return(sapply(x, function(z)
    LbSplineSingle(z, knots, k)))
}

QbSpline <- function(x, knots, k = 1) {
  QbSplineSingle <- function(x, knots, k = 1) {
    k0 <- knots[k]
    k1 <- knots[k + 1]
    k2 <- knots[k + 2]
    k3 <- knots[k + 3]
    f1 <- function(x)
      ((x - k0) ^ 2) / ((k2 - k0) * (k1 - k0))
    f2 <- function(x) {
      term1 <- ((x - k0) * (k2 - x)) / ((k2 - k0) * (k2 - k1))
      term2 <- ((x - k1) * (k3 - x)) / ((k3 - k1) * (k2 - k1))
      return(term1 + term2)
    }
    f3 <- function(x)
      ((k3 - x) ^ 2) / ((k3 - k1) * (k3 - k2))
    if ((x > k0) && (x <= k1)) {
      return(f1(x))
    }
    if ((x > k1) && (x <= k2)) {
      return(f2(x))
    }
    if ((x > k2) && (x <= k3)) {
      return(f3(x))
    }
    return(0)
  }
  return(sapply(x, function(z)
    QbSplineSingle(z, knots, k)))
}

IZbSpline <- function(x, knots, k = 1) {
  IZbSplineSingle <- function(x, knots, k = 1) {
    k0 <- knots[k]
    k1 <- knots[k + 1]
    if (x <= k0) {
      return (0)
    }
    if ((x > k0) && (x <= k1)) {
      return((x - k0) / (k1 - k0))
    }
    if (x > k1) {
      return (1)
    }
  }
  return(sapply(x, function(z)
    IZbSplineSingle(z, knots, k)))
}

ILbSpline <- function(x, knots, k = 1) {
  ILbSplineSingle <- function(x, knots, k = 1) {
    k0 <- knots[k]
    k1 <- knots[k + 1]
    k2 <- knots[k + 2]
    f1 <- function(x)
      ((x - k0) ^ 2) / ((k1 - k0) * (k2 - k0))
    f2 <-
      function(x)
        ((x - k0) / (k2 - k0)) + ((x - k1) * (k2 - x)) / ((k2 - k0) * (k2 - k1))
    if (x <= k0) {
      return (0)
    }
    if ((x > k0) && (x <= k1)) {
      return(f1(x))
    }
    if ((x > k1) && (x <= k2)) {
      return(f2(x))
    }
    if (x > k2) {
      return (1)
    }
    return(0)
  }
  return(sapply(x, function(z)
    ILbSplineSingle(z, knots, k)))
}
