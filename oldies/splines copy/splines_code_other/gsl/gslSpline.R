dyn.load("gslSpline.so")

gslSpline <- function (x, k, br) {
  nbr <- length (br)
  nrs <- k + nbr - 2
  res <-
    .C(
      "BSPLINE",
      as.double (x),
      as.integer (k),
      as.integer (nbr),
      as.double (br),
      results = as.double (rep(0.0, nrs))
    )
  return (res$results)
}
