dyn.load ("dbspvt.so")

bsplvt <- function (x, order, knots, index = 1) {
  left <- which.min (x >= knots)
  biatx <- dbiatx <- rep (0, order)
  h <-
    .Fortran (
      "DBSPVT",
      as.double (knots),
      as.integer (order),
      as.integer (index),
      as.double(x),
      as.integer (left),
      as.double (biatx),
      as.double (dbiatx)
    )
  return (h)
}
