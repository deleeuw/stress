
jacobi <- function (a, itmax = 100, eps = 1e-10) {
  m <- length (a)
  n <- (sqrt (1 + 8 * m) - 1) / 2
  i <- 1:n
  j <- (-(i ^ 2) / 2) + (n + (3 / 2)) * i - n
  h <-
    .C(
      "jacobiC",
      as.integer (n),
      a = as.double (a),
      e = as.double (rep (0, n * n)),
      as.double (rep(0, n)),
      as.double (rep(0, n)),
      as.integer (itmax),
      as.double (eps)
    )
  return (list (values = h$a[j], vectors = matrix (h$e, n, n)))
}
