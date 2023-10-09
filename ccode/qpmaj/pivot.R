dyn.load ("pivot.so")

pivotOneRC <- function (a, k, type = 3, eps = 1e-10) {
  n <- nrow (a)
  m <- ncol (a)
  h <-
    .C(
      "pivotOneC",
      a = as.double (a),
      k = as.integer(k),
      n = as.integer (n),
      m = as.integer (m),
      type = as.integer (type),
      refuse = as.integer (0),
      eps = as.double (eps)
    )
  return (list(a = matrix(h$a, n, m), refuse = h$refuse))
}

pivotRC <- function (a, ind, type = 3, eps = 1e-10) {
  n <- nrow (a)
  m <- ncol (a)
  p <- length (ind)
  done <- rep (0, p)
  h <-
    .C(
      "pivotC",
      a = as.double (a),
      ind = as.integer (ind),
      jnd = as.integer (rep (0, p)),
      done = as.integer (rep (0, p)),
      n = as.integer (n),
      m = as.integer (m),
      p = as.integer (p),
      type = as.integer (type),
      skip = as.integer (rep (0, p)),
      eps = as.double (eps)
    )
  return (list(a = matrix(h$a, n, m), pivot = h$jnd, skip = h$skip))
}
