dyn.load("smatrix.so")

symmetricCholesky <- function (a, eps = 1e-10) {
  m <- length (a)
  n <- (sqrt (1 + 8 * m) - 1) / 2
  h <- .C(
    "scholesky",
    n = as.integer(n),
    a = as.double (a),
    singular = as.integer(0),
    indefinite = as.integer(0),
    det = as.double (0),
    eps = as.double(eps)
  )
  return (h)
}

pivotCholesky <- function (a, eps = 1e-10) {
  m <- length (a)
  n <- (sqrt (1 + 8 * m) - 1) / 2
  h <- .C(
    "pcholesky",
    n = as.integer(n),
    a = as.double (a),
    det = as.double (0),
    rank = as.integer(0),
    order = as.integer(1:n),
    indefinite = as.integer(0),
    eps = as.double(eps)
  )
  return (h)
}

symmetricSweep <- function (a, eps = 1e-10) {
  m <- length (a)
  n <- (sqrt (1 + 8 * m) - 1) / 2
  h <- .C(
    "ssweep",
    n = as.integer(n),
    a = as.double (a),
    det = as.double (0),
    rank = as.integer(0),
    singular = as.integer(0),
    eps = as.double (eps)
  )
  return (h)
}

pivotSweep <- function (a, eps = 1e-10) {
  m <- length (a)
  n <- (sqrt (1 + 8 * m) - 1) / 2
  h <- .C(
    "psweep",
    n = as.integer(n),
    a = as.double (a),
    det = as.double (0),
    rank = as.integer(0),
    order = as.integer(1:n),
    singular = as.integer(0),
    eps = as.double (eps)
  )
  return (h)
}


triangle2matrix <- function (x) {
  m <- length (x)
  n <- round ((sqrt (1 + 8 * m) - 1) / 2)
  h <-
    .C("trimat", as.integer (n), as.double (x), as.double (rep (0, n * n)))
  return (matrix(h[[3]], n, n))
}

matrix2triangle <- function (x) {
  n <- dim(x)[1]
  m <- n * (n + 1) / 2
  h <-
    .C("mattri", as.integer (n), as.double (x), as.double (rep (0, m)))
  return (h[[3]])
}

triangle2triangle <- function (x) {
  m <- length (x)
  n <- round ((sqrt (1 + 8 * m) - 1) / 2)
  h <-
    .C("tritri", as.integer (n), as.double (x), as.double (rep (0, n * n)))
  return (matrix(h[[3]], n, n))
}

matrix2print <- function (x, w = 15, p = 10) {
  n <- nrow (x)
  m <- ncol (x)
  h <-
    .C(
      "primat",
      as.integer(n),
      as.integer(m),
      as.integer(w),
      as.integer(p),
      as.double (x)
    )
}
