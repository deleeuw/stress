

dposvR <- function (a, b) {
  n <- nrow (a)
  m <- max (ncol (b), 1)
  h <-
    .C("dposv",
       as.integer(n),
       as.integer(m),
       as.double (a),
       as.double (b))
  if (is.null(ncol(b))) {
    return (h[[4]])
  } else {
    return (matrix(h[[4]], n, m))
  }
}

dsyevdR <- function (a) {
  n <- nrow (a)
  x <- matrix (0, n, n)
  h <- .C("dsyevd", as.integer(n), as.double (a), as.double (x))
  return (list(values = h[[3]][1:n], vectors = matrix(h[[2]], n, n)))
}

dgeqrR <- function (x) {
  n <- nrow (x)
  m <- ncol (x)
  h <-
    .C("dgeqrf",
       as.integer(n),
       as.integer(m),
       z = as.double (x),
       tau = as.double (rep (0, m)))
  return (list (z = matrix(h$z, n, m), tau = h$tau))
}

dorgqrR <- function (z, tau) {
  n <- nrow (z)
  m <- ncol (z)
  h <-
    .C("dorgqr",
       as.integer(n),
       as.integer(m),
       as.double (z),
       as.double (tau))
  return (matrix(h[[3]], n, m))
}

dorthoR <- function(x) {
  n <- nrow (x)
  m <- ncol (x)
  h <-
    .C("dortho",
       as.integer(n),
       as.integer(m),
       as.double (x),
       as.double(rep(0, m)))
  return (matrix(h[[3]], n, m))
}
