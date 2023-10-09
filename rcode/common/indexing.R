
kron <- function (i, j) {
  return (ifelse (i == j, 1, 0))
}

ein <- function (i, n) {
  return (ifelse (i == 1:n, 1, 0))
}

aijn <- function (i, j, n) {
  dif <- ein (i, n) - ein (j, n)
  return (outer (dif, dif))
}

jmat <- function (n) {
  return (diag(n) - 1 / n)
}

ccen <- function (x) {
  return (apply (x, 2, function (y)
    y - mean (y)))
}

repList <- function(x, n) {
  z <- list()
  for (i in 1:n)
    z <- c(z, list(x))
  return(z)
}

rcen <- function (x) {
  return (t (apply (x, 1, function (y)
    y - mean (y))))
}

dcen <- function (x) {
  return (ccen (rcen (x)))
}

wdef <- function (n) {
  return (1 - diag (n))
}


lower_triangle <- function (x) {
  n <- nrow (x)
  return (x[outer(1:n, 1:n, ">")])
}

fill_symmetric <- function (x) {
  m <- length (x)
  n <- (1 + sqrt (1 + 8 * m)) / 2
  d <- matrix (0, n, n)
  d[outer(1:n, 1:n, ">")] <- x
  return (d + t(d))
}

