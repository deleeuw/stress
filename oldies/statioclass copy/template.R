
stress <- function (x, delta, w = NULL) {
  x <- as.matrix(x)
  n <- nrow (x)
  if (is.null (w)) {
    w <- as.dist(matrix(1, n, n))
  }
  d <- dist (x)
  r <- sum (w * d * delta)
  e <- sum (w * d ^ 2)
  return (1 - r + e / 2)
}


