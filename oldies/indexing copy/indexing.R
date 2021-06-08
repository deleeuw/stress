
fArrayFirstR <- function (cell, dimension) {
  rank <- length (cell)
  k <- cumprod (c(1, dimension))[1:rank]
  return (1 + sum ((cell - 1) * k))
}

fSupSymIncreasingFirstR <- function (cell) {
  f <- function (r)
    choose (r + (cell[r] - 1) - 1, r)
  rank <- length (cell)
  cell <- sort (cell)
  return (1 + sum (sapply (1:rank, f)))
}

fArrayFirstInverseR <- function(index, dimension) {
  rank <- length (dimension)
  b <- cumprod (c(1, dimension))[1:rank]
  r <- rep(0, length(b))
  s <- index - 1
  for (j in rank:1) {
    r[j] <- s %/% b[j]
    s <- s - r[j] * b[j]
  }
  return(1 + r)
}

fSupSymIncreasingFirstInverseR <- function (dimension, rank, index) {
  last.true <- function (x) {
    w <- which (x)
    return (w[length(w)])
  }
  a <- rep (0, rank)
  v <- index - 1
  for (k in rank:1) {
    s <- choose (k + (0:(dimension - 1)) - 1, k)
    u <- last.true (v >= s)
    a[k] <- u
    v <- v - s[u]
  }
  return (a)
}
