f <- function (i, j, n) {
  if (j >= i) {
    stop ("no such element")
  }
  return (i + (j - 1) * n - choose (j + 1, 2))
}

g <- function (k, n) {
  if (k > 0.5 * n * (n - 1)) {
    stop("no such element")
  }
  j <- m <- 1
  while (k >= ((j * n) - m + 1)) {
    j <- j + 1
    m <- m + j
  }    
  i <- k - (j - 1) * n + m
  return (c(i, j))
}

testf <- function (n) {
  for (j in 1:(n-1)) {
    for (i in (j + 1):n) {
      print(f(i, j, n))
    }
  }
}

testg <- function (n) {
  for (i in 1:(n*(n-1)/2)) {
    print (g(i, n))
  }
}
