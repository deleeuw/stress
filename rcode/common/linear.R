## invert with pivot
## invert with bordering
## solve wth bordering
## solve with pivoting
## jacobi
## QR

gramy <- function (y, v) {
  r <- length (y)
  s <- sum (y[[1]] * (v %*% y[[1]]))
  y[[1]] <- y[[1]] / sqrt (s)
  for (j in 2:r) {
    for (i in 1:(j - 1)) {
      s <- sum (y[[i]] * (v %*% y[[j]]))
      y[[j]] <- y[[j]] - s * y[[i]]
    }
    s <- sum (y[[j]] * v %*% y[[j]])
    y[[j]] <- y[[j]] / sqrt (s)
  }
  return (y)
}



hinv <- function(x) {
  return (apply (x, c(1, 2), function (a)
    ifelse (a == 0, 0, 1 / a)))
}


circular <- function (n) {
  x <- seq (0, 2 * pi, length = n + 1)
  z <- matrix (0, n + 1, 2)
  z[, 1] <- sin (x)
  z[, 2] <- cos (x)
  return (z[-1, ])
}

direct_sum <- function (x) {
  n <- length (x)
  nr <- sapply (x, nrow)
  nc <- sapply (x, ncol)
  s <- matrix (0, sum (nr), sum (nc))
  k <- 0
  l <- 0
  for (j in 1:n) {
    s[k + (1:nr[j]), l + (1:nc [j])] <- x[[j]]
    k <- k + nr[j]
    l <- l + nc[j]
  }
  return (s)
}


