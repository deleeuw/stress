bmat <- function (a, b, x, y, delta) {
  bm <- matrix (0, 2, 2)
  hm <- matrix (0, 2, 2)
  z <- c(a, b)
  for (i in 1:4) {
    for (j in 1:4) {
      if (i == j) next
      uij <- uu (i, j, x, y)
      uz <- drop (uij %*% z)
      dij <- sqrt (sum (uij * outer (z, z)))
      bm <- bm + (delta[i,j] / dij) * uij
      hm <- hm + (delta[i,j] / dij) * (uij - outer (uz, uz) / sum (z * uz))
    }
  }
  return (list (b = bm, h = hm))
}

stress <- function (a, b, x, y, delta) {
  z <- c (a, b)
  bm <- bmat (a, b, x, y, delta)$b
  return (1 + sum(z ^ 2) / 2 - sum (z * bm %*% z))
}

rho <- function (a, b, x, y, delta) {
  z <- c (a, b)
  bm <- bmat (a, b, x, y, delta)$b
  return (sum (z * bm %*% z))
}

vv <- function (i, j, x, y) {
  a <- matrix (0, 2, 2)
  a[1, 1] <- sum ((x[i, ]- x[j,]) ^ 2)
  a[2, 2] <- sum ((y[i, ]- y[j,]) ^ 2)
  a[1, 2] <- a[2, 1] <- sum ((x[i, ]- x[j,]) * (y[i, ]- y[j, ]))
  return (a)
}

uu <- function (i, j, x, y) {
  n <- nrow (x)
  asum <- 2 * n * matrix (c (sum(x ^ 2), sum (x * y), sum (x * y), sum (y ^ 2)), 2, 2)
  csum <- solve (chol (asum))
  return (t(csum) %*% vv (i, j, x, y) %*% csum)
}

smacof <- function (a, b, x, y, delta, eps = 1e-10, itmax = 1000, verbose = TRUE) {
  zold <- c(a,b)
  bold <- bmat (a, b, x, y, delta)$b
  fold <- 1 + sum(zold ^ 2) / 2 - sum (zold * bold %*% zold)
  itel <- 1
  repeat {
    znew <- drop (bold %*% zold)
    bhmt <- bmat (znew[1], znew[2], x, y, delta)
    bnew <- bhmt$b
    fnew <- 1 + sum(znew ^ 2) / 2 - sum (znew * bnew %*% znew)
    if (verbose) {
      cat (
        formatC (itel, width = 4, format = "d"),
        formatC (
          fold,
          digits = 10,
          width = 13,
          format = "f"
        ),
        formatC (
          fnew,
          digits = 10,
          width = 13,
          format = "f"
        ),
        "\n"
      )
    }
    if ((itel == itmax) || (fold - fnew) < eps)
      break ()
    itel <- itel + 1
    fold <- fnew
    zold <- znew
    bold <- bnew
  }
  return (list (stress = fnew, theta = znew, itel= itel, b = bnew, g = znew - bnew %*% znew, h = diag(2) - bhmt$h))
}


newton <- function (a, b, x, y, delta, eps = 1e-10, itmax = 1000, verbose = TRUE) {
  zold <- c(a,b)
  bhmt <- bmat (a, b, x, y, delta)
  bold <- bhmt$b
  hold <- diag(2) - bhmt$h
  fold <- 1 + sum(zold ^ 2) / 2 - sum (zold * bold %*% zold)
  itel <- 1
  repeat {
    znew <- drop (solve (hold, bold %*% zold))
    bhmt <- bmat (znew[1], znew[2], x, y, delta)
    bnew <- bhmt$b
    hnew <- diag(nrow(bnew)) - bhmt$h
    fnew <- 1 + sum(znew ^ 2) / 2 - sum (znew * bnew %*% znew)
    if (verbose) {
      cat (
        formatC (itel, width = 4, format = "d"),
        formatC (
          fold,
          digits = 10,
          width = 13,
          format = "f"
        ),
        formatC (
          fnew,
          digits = 10,
          width = 13,
          format = "f"
        ),
        "\n"
      )
    }
    if ((itel == itmax) || abs (fold - fnew) < eps)
      break ()
    itel <- itel + 1
    fold <- fnew
    zold <- znew
    bold <- bnew
    hold <- hnew
  }
  return (list (stress = fnew, theta = znew, itel = itel, b = bnew, g = znew - bnew %*% znew, h = hnew))
}

mprint <- function (x, d = 2, w = 5) {
  print (noquote (formatC (x, di = d, wi = w, fo = "f")))
}
