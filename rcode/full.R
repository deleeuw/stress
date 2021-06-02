library(MASS)

torgerson <- function(delta, p = 2) {
  doubleCenter <- function(x) {
    n <- dim(x)[1]
    m <- dim(x)[2]
    s <- sum(x) / (n * m)
    xr <- rowSums(x) / m
    xc <- colSums(x) / n
    return((x - outer(xr, xc, "+")) + s)
  }
  z <-
    eigen(-doubleCenter((as.matrix (delta) ^ 2) / 2), symmetric = TRUE)
  v <- pmax(z$values, 0)
  return(z$vectors[, 1:p] %*% diag(sqrt(v[1:p])))
}


makeA <- function (n) {
  m <- n * (n - 1) / 2
  a <- list()
  for (j in 1:(n - 1))
    for (i in (j + 1):n) {
      d <- ein (i, n) -ein (j, n)
      e <- outer (d, d)
      a <- c(a, list (e))
    }
  return (a)
}

makeD <- function (a, x) {
  return (sapply (a, function (z)
    sqrt (sum (x * (
      z %*% x
    )))))
}

makeB <- function (w, delta, d, a) {
  n <- length (a)
  m <- nrow (a[[1]])
  b <- matrix (0, m , m)
  for (i in 1:n)
    b <- b + w[i] * (delta[i] / d[i]) * a[[i]]
  return (b)
}

makeV <- function (w, a) {
  n <- length (a)
  m <- nrow (a[[1]])
  v <- matrix (0, m, m)
  for (i in 1:n)
    v <- v + w[i] * a[[i]]
  return (v)
}

inBetween <- function (alpha, beta, x, y, w, delta, a) {
  z <- alpha * x + beta * y
  d <- makeD (a, z)
  return (sum (w * (delta - d) ^ 2))
}

biBase <- function (x, y, a) {
  biBi <- function (x, y, v) {
    a11 <- sum (x * (v %*% x))
    a12 <- sum (x * (v %*% y))
    a22 <- sum (y * (v %*% y))
    return (matrix (c(a11, a12, a12, a22), 2, 2))
  }
  return (lapply (a, function (u)
    biBi (x, y, u)))
}

fullMDS <-
  function (delta,
            w = rep (1, length (delta)),
            xini,
            a,
            itmax = 100,
            eps = 1e-6,
            verbose = TRUE) {
    m <- length (a)
    v <- makeV (w, a)
    vv <- ginv (v)
    xold <- xini
    dold <- makeD (a, xini)
    sold <- sum ((delta - dold) ^ 2)
    bold <- makeB (w, delta, dold, a)
    itel <- 1
    repeat {
      xnew <- vv %*% bold %*% xold
      dnew <- makeD (a, xnew)
      bnew <- makeB (w, delta, dnew, a)
      snew <- sum ((delta - dnew) ^ 2)
      if (verbose) {
        cat (
          formatC (itel, width = 4, format = "d"),
          formatC (
            sold,
            digits = 10,
            width = 13,
            format = "f"
          ),
          formatC (
            snew,
            digits = 10,
            width = 13,
            format = "f"
          ),
          "\n"
        )
      }
      if ((itel == itmax) || (abs(sold - snew) < eps))
        break
      itel <- itel + 1
      xold <- xnew
      dold <- dnew
      sold <- snew
      bold <- bnew
    }
    return (list (
      x = xnew,
      d = dnew,
      delta = delta,
      s = snew,
      b = bnew,
      v = v,
      itel = itel
    ))
  }
