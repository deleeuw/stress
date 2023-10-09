tau <- function (x) {
  return (- 0.5 * dcen (x))
}

kappa <- function (x) {
  return (outer (diag (x), diag (x), "+") - 2 * x)
}

fcmds <-
  function (delta,
            xold,
            ninner = 1,
            itmax = 100,
            eps = 1e-6,
            verbose = TRUE) {
    itel <- 0
    p <- ncol (xold)
    xold <- apply (xold, 2, function (x)
      x - mean(x))
    xold <- qr.Q (qr (xold))
    repeat {
      xinn <- xold
      for (i in 1:ninner) {
        xnew <- delta %*% xinn
        xnew <- -apply (xnew, 2, function (x)
          x - mean (x)) / 2
        xinn <- xnew
        itel <- itel + 1
      }
      qnew <- qr (xnew)
      xnew <- qr.Q (qnew)
      rnew <- qr.R (qnew)
      epsi <- 2 * p - 2 * sum (svd (crossprod (xold, xnew))$d)
      if (verbose) {
        cat(
          "itel ",
          formatC (
            itel,
            digits = 4,
            width = 6,
            format = "d"
          ),
          "epsi ",
          formatC (
            epsi,
            digits = 10,
            width = 15,
            format = "f"
          ),
          "\n"
        )
      }
      if ((epsi < eps) || (itel == itmax))
        break
      xold <- xnew
    }
    return (list (x = xnew, r = rnew, itel = itel))
  }

treq <- function (x) {
  n <- nrow (d)
  m <- -Inf
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      for (k in 1:n) {
        if ((k == i) || (k == j))
          next
        m <- max (m, x[i, j] - (x[i, k] + x[k, j]))
      }
    }
  }
  return (m)
}

acbound <- function (d) {
  n <- nrow (d)
  s <- qr.Q (qr (cbind (1, matrix (rnorm (
    n * (n - 1)
  ), n, n - 1))))
  k <- tau (d * d)
  l <- 2 * tau (d)
  m <- jmat (n) / 2
  ma <- -Inf
  for (i in 2:n)
    for (j in 1:(i - 1)) {
      v <- solve(polynomial(c(k[i, j], l[i, j], m[i, j])))
      ma <- max(ma, max(v))
    }
  return (list(ma = ma, mw = k + ma * l + m * ma ^ 2))
}

aceval <- function (d, bnd = c(-10, 10)) {
  n <- nrow (d)
  k <- tau (d * d)
  l <- 2 * tau (d)
  m <- jmat (n) / 2
  s <- qr.Q (qr (cbind (1, matrix (rnorm (
    n * (n - 1)
  ), n, n - 1))))
  kc <- (crossprod (s, k) %*% s)[-1,-1]
  lc <- (crossprod (s, l) %*% s)[-1,-1]
  mc <- (crossprod (s, m) %*% s)[-1,-1]
  a <- seq(bnd[1], bnd[2], length = 1000)
  b <- rep(0, 1000)
  for (i in 1:1000) {
    ww <- kc + lc * a[i] + mc * (a[i] ^ 2)
    b[i] <- min (eigen(ww)$values)
  }
  return (list(a = a, b = b))
}

acqep <- function(d) {
  n <- nrow (d)
  k <- tau (d * d)
  l <- 2 * tau (d)
  m <- jmat (n) / 2
  s <- qr.Q (qr (cbind (1, matrix (rnorm (
    n * (n - 1)
  ), n, n - 1))))
  nn <- n - 1
  ns <- 1:nn
  kc <- (crossprod (s, k) %*% s)[-1,-1]
  lc <- (crossprod (s, l) %*% s)[-1,-1]
  mc <- (crossprod (s, m) %*% s)[-1,-1]
  ma <- matrix(0, 2 * nn, 2 * nn)
  ma[ns, nn + ns] <- diag (n - 1)
  ma[nn + ns, ns] <- -2 * kc
  ma[nn + ns, nn + ns] <- -2 * lc
  return (list (ma = ma, me = eigen(ma)$values))
}
