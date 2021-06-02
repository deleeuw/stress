dummy <- function () {
  set.seed(12345)
  
  x <- matrix (rnorm(16), 8, 2)
  x <- apply (x, 2, function (x)
    x - mean (x))
  y <- matrix (rnorm(10), 5, 2)
  a <- rowSums (x ^ 2)
  b <- rowSums (y ^ 2)
  d <- sqrt (outer(a, b, "+") - 2 * tcrossprod (x, y))
  
  set.seed (12345)
  x <- matrix(rnorm(10), 5, 2)
  x <- apply (x, 2, function (x)
    x - mean (x))
  x <- qr.Q(qr(x))
  y <- matrix(rnorm(12), 6, 2)
  v <- apply (y, 2, mean)
  print (v)
  dx <- diag(tcrossprod(x))
  dy <- diag(tcrossprod(y))
  xy <- tcrossprod(x, y)
  dd <- outer(dx, dy, "+") - 2 * xy
  j5 <- diag(5) - 1 / 5
  j6 <- diag(6) - 1 / 6
  dc <- -(j5 %*% dd %*% j6) / 2
  sv <- svd (dc, nu = 2, nv = 2)
  xs <- sv$u
  ys <- sv$v %*% diag (sv$d[1:2])
  tt <- crossprod (x, xs)
  dk <- diag (tcrossprod(xs))
  dl <- diag (tcrossprod(ys))
  dr <- dd - outer (dk, dl, "+") + 2 * tcrossprod (xs, ys)
  print (dr)
}

schoenemann <- function (delta, p) {
  n <- nrow (delta)
  m <- ncol (delta)
  l <- p * (p + 1) / 2
  d <- delta ^ 2
  e <- torgerson (d)
  q <- svd (e, nu = p, nv = p)
  g <- q$u
  h <- q$v %*% diag (q$d[1:p])
  f <- d + 2 * tcrossprod(g, h)
  a <- apply (ccen (f), 1, mean)
  r <- matrix (0, n, l)
  k <- 1
  for (i in 1:p) {
    for (j in 1:i) {
      if (i == j) {
        r[, k] <- g[, i] ^ 2
      } else {
        r[, k] <- 2 * g[, i] * g[, j]
      }
      k <- k + 1
    }
  }
  lhs <- cbind (ccen(r), ccen (-2 * g))
  b <- lm.fit (lhs, a)$coefficients
  k <- 1
  s <- matrix (0, p, p)
  for (i in 1:p) {
    for (j in 1:i) {
      if (i == j) {
        s[i, i] = b[k]
      } else {
        s[i, j] <- s[j, i] <- b[k]
      }
      k <- k + 1
    }
  }
  e <- eigen (s)
  f <- e$values
  if (min(f) < 0) {
    stop ("Negative eigenvalue, cannot proceed")
  }
  t <- e$vectors %*% diag (sqrt (f))
  v <- solve (t, b[-(1:l)])
  x <- g %*% t
  y <-
    h %*% (e$vectors %*% diag (1 / sqrt (f))) + matrix (v, m, p, byrow = TRUE)
  return (list (x = x, y = y))
}

unfoldals <- function (offdiag) {
  n <- nrow (offdiag)
  m <- ncol (offdiag)
  dd <- offdiag ^ 2
  delta <- matrix (0, n + m, n + m)
  delta[1:n, n + (1:m)] <- dd
  delta <- pmax(delta, t(delta))
  cc <-
    dd - outer (rowSums(dd) / m, colSums (dd) / n, "+") + sum(dd) / (n * m)
  sc <- svd (-cc / 2)
  lb <- diag (sqrt(sc$d))
  xold <- sc$u %*% lb
  yold <- sc$v %*% lb
  zold <- rbind (xold, yold)
  lold <- rowSums (zold ^ 2)
  dold <- outer (lold, lold, "+") - 2 * tcrossprod (zold)
  
}

teqbounds <- function (offdiag) {
  n <- nrow (offdiag)
  m <- ncol (offdiag)
  a <- matrix (0, n, n)
  b <- matrix (0, m, m)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      smin = Inf
      smax = -Inf
      for (k in 1:m) {
        smin = min (smin, offdiag[i, k] + offdiag[j, k])
        smax = max (smax, abs (offdiag[i, k] - offdiag[j, k]))
      }
      a[i, j] <- a[j, i] <- (smin + smax) / 2
    }
  }
  for (i in 2:m) {
    for (j in 1:(i - 1)) {
      smin = Inf
      smax = -Inf
      for (k in 1:n) {
        smin = min (smin, offdiag[k, i] + offdiag[k, j])
        smax = max (smax, abs (offdiag[k, i] - offdiag[k, j]))
      }
      b[i, j] <- b[j, i] <- (smin + smax) / 2
    }
  }
  return (list (a = a, b = b))
}
