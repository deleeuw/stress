dommy <- function () {
  ya <- matrix(c(1, -1, 1, -1, 0, 1, 1, -1, -1, 0), 5, 2)
  yy <- seq (0, 2 * pi, length = 6)[1:5]
  yb <- cbind (sin (yy), cos (yy))
  ya <- apply(ya, 2, function (x)
    x - mean(x))
  yb <- apply(yb, 2, function (x)
    x - mean(x))
  y1 <- ya / sqrt (5 * sum (ya ^ 2))
  y2 <- yb / sqrt (5 * sum (yb ^ 2))
  deq <- as.dist (1 - diag(5))
  deq <- deq / sqrt (sum (deq ^ 2))
  y1 <- sum (dist (y1) * deq) * y1
  y2 <- sum (dist (y2) * deq) * y2
}

twostress <- function (deq, y1, y2, a, b) {
  d <- dist (a * y1 + b * y2)
  eta2 <- sum (d ^ 2)
  rho <- sum (d * deq)
  stress <- 1 - 2 * rho + eta2
  return(list(
    eta2 = eta2,
    rho = rho,
    stress = stress
  ))
}

zeroes <- function (y1, y2) {
  n <- nrow (y1)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      yy <- rbind (y1[i, ] - y1[j, ], y2[i, ] - y2[j, ])
      ee <- eigen(tcrossprod(yy))$values
      print (c(i, j, ee))
    }
  }
}

pairme <- function (x, y) {
  n <- nrow(x)
  m <- ncol(x)
  z <- matrix (0, 2, m)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      z[1,] <- x[i,] - x[j,]
      z[2,] <- y[i,] - y[j,]
      s <- svd (z)
      a <- s$d[2]
      b <- s$u[, 2]
      if (a < 1e-10) {
        cat (
          formatC(
            i,
            format = "d",
            digits = 2,
            width = 4
          ),
          formatC(
            j,
            format = "d",
            digits = 2,
            width = 4
          ),
          formatC(
            c(a, b),
            format = "f",
            digits = 6,
            width = 10
          ),
          "\n"
        )
      }
    }
  }
}

dummy <- function () {
  set.seed(12345)
  x <- matrix(rnorm(10), 5, 2)
  x <- apply (x, 2, function(x)
    x - mean(x))
  delta <- dist(x)
  d <- dist (x)
  eps <- (-500:500) / 100
  sy <- rep (0, 1001)
  plot (
    0,
    0,
    xlim = c(-5, 5),
    ylim = c(0, 20),
    xlab = "epsilon",
    ylab = "stress",
    type = "n"
  )
  for (i in 1:5) {
    for (j in 1:2) {
      for (k in 1:1001) {
        y <- x
        y[i, j] <- x[i, j] + eps[k]
        dy <- dist (y)
        sy[k] <- sum ((delta - dy) ^ 2)
      }
      lines (eps, sy, lwd = 2, col = "RED")
    }
  }
}


bmat2 <- function (a, b, x, y, delta) {
  bm <- matrix (0, 2, 2)
  hm <- matrix (0, 2, 2)
  z <- c(a, b)
  for (i in 1:4) {
    for (j in 1:4) {
      if (i == j)
        next
      uij <- uu (i, j, x, y)
      uz <- drop (uij %*% z)
      dij <- sqrt (sum (uij * outer (z, z)))
      bm <- bm + (delta[i, j] / dij) * uij
      hm <-
        hm + (delta[i, j] / dij) * (uij - outer (uz, uz) / sum (z * uz))
    }
  }
  return (list (b = bm, h = hm))
}

stress2 <- function (a, b, x, y, delta) {
  z <- c (a, b)
  bm <- bmat2 (a, b, x, y, delta)$b
  return (1 + sum(z ^ 2) / 2 - sum (z * bm %*% z))
}

rho2 <- function (a, b, x, y, delta) {
  z <- c (a, b)
  bm <- bmat2 (a, b, x, y, delta)$b
  return (sum (z * bm %*% z))
}

vv <- function (i, j, x, y) {
  a <- matrix (0, 2, 2)
  a[1, 1] <- sum ((x[i,] - x[j, ]) ^ 2)
  a[2, 2] <- sum ((y[i,] - y[j, ]) ^ 2)
  a[1, 2] <- a[2, 1] <- sum ((x[i,] - x[j, ]) * (y[i,] - y[j,]))
  return (a)
}

uu <- function (i, j, x, y) {
  n <- nrow (x)
  asum <-
    2 * n * matrix (c (sum(x ^ 2), sum (x * y), sum (x * y), sum (y ^ 2)), 2, 2)
  csum <- solve (chol (asum))
  return (t(csum) %*% vv (i, j, x, y) %*% csum)
}

smacof2 <-
  function (a,
            b,
            x,
            y,
            delta,
            eps = 1e-10,
            itmax = 1000,
            verbose = TRUE) {
    zold <- c(a, b)
    bold <- bmat2 (a, b, x, y, delta)$b
    fold <- 1 + sum(zold ^ 2) / 2 - sum (zold * bold %*% zold)
    itel <- 1
    repeat {
      znew <- drop (bold %*% zold)
      bhmt <- bmat2 (znew[1], znew[2], x, y, delta)
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
    return (
      list (
        stress = fnew,
        theta = znew,
        itel = itel,
        b = bnew,
        g = znew - bnew %*% znew,
        h = diag(2) - bhmt$h
      )
    )
  }


newton2 <-
  function (a,
            b,
            x,
            y,
            delta,
            eps = 1e-10,
            itmax = 1000,
            verbose = TRUE) {
    zold <- c(a, b)
    bhmt <- bmat2 (a, b, x, y, delta)
    bold <- bhmt$b
    hold <- diag(2) - bhmt$h
    fold <- 1 + sum(zold ^ 2) / 2 - sum (zold * bold %*% zold)
    itel <- 1
    repeat {
      znew <- drop (solve (hold, bold %*% zold))
      bhmt <- bmat2 (znew[1], znew[2], x, y, delta)
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
    return (list (
      stress = fnew,
      theta = znew,
      itel = itel,
      b = bnew,
      g = znew - bnew %*% znew,
      h = hnew
    ))
  }
