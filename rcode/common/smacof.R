

smacofNormW <- function(w) {
  return (w / sum(w))
}

smacofNormDelta <- function(w, delta) {
  return (delta / sqrt(sum (w * delta ^ 2)))
}

smacofNormXD <- function(w, delta, xold) {
  x <- apply (xold, 2, function(x)
    x - mean(x))
  d <- as.matrix(dist(x))
  s <- sum (w * delta * d) / sum (w * d ^ 2)
  return (list(x = x * s, d = d * s))
}

smacofLossR <- function (d, w, delta) {
  return (sum (w * (delta - d) ^ 2) / 2)
}

smacofBmatR <- function (d, w, delta) {
  dd <- ifelse (d == 0, 0, 1 / d)
  b <- -dd * w * delta
  diag (b) <- -rowSums (b)
  return(b)
}

smacofVmatR <- function (w) {
  v <- -w
  diag(v) <- -rowSums(v)
  return (v)
}

smacofGuttmanR <- function (x, b, vinv) {
  return (vinv %*% b %*% x)
}

smacofGradientR <- function (x, b, v) {
  return (2 * ((v - b) %*% x))
}

smacofHmatR <- function (x, b, v, d, w, delta) {
  n <- nrow (x)
  p <- ncol (x)
  r <- n * p
  h <- matrix (0, r, r)
  dd <- ifelse (d == 0, 0, 1 / d)
  cc <- w * delta * (dd ^ 3)
  for (s in 1:p) {
    ns <- (s - 1) * n + 1:n
    for (t in 1:s) {
      nt <- (t - 1) * n + 1:n
      cst <- matrix (0, n, n)
      for (i in 1:n) {
        for (j in 1:n) {
          cst[i, j] <- cc[i, j] * (x[i, s] - x[j, s]) * (x[i, t] - x[j, t])
        }
      }
      cst <- -cst
      diag(cst) <- -rowSums(cst)
      if (s == t) {
        h[ns, ns] <- b - cst
      } else {
        h[ns, nt] <- -cst
        h[nt, ns] <- -cst
      }
    }
  }
  return (h)
}

smacofHessianR <- function (x, b, v, d, w, delta) {
  n <- nrow (x)
  p <- ncol (x)
  h <- -smacofHmatR (x, b, v, d, w, delta)
  for (s in 1:p) {
    nn <- (s - 1) * n + 1:n
    h[nn, nn] <- h[nn, nn] + v
  }
  return(h)
}

smacofDerGuttmanR <- function(x, b, vinv, d, w, delta) {
  n <- nrow (x)
  p <- ncol (x)
  h <- smacofHmatR (x, b, v, d, w, delta)
  for (s in 1:p) {
    ns <- (s - 1) * n + 1:n
    for (t in 1:s) {
      nt <- (t - 1) * n + 1:n
      h[ns, nt] <- vinv %*% h[ns, nt]
    }
  }
  return(h)
}

smacofInitialR <- function (delta, p) {
  n <- nrow(delta)
  delta <- delta ^ 2
  rw <- rowSums (delta) / n
  sw <- sum (delta) / (n ^ 2)
  h <- -(delta - outer (rw, rw, "+") + sw) / 2
  e <- eigen (h)
  ea <- e$values
  ev <- e$vector
  ea <- ifelse (ea > 0, sqrt (abs(ea)), 0)[1:p]
  return (ev[, 1:p] %*% diag (ea))
}

smacofRandomStart <- function (w, delta, n, p) {
  x <- matrix(rnorm(n * p), n, p)
  x <- apply (x, 2, function(x)
    x - mean(x))
  d <- as.matrix(dist (x))
  a <- sum (w * delta * d) / sum (w * d ^ 2)
  return (a * x)
}

smacofVinvR <- function (v) {
  e <- 1 / nrow(v)
  return (solve (v + e) - e)
}

smacofR <-
  function (w,
            delta,
            p,
            xold = smacofInitialR(delta, p),
            xstop = FALSE,
            itmax = 1000,
            eps = 1e-10,
            verbose = TRUE) {
    labels = c("itel", "eiff", "sold", "snew")
    digits = c(4, 10, 10, 10)
    widths = c(6, 15, 15, 15)
    format = c("d", "f", "f", "f")
    n <- dim(delta)[1]
    itel <- 1
    w <- smacofNormW(w)
    delta <- smacofNormDelta(w, delta)
    xdold <- smacofNormXD(w, delta, xold)
    xold <- xdold$x
    dold <- xdold$d
    sold <- smacofLossR (dold, w, delta)
    bold <- smacofBmatR (dold, w, delta)
    vmat <- smacofVmatR (w)
    vinv <- smacofVinvR (vmat)
    repeat {
      xnew <- smacofGuttmanR (xold, bold, vinv)
      dnew <- as.matrix (dist (xnew))
      bnew <- smacofBmatR (dnew, w, delta)
      snew <- smacofLossR (dnew, w, delta)
      if (xstop) {
        eiff <- max (abs (xold - xnew))
      } else {
        eiff <- sold - snew
      }
      if (verbose) {
        values = c(itel, eiff, sold, snew)
        iterationWrite (labels, values, digits, widths, format)
      }
      if ((eiff < eps) || (itel == itmax)) {
        break
      }
      itel <- itel + 1
      xold <- xnew
      bold <- bnew
      dold <- dnew
      sold <- snew
    }
    return (
      list (
        x = xnew,
        d = dnew,
        b = bnew,
        g = smacofGradientR(xnew, bnew, vmat),
        h = smacofHessianR(xnew, bnew, vmat, dnew, w, delta),
        s = snew,
        itel = itel
      )
    )
  }
