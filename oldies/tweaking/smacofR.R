library(MASS)

smacofLossR <- function (d, w, delta) {
  return (sum (w * (delta - d) ^ 2) / 4)
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
  return ((v - b) %*% x)
}

smacofHmatR <- function (x, b, v, d, w, delta) {
  n <- nrow (x)
  p <- ncol (x)
  r <- n * p
  h <- matrix (0, r, r)
  dd <- ifelse (d == 0, 0, 1 / d)
  cc <- w * delta * (dd ^ 3)
  for (s in 1:p) {
    for (t in 1:s) {
      cst <- matrix (0, n, n)
      for (i in 1:n) {
        for (j in 1:n) {
          cst[i, j] <- cc[i, j] * (x[i, s] - x[j, s]) * (x[i, t] - x[j, t])
        }
      }
      cst <- -cst
      diag(cst) <- -rowSums(cst)
      if (s == t) {
        h[(s - 1) * n + 1:n, (s - 1) * n + 1:n] <- b - cst
      } else {
        h[(s - 1) * n + 1:n, (t - 1) * n + 1:n] <- -cst
        h[(t - 1) * n + 1:n, (s - 1) * n + 1:n] <- -cst
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
    h[(s - 1) * n + 1:n, (s - 1) * n + 1:n] <-
      h[(s - 1) * n + 1:n, (s - 1) * n + 1:n] + v
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

smacofR <-
  function (w,
            delta,
            p,
            xold = smacofInitialR(delta, p),
            itmax = 100,
            eps = 1e-10,
            verbose = FALSE) {
    itel <- 1
    dold <- as.matrix (dist (xold))
    sold <- smacofLossR (dold, w, delta)
    bold <- smacofBmatR (dold, w, delta)
    vinv <- ginv (smacofVmatR (w))
    repeat {
      xnew <- smacofGuttmanR (xold, bold, vinv)
      eiff <- max (abs (xold - xnew))
      dnew <- as.matrix (dist (xnew))
      bnew <- smacofBmatR (dnew, w, delta)
      snew <- smacofLossR (dnew, w, delta)
      if (verbose) {
        cat(
          "itel ",
          formatC(itel, width = 4, format = "d"),
          "eiff ",
          formatC(
            eiff,
            width = 15,
            digits = 10,
            format = "f"
          ),
          "sold ",
          formatC(
            sold,
            width = 15,
            digits = 10,
            format = "f"
          ),
          "snew ",
          formatC(
            snew,
            width = 15,
            digits = 10,
            format = "f"
          ),
          "\n"
        )
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
    return (list (
      x = xnew,
      d = dnew,
      b = bnew,
      s = snew,
      itel = itel
    ))
  }
