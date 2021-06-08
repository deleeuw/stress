
smacofLossRC <- function (d, w, delta) {
  m <- length (d)
  h <-
    .C(
      "smacofLossC",
      as.double (d),
      as.double (w),
      as.double (delta),
      as.integer (m),
      stress = as.double (0)
    )
  return (h$stress)
}

smacofDistRC <- function (x, p) {
  n <- round (length (x) / p)
  m <- n * (n - 1) / 2
  h <-
    .C("smacofDistC",
       as.double (x),
       as.integer (n),
       as.integer(p),
       dist = as.double (rep(0, m)))
  return (h$dist)
}

smacofInitialRC <- function (delta, p) {
  m <- length (delta)
  n <- round ((1 + sqrt (1 + 8 * m)) / 2)
  r <- n * (n + 1) / 2
  h <-
    .C(
      "smacofInitialC",
      as.double(delta),
      as.integer(n),
      as.integer(p),
      as.double (rep(0, n)),
      as.double (rep(0, r)),
      as.double (rep(0, n * n)),
      as.double (rep (0, n)),
      x = as.double (rep (0, n * p))
    )
  return (h$x)
}

smacofBmatRC <- function (d, w, delta) {
  m <- length (w)
  n <- round ((1 + sqrt (1 + 8 * m)) / 2)
  r <- n * (n + 1) / 2
  h <-
    .C(
      "smacofBmatC",
      as.double(d),
      as.double(w),
      as.double(delta),
      as.integer (n),
      bmat = as.double(rep(0, r))
    )
  return (h$bmat)
}

smacofVmatRC <- function (w) {
  m <- length (w)
  n <- round ((1 + sqrt (1 + 8 * m)) / 2)
  r <- n * (n + 1) / 2
  h <-
    .C("smacofVmatC",
       as.double(w),
       as.integer (n),
       vmat = as.double(rep(0, r)))
  return (h$vmat)
}

smacofGuttmanRC <- function (x, b, vinv) {
  m <- length (b)
  n <- round ((sqrt (1 + 8 * m) - 1) / 2)
  r <- length (x)
  p <- round (r / n)
  h <-
    .C(
      "smacofGuttmanC",
      as.double(x),
      as.double (b),
      as.double (vinv),
      as.integer(n),
      as.integer(p),
      as.double (rep(0, r)),
      y = as.double (rep(0, r))
    )
  return (h$y)
}

smacofGradientRC <- function (x, b, v) {
  m <- length (b)
  n <- round ((sqrt (1 + 8 * m) - 1) / 2)
  r <- length (x)
  p <- round (r / n)
  h <-
    .C(
      "smacofGradientC",
      as.double(x),
      as.double (b),
      as.double (v),
      as.integer(n),
      as.integer(p),
      y = as.double (rep(0, r))
    )
  return (h$y)
}

smacofHmatRC <- function (x, b, v, d, w, delta) {
  m <- length (w)
  n <- round ((1 + sqrt (1 + 8 * m)) / 2)
  q <- n * (n + 1) / 2
  r <- length (x)
  p <- round (r / n)
  h <-
    .C(
      "smacofHmatC",
      as.double(x),
      as.double (b),
      as.double (v),
      as.double (w),
      as.double (delta),
      as.double (d),
      as.integer(n),
      as.integer(p),
      as.double (rep (0, q)),
      hmat = as.double (rep (0, r * (r + 1) / 2))
    )
  return (h$hmat)
}

smacofHessianRC <- function (x, b, v, d, w, delta) {
  m <- length (w)
  n <- round ((1 + sqrt (1 + 8 * m)) / 2)
  q <- n * (n + 1) / 2
  r <- length (x)
  p <- round (r / n)
  h <-
    .C(
      "smacofHessianC",
      as.double(x),
      as.double (b),
      as.double (v),
      as.double (w),
      as.double (delta),
      as.double (d),
      as.integer(n),
      as.integer(p),
      as.double (rep (0, q)),
      hmat = as.double (rep (0, r * (r + 1) / 2))
    )
  return (h$hmat)
}

smacofUpdateRC <- function (x, w, delta, b, vinv) {
  m <- length (b)
  n <- round ((sqrt (1 + 8 * m) - 1) / 2)
  r <- length (x)
  p <- round (r / n)
  q <- n * (n - 1) / 2
  h <-
    .C(
      "smacofUpdateC",
      as.double (x),
      as.double (w),
      as.double (delta),
      as.double (vinv),
      as.integer (n),
      as.integer (p),
      dnew = as.double (rep(0, q)),
      bnew = as.double (b),
      work = as.double (rep (0, n * p)),
      snew = as.double (0),
      xnew = as.double (rep (0, n * p))
    )
  return (h)
}

smacofUpdateRC64 <- function (x, w, delta, b, vinv) {
  m <- length (b)
  n <- round ((sqrt (1 + 8 * m) - 1) / 2)
  r <- length (x)
  p <- round (r / n)
  q <- n * (n - 1) / 2
  h <-
    .C64(
      "smacofUpdateC",
      SIGNATURE = c(rep("double", 4), rep("integer", 2), rep("double", 5)),
      x = x,
      w = w,
      delta = delta,
      vinv = vinv,
      n = n,
      p = p,
      dnew = numeric_dc (q),
      bnew = as.double (b),
      work = numeric_dc (n * p),
      snew = numeric_dc (1),
      xnew = numeric_dc (n * p),
      INTENT = c(rep("r", 6), "w", "rw", rep ("w", 3)),
      NAOK = TRUE
    )
  return (h)
}

smacofRC <-
  function (w,
            delta,
            p,
            xold = smacofInitialRC(delta, p),
            itmax = 100,
            eps = 1e-10,
            verbose = FALSE) {
    itel <- 1
    dold <- smacofDistRC (xold, p)
    sold <- smacofLossRC (dold, w, delta)
    bold <- smacofBmatRC (dold, w, delta)
    vinv <- mpowerRC (smacofVmatRC (w), -1)
    repeat {
      xnew <- smacofGuttmanRC (xold, bold, vinv)
      eiff <- max (abs (xold - xnew))
      dnew <- smacofDistRC (xnew, p)
      bnew <- smacofBmatRC (dnew, w, delta)
      snew <- smacofLossRC (dnew, w, delta)
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

smacofRCU <-
  function (w,
            delta,
            p,
            xold = smacofInitialRC(delta, p),
            itmax = 100,
            eps = 1e-10,
            verbose = FALSE) {
    itel <- 1
    dold <- smacofDistRC (xold, p)
    sold <- smacofLossRC (dold, w, delta)
    bold <- smacofBmatRC (dold, w, delta)
    vinv <- mpowerRC (smacofVmatRC (w), -1)
    repeat {
      h <- smacofUpdateRC(xold, w, delta, bold, vinv)
      xnew <- h$xnew
      dnew <- h$dnew
      snew <- h$snew
      bnew <- h$bnew
      eiff <- max (abs (xold - xnew))
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

smacofRCU64 <-
  function (w,
            delta,
            p,
            xold = smacofInitialRC(delta, p),
            itmax = 100,
            eps = 1e-10,
            verbose = FALSE) {
    itel <- 1
    dold <- smacofDistRC (xold, p)
    sold <- smacofLossRC (dold, w, delta)
    bold <- smacofBmatRC (dold, w, delta)
    vinv <- mpowerRC (smacofVmatRC (w), -1)
    repeat {
      h <- smacofUpdateRC64(xold, w, delta, bold, vinv)
      xnew <- h$xnew
      dnew <- h$dnew
      snew <- h$snew
      bnew <- h$bnew
      eiff <- max (abs (xold - xnew))
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
