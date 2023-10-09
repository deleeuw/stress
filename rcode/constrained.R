

pcircsmacof <-
  function (delta,
            w = wdef (nrow (delta)),
            p = 2,
            x = smacofInitialR (delta, p),
            itmax = 1000,
            eps = 1e-6,
            verbose = TRUE) {
    labels = c("itel", "sold", "snew")
    digits = c(4, 10, 10)
    widths = c(6, 15, 15)
    format = c("d", "f", "f")
    n <- nrow (x)
    p <- ncol (x)
    xold <- x / sqrt (rowSums (x ^ 2))
    dold <- as.matrix (dist (xold))
    v <- smacofVmatR (w)
    e <- max (eigen (v, only.values = TRUE)$values)
    vinv <- ginv(v)
    itel <- 1
    sold <- smacofLossR (dold, w, delta)
    repeat {
      b <- smacofBmatR (dold, w, delta)
      xgut <- smacofGuttmanR(xold, b, vinv)
      xtar <- xold + v %*% (xgut - xold) / e
      xlen <- sqrt (rowSums (xtar ^ 2))
      xrad <- mean (xlen)
      xnew <- (xtar / xlen) * xrad
      dnew <- as.matrix (dist(xnew))
      snew <- smacofLossR (dnew, w, delta)
        if (verbose) {
          values = c(itel, sold, snew)
          iterationWrite (labels, values, digits, width, format)
        }
        if (((sold - snew) < eps) || (itel == itmax)) {
          break
        }
        itel <- itel + 1
        xold <- xnew
        sold <- snew
      }
      return (list (
        x = xnew,
        d = dnew,
        stress = snew,
        radius = xrad,
        itel = itel
      ))
  }


pellipsmacof <-
  function (delta,
            w = wdef (nrow (delta)),
            p = 2,
            x = smacofInitialR (delta, p),
            itmax = 1000,
            eps = 1e-6,
            verbose = TRUE) {
    labels = c("itel", "sold", "smid", "snew")
    digits = c(4, 10, 10, 10)
    widths = c(6, 15, 15, 15)
    format = c("d", "f", "f", "f")
    n <- nrow (x)
    p <- ncol (x)
    yold <- x / sqrt (rowSums (x ^ 2))
    xlbd <- rep (1, p)
    xold <- yold %*% diag (xlbd)
    dold <- as.matrix (dist (xold))
    v <- smacofVmatR (w)
    e <- max (eigen (v, only.values = TRUE)$values)
    vinv <- ginv(v)
    itel <- 1
    sold <- smacofLoss(dold, w, delta)
    repeat {
      b <- smacofBmatR (dold, w, delta)
      xgut <- smacofGuttmanR(xold, b, vinv)
      for (s in 1:p) {
        xlbd[s] <-
          sum (xgut[, s] * (v %*% yold[, s])) / sum (yold[, s] * (v %*% yold[, s]))
      }
      xmid <- yold %*% diag (xlbd)
      dmid <- as.matrix (dist (xmid))
      smid <- sum (w * (delta - dmid) ^ 2) / 2
      mlbd <- max (xlbd ^ 2)
      ytar <-
        yold + v %*% ((xgut %*% diag (1 / xlbd)) - yold) %*% diag (xlbd ^ 2) / (e * mlbd)
      ylen <- sqrt (rowSums (ytar ^ 2))
      ynew <- ytar / ylen
      xnew <- ynew %*% diag (xlbd)
      dnew <- as.matrix (dist(xnew))
      snew <- sum (w * (delta - dnew) ^ 2) / 2
      if (verbose) {
        values = c(itel, sold, smid, snew)
        iterationWrite (labels, values, digits, width, format)
      }
      if (((sold - snew) < eps) || (itel == itmax)) {
        break
      }
      itel <- itel + 1
      xold <- xnew
      yold <- ynew
      sold <- snew
    }
    return (list (
      x = xnew,
      d = dnew,
      stress = snew,
      axes = xlbd,
      itel = itel
    ))
  }

dcircsmacof <-
  function (delta,
            w = wdef (nrow (delta)),
            p = 2,
            x = smacofInitialR (delta, p),
            pen = 1,
            itmax = 1000,
            eps = 1e-6,
            verbose = TRUE) {
    labels = c("itel", "sold", "smid", "snew")
    digits = c(4, 10, 10, 10)
    widths = c(6, 15, 15, 15)
    format = c("d", "f", "f", "f")
    n <- nrow (x)
    xold <-
      rbind (0, x / sqrt (rowSums(x ^ 2)))
    dold <- as.matrix (dist (xold))
    w <- rbind (pen, cbind (pen, w))
    delta <- rbind (1, cbind (1, delta))
    w[1, 1] <- delta [1, 1] <- 0
    v <- smacofVmatR (w)
    vinv <- ginv(v)
    itel <- 1
    sold <- smacofLossR(dold, w, delta)
    repeat {
      b <- smacofBmatR(dold, w, delta)
      xnew <- smacofGuttmanR(xold, b, vinv)
      dnew <- as.matrix (dist (xnew))
      smid <- smacofLossR(dnew, w, delta)
      a <- sum (dnew[1,]) / n
      delta[1,] <- delta[, 1] <- a
      delta[1, 1] <- 0
      snew <- smacofLossR(dnew, w, delta)
      if (verbose) {
        values = c(itel, sold, smid, snew)
        iterationWrite (labels, values, digits, width, format)
      }
      if (((sold - snew) < eps) || (itel == itmax)) {
        break
      }
      itel <- itel + 1
      xold <- xnew
      sold <- snew
    }
    return (list (
      x = xnew,
      d = dnew,
      a = a,
      stress = snew,
      itel = itel
    ))
  }

dellipsmacof <-
  function (delta,
            w = wdef (nrow (delta)),
            p = 2,
            x = smacofInitialR (delta, p),
            pen = 1,
            itmax = 1000,
            eps = 1e-6,
            verbose = TRUE) {
    labels = c("itel", "sold", "smid", "snew")
    digits = c(4, 10, 10, 10)
    widths = c(6, 15, 15, 15)
    format = c("d", "f", "f", "f")
    n <- nrow (x)
    set.seed(12345)
    focal <- rnorm(2)
    xold <-
      rbind (focal, -focal, x / sqrt (rowSums(x ^ 2)))
    dold <- as.matrix (dist (xold))
    w <- rbind (pen, pen, cbind (pen, pen, w))
    delta <- rbind (1, 1, cbind (1, 1, delta))
    w[1:2, 1:2] <- delta [1:2, 1:2] <- 0
    v <- smacofVmatR (w)
    vinv <- ginv(v)
    itel <- 1
    sold <- smacofLossR(dold, w, delta)
    repeat {
      b <- smacofBmatR (dold, w, delta)
      xnew <- smacofGuttmanR (xold, b, vinv)
      dnew <- as.matrix (dist (xnew))
      smid <- smacofLossR(dnew, w, delta)
      dsub <- dnew[1:2, 2 + (1:n)]
      asub <- sum (dsub) / (2 * n)
      dsub <- ccen (dsub) + asub
      delta[1:2, 2 + (1:n)] <- dsub
      delta[2 + (1:n), 1:2] <- t(dsub)
      snew <- smacofLossR(dnew, w, delta)
      if (verbose) {
        values = c(itel, sold, smid, snew)
        iterationWrite (labels, values, digits, width, format)
      }
      if (((sold - snew) < eps) || (itel == itmax)) {
        break
      }
      itel <- itel + 1
      xold <- xnew
      sold <- snew
    }
    return (list (
      pen = pen,
      x = xnew,
      d = dnew,
      stress = snew,
      itel = itel
    ))
  }
