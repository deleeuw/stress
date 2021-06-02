smacofEpsR <-
  function (w,
            delta,
            p,
            xold = smacofInitialR(delta, p),
            eold = .1,
            eopt = FALSE,
            xstop = FALSE,
            itmax = 1000,
            eps = 1e-10,
            verbose = TRUE) {
    labels = c("itel", "eiff", "sold", "snew")
    digits = c(4, 10, 10, 10)
    widths = c(6, 15, 15, 15)
    format = c("d", "f", "f", "f")
    n <- nrow (xold)
    j <- diag(n) - (1 / n)
    itel <- 1
    xold <- cbind(xold, eold * j)
    dold <- as.matrix (dist (xold))
    sold <- smacofLossR (dold, w, delta)
    bold <- smacofBmatR (dold, w, delta)
    vmat <- smacofVmatR (w)
    vinv <- ginv (vmat)
    repeat {
      xnew <- smacofGuttmanR (xold, bold, vinv)
      if (eopt) {
        enew <- sum (diag (xnew[,-(1:p)])) / (n - 1)
      } else {
        enew <- eold
      }
      xnew[, -(1:p)] <- enew * j
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
    return (list (
      x = xnew,
      d = dnew,
      b = bnew,
      e = enew,
      g = smacofGradientR(xnew, bnew, vmat),
      h = smacofHessianR(xnew, bnew, vmat, dnew, w, delta),
      s = snew,
      itel = itel
    ))
  }
