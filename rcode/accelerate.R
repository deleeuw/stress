library(MASS)

smacofRelaxR <-
  function (w,
            delta,
            p,
            xold = smacofInitialR(delta, p),
            strategy = 1,
            renormalize = 0,
            xstop = TRUE,
            itmax = 1000,
            eps = 1e-10,
            verbose = FALSE) {
    labels = c("itel", "eiff", "sold", "snew")
    digits = c(4, 10, 10, 10)
    widths = c(6, 15, 15, 15)
    format = c("d", "f", "f", "f")
    itel <- 1
    dold <- as.matrix (dist (xold))
    sold <- smacofLossR (dold, w, delta)
    bold <- smacofBmatR (dold, w, delta)
    vmat <- smacofVmatR (w)
    vinv <- ginv (vmat)
    repeat {
      znew <- smacofGuttmanR (xold, bold, vinv)
      if (strategy == 1) {
        xnew = znew;
      }
      if (strategy == 2) {
        xnew = 2 * znew - xold
      }
      if (strategy == 3) {
        zmid <- znew
        dmid <- as.matrix (dist (zmid))
        bmid <- smacofBmatR(dmid, w, delta)
        xnew <- smacofGuttmanR(zmid, bmid, vinv)
      }
      if (strategy == 4) {
        zmid <- 2 * znew - xold
        dmid <- as.matrix (dist (zmid))
        bmid <- smacofBmatR(dmid, w, delta)
        xnew <- smacofGuttmanR(zmid, bmid, vinv)
      }
      if (strategy == 5) {
        zmid <- 2 * znew - xold
        dmid <- as.matrix (dist (zmid))
        bmid <- smacofBmatR(dmid, w, delta)
        znxt <- smacofGuttmanR(zmid, bmid, vinv)
        xnew <- 2 * znxt - zmid
      }
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
        iterationWrite (labels, values, digits, width, format)
      }
       if ((eiff < eps) || (itel == itmax)) {
        break
      }
      itel <- itel + 1
      if (itel == itmax) {
        xaux <- xold
      }
      xold <- xnew
      bold <- bnew
      dold <- dnew
      sold <- snew
    }
    if (renormalize == 1) {
      lbd <- sum (dnew * delta) / sum (dnew ^ 2)
      xnew <- xnew * lbd
      dnew <- dnew * lbd
      snew <- smacofLossR (dnew, w, delta)
    }
    if (renormalize == 2) {
      xnew <- smacofGuttmanR(xnew, bnew, vinv)
      dnew <- as.matrix(dist(xnew))
      snew <- smacofLossR(dnew, w, delta)
    }
    return (list (
      x = xnew,
      d = dnew,
      b = bnew,
      g = smacofGradientR(xnew, bnew, vmat),
      h = smacofHessianR(xnew, bnew, vmat, dnew, w, delta),
      s = snew,
      itel = itel
    ))
  }
