/*
qpmaj <-
  function (z,
            v = diag (length (z)),
            a = diff (diag (length(z))),
            b = ifelse (!is.null(a), rep(0, nrow(a)), NULL),
            c = NULL,
            d = ifelse (!is.null(c), rep(0, nrow(c)), NULL),
            h = NULL,
            itmax = 1000,
            eps = 1e-15,
            verbose = FALSE) {
    labs <- c("itel", "fold", "fnew")
    digs <- c(0, 6, 6)
    wids <- c(3, 10, 10)
    fors <- c("d", "f", "f")
    if (is.null(h)) {
      w <- v
      y <- z
      rsum <- 0
    } else {
      w <- crossprod(h, v %*% h)
      y <- drop (solve (w, crossprod (h, v %*% z)))
      rsum <- sum (z * (v %*% (z - h %*% y))) / 2
    }
    winv <- solve (w)
    if (!is.null(a)) {
      nin <- nrow(a)
      dualaa <- a %*% winv %*% t(a)
      feasa <- drop ((a %*% y) - b)
    }
    if (!is.null(c)) {
      neq <- nrow (c)
      dualcc <- c %*% winv %*% t(c)
      feasc <- drop((c %*% y) - d)
    }
    if ((!is.null(a)) && (!is.null(c))) {
      dualac <- a %*% winv %*% t(c)
      feas <- c(feasa, feasc)
      dual <- rbind(cbind(dualaa, dualac), cbind(t(dualac), dualcc))
    }
    if ((!is.null(a)) && (is.null(c))) {
      feas <- feasa
      dual <- dualaa
    }
    if ((is.null(a)) && (!is.null(c))) {
      feas <- feasc
      dual <- dualcc
    }
    vmax <- max(eigen(dual)$values)
    itel <- 1
    lold <- rep (0, nrow (dual))
    fold <-
      -(sum (lold * drop (dual %*% lold)) / 2 +
          sum (lold * feas))
    repeat {
      lnew <- lold - (drop(dual %*% lold) + feas) / vmax
      if (!is.null(a)) {
        lnew[1:nin] <- pmax(lnew[1:nin], 0)
      }
      fnew <-
        -(sum (lnew * drop (dual %*% lnew)) / 2 + sum (lnew * feas))
      if (verbose) {
        vals <- c(itel, fold, fnew)
        iterWrite (labs, vals, digs, wids, fors)
      }
      if ((itel == itmax) || ((fnew - fold) < eps)) {
        break
      }
      fold <- fnew
      lold <- lnew
      itel <- itel + 1
    }
    fdua <- fnew
    if ((!is.null(c) && (!is.null(a)))) {
      x <- y + drop (winv %*% drop(cbind(t(a), t(c)) %*% lnew))
      lb <- lnew[1:nin]
      mu <- lnew[-(1:nin)]
    }
    if ((is.null(c) && (!is.null(a)))) {
      x <- y + drop (winv %*% drop (t(a) %*% lnew))
      lb <- lnew
    }
    if ((!is.null(c) && (is.null(a)))) {
      x <- y + drop (winv %*% drop (t(c) %*% lnew))
      mu <- lnew
    }
    fprs <- sum ((x - y) * drop (w %*% (x - y))) / 2
    out <- list(x = x,
                fprimal = fprs,
                fdual = fdua)
    if (!is.null(h)) {
      out <-
        list.append(out, ftotal = fprs + rsum, predict = drop (h %*% x))
    }
    if (!is.null(a)) {
      out <-
        list.append(out, lambda = lb, inequalities = drop (a %*% x - b))
    }
    if (!is.null(c)) {
      out <- list.append(out, mu = mu, equations = drop (c %*% x - d))
    }
    return (list.append(out, itel = itel))
  }

*/

/*
    qpmaj minimizes (z-Hx)'V(z-Hx) over Ax≤b and Cx=d. The parts involving
    can be missing.
*/

void qpmaj(double *z, double *v, double *a, double *b, double *c, double *d,
           double *h, int n, int m, int p, int q, int itmax, double eps, bool verbose) {

    if (h == NULL) {
      (double *) w = calloc ((size_t) )

    }
}