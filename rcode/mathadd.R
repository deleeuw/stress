

lsuw <- function (y,
                  w,
                  proj,
                  xold = rep (1, length (y)),
                  v = max (eigen (w)$values)  * diag (length (y)),
                  itmax = 100,
                  eps = 1e-6,
                  verbose = FALSE,
                  add = 1e-6) {
  f <- function (x, y, w) {
    return (sum ((x - y) * w %*% (x - y)))
  }
  n <- length (y)
  labels <- c("itel", "fold", "fnew")
  digits <- c(0, 6, 6)
  widths <- c(3, 10, 10)
  formats <- c("d", "f", "f")
  fold <- f (xold, y, w)
  itel <- 1
  repeat {
    u <- drop (solve (v, w %*% (xold - y)))
    xnew <- proj (xold - u, v)
    fnew <- f (xnew, y, w)
    if (verbose) {
      values <- c(itel, fold, fnew)
      iterWrite (labels, values, digits, widths, formats)
    }
    if ((itel == itmax) || ((fold - fnew) < eps)) {
      break
    }
    fold <- fnew
    xold <- xnew
    itel <- itel + 1
  }
  return (list (x = xnew, f = fnew, itel = itel))
}

projeq <- function (x, v) {
  s <- sum (v)
  h <- sum (x * rowSums (v))
  return (rep (h / s, length (x)))
}

projplus <- function (x, v) {
  if (!all(v == diag(diag(v)))) {
    stop ("V must be diagonal")
  }
  if (min (diag (v)) < 0) {
    stop ("V must be positive semidefinite")
  }
  return (pmax(x, 0))
}

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


checkIncreasing <- function (innerknots, lowend, highend) {
  h <- .C(
    "checkIncreasing",
    as.double (innerknots),
    as.double (lowend),
    as.double (highend),
    as.integer (length (innerknots)),
    fail = as.integer (0)
  )
  return (h$fail)
}

extendPartition <-
  function (innerknots,
            multiplicities,
            order,
            lowend,
            highend) {
    ninner <- length (innerknots)
    kk <- sum(multiplicities)
    nextended <- kk + 2 * order
    if (max (multiplicities) > order)
      stop ("multiplicities too high")
    if (min (multiplicities) < 1)
      stop ("multiplicities too low")
    if (checkIncreasing (innerknots, lowend, highend))
      stop ("knot sequence not increasing")
    h <-
      .C(
        "extendPartition",
        as.double (innerknots),
        as.integer (multiplicities),
        as.integer (order),
        as.integer (ninner),
        as.double (lowend),
        as.double (highend),
        knots = as.double (rep (0, nextended))
      )
    return (h)
  }

bisect <-
  function (x,
            knots,
            lowindex = 1,
            highindex = length (knots)) {
    h <- .C(
      "bisect",
      as.double (x),
      as.double (knots),
      as.integer (lowindex),
      as.integer (highindex),
      index = as.integer (0)
    )
    return (h$index)
  }

bsplines <- function (x, knots, order) {
  if ((x > knots[length(knots)]) || (x < knots[1]))
    stop ("argument out of range")
  h <- .C(
    "bsplines",
    as.double (x),
    as.double (knots),
    as.integer (order),
    as.integer (length (knots)),
    index = as.integer (0),
    q = as.double (rep(0, order))
  )
  return (list (q = h$q, index = h$ind))
}


bsplineBasis <- function (x, knots, order) {
  n <- length (x)
  k <- length (knots)
  m <- k - order
  result <- rep (0, n * m)
  h <- .C(
    "bsplineBasis",
    as.double (x),
    as.double (knots),
    as.integer (order),
    as.integer (k),
    as.integer (n),
    result = as.double (result)
  )
  return (matrix (h$result, n, m))
}

isplineBasis <-  function (x, knots, order) {
  n <- length (x)
  k <- length (knots)
  m <- k - order
  result <- rep (0, n * m)
  h <- .C(
    "isplineBasis",
    as.double (x),
    as.double (knots),
    as.integer (order),
    as.integer (k),
    as.integer (n),
    result = as.double (result)
  )
  return (matrix (h$result, n, m))
}
