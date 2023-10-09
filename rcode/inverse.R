
imdsSolver <- function (tau) {
  radius <- sqrt (((tau - 3) ^ 2) / 3)
  center <- c(tau, tau) / 3
  a <- matrix (c(1, .5, .5, 1), 2, 2)
  draw_ellipse (
    center,
    radius,
    a,
    col = "RED",
    lwd = 2,
    xlim = c(0, tau),
    ylim = c(0, tau),
    xlab = "alpha",
    ylab = "beta"
  )
  lines (matrix(c(0, tau, tau, 0), 2, 2), col = "BLUE", lwd = 2)
  abline(h = 0)
  abline(v = 0)
}

bs <- function () {
  z <- matrix(c(1, 1, 1, 1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, -1, 1), 4 , 4) / 2
  a <- as.list (1:6)
  k <- 1
  for (i in 1:3) {
    for (j in (i + 1):4) {
      a[[k]] <- crossprod (z, aijn(i, j, 4) %*% z)
      k <- k + 1
    }
  }
  return(a)
}

imdsChecker <- function (a) {
  aa <- bs()
  bb <- matrix(0, 4, 4)
  for (i in 1:6) {
    bb <- bb + a[i] * aa[[i]]
  }
  return (bb)
}

inverseMDS <- function (x) {
  n <- nrow (x)
  m <- ncol (x)
  x <- apply (x, 2, function (y)
    y - mean (y))
  nm <- n - (m + 1)
  kk <- cbind (1, x, matrix (rnorm (n * nm), n , nm))
  kperp <- as.matrix (qr.Q (qr (kk))[,-(1:(m + 1))])
  dd <- as.matrix (dist (x))
  k <- 1
  base <- matrix (0, n * (n - 1) / 2, nm * (nm + 1) / 2)
  for (i in 1:nm) {
    for (j in 1:i) {
      oo <- outer (kperp[, i], kperp[, j])
      if (j != i) {
        oo <- oo + t(oo)
      }
      base[, k] <- lower_triangle (dd + (1 - oo))
      k <- k + 1
      print (c(i, j, k))
    }
  }
  return (base = cbind (lower_triangle (dd), base))
}

inversePlus <- function (base, affine = TRUE) {
  if (affine) {
    hrep <- makeH (
      a1 = d2q (-base),
      b1 = d2q (rep (0, nrow (base))),
      a2 = d2q (rep (1, ncol (base))),
      b2 = d2q (1)
    )
  } else {
    hrep <- makeH (a1 = d2q (-base), b1 = d2q (rep (0, nrow (base))))
  }
  vrep <- scdd (hrep)
  hrep <- q2d (hrep)
  vrep <- q2d (vrep$output)
  pr <- tcrossprod (hrep[, -c(1, 2)], vrep[, -c(1, 2)])[-1, ]
  return (list (
    base = base,
    hrep = hrep,
    vrep = vrep,
    pr = pr
  ))
}

twoPoints <- function (x, y, w = 1 - diag (nrow (x))) {
  dx <- lower_triangle (as.matrix (dist (x)))
  dy <- lower_triangle (as.matrix (dist (y)))
  w <- lower_triangle (w)
  gx <- makeG (x)
  gy <- makeG (y)
  hx <- (dx / w) * gx
  hy <- (dy / w) * gy
  lxy <- lm.fit (cbind (hx,-hy), dx - dy)
  lxx <- lxy$coefficients[1:ncol(hx)]
  lyy <- lxy$coefficients[-(1:ncol(hx))]
  return (list(
    delta1 = dx - hx %*% lxx,
    delta2 = dy - hy %*% lyy,
    res = sum (abs(lxy$residuals)),
    rank = lxy$rank
  ))
}

second_partials_stress <-
  function (x, delta, w = wdef (nrow (x))) {
    n <- nrow (x)
    p <- ncol (x)
    d <- as.matrix (dist (x))
    fac <- (w * delta) / (d + diag (n))
    dd <- d * d
    v <- smacofVmatR (w)
    deri <- direct_sum (repList (v, p))
    xx <- as.vector (x)
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        aa <- direct_sum (repList (aijn (i, j, n), p))
        ax <- drop (aa %*% xx)
        deri <- deri - fac[i, j] * (aa - outer (ax, ax) / dd[i, j])
      }
    }
    return (deri)
  }

second_partials_numerical <-
  function (x, delta, w = wdef (nrow (x))) {
    stress <- function (x, delta, w) {
      n <- nrow (delta)
      p <- length (x) / n
      d <- as.matrix(dist(matrix (x, n, p)))
      res <- delta - d
      return (sum (w * res * res) / 2)
    }
    return (hessian (stress, x, delta = delta, w = w))
  }

cleanUp <- function (a, eps = 1e-3) {
  nv <- nrow (a)
  ind <- rep (TRUE, nv)
  for (i in 1:(nv - 1)) {
    xx <- a[i, ]
    for (j in (i + 1):nv) {
      if (!ind[j])
        next
      yy <- a[j, ]
      mm <- max (abs (xx - yy))
      if (mm < eps)
        ind[j] <- FALSE
    }
  }
  return (ind)
}

bruteForce <- function (a, b, eps = 1e-3) {
  n <- nrow (a)
  m <- ncol (a)
  cb <- combn (n, m)
  n1 <- ncol (cb)
  ind <- rep(TRUE, n1)
  ht <- numeric()
  for (i in 1:n1) {
    gg <- a[cb[, i],]
    bg <- b[cb[, i]]
    qg <- qr(gg)
    if (qg$rank < m) {
      ind[i] <- FALSE
      next
    }
    hh <- solve (qg, bg)
    hg <- drop (a %*% hh)
    if (min (b - hg) < -eps) {
      ind[i] <- FALSE
      next
    }
    ht <- c(ht, hh)
  }
  n2 <- sum (ind)
  ht <- matrix (ht, m, n2)
  ind <-
    .C (
      "cleanup",
      as.double(ht),
      as.integer(n2),
      as.integer(m),
      as.integer(rep(1, n2)),
      as.double (eps)
    )[[4]]
  n3 <- sum (ind)
  return (list (
    x = t(ht)[which(ind == 1),],
    n1 = n1,
    n2 = n2,
    n3 = n3
  ))
}

bruteForceOne <- function (a, b, p, q, v, eps = 1e-3) {
  n <- nrow (a)
  m <- ncol (a)
  ind <- which ((q - v %*% p) > -eps)
  v <- v[ind,]
  cb <- combn (n, m - 1)
  n1 <- ncol (cb)
  ind <- rep(TRUE, n1)
  ht <- numeric()
  for (i in 1:n1) {
    gg <- rbind (a[cb[, i],], p)
    bg <- c (b[cb[, i]], q)
    qg <- qr(gg)
    if (qg$rank < m) {
      ind[i] <- FALSE
      next
    }
    hh <- solve (qg, bg)
    hg <- drop (a %*% hh)
    if (min (b - hg) < -eps) {
      ind[i] <- FALSE
      next
    }
    ht <- c(ht, hh)
  }
  n2 <- sum (ind)
  ht <- t (matrix (ht, m, n2))
  ht <- rbind (v, ht)
  ind <- cleanUp (ht, eps)
  print (ind)
  n3 <- sum (ind)
  return (list (
    x = ht[ind,],
    n1 = n1,
    n2 = n2,
    n3 = n3
  ))
}

rankTest <- function (x, a, b, eps = 1e-3) {
  h <- drop (a %*% x)
  ind <- which (abs (h - b) < eps)
  r <- qr (a[ind, ])$rank
  f <- min (b - h) > -eps
  return (list (rank = r, feasibility = f))
}

makeDC <- function (x) {
  y <- -x
  diag(y) <- -rowSums (y)
  return (y)
}

bmat <- function (delta, w, d) {
  n <- nrow (w)
  dd <- ifelse (d == 0, 0, 1 / d)
  return (makeDC (w * delta * dd))
}

smacof <-
  function (delta,
            w,
            xini,
            eps = 1e-6,
            itmax = 100,
            verbose = TRUE) {
    n <- nrow (xini)
    xold <- xini
    dold <- as.matrix (dist (xold))
    sold <- sum (w * (delta - dold) ^ 2) / 2
    itel <- 1
    v <- ginv (makeDC (w))
    repeat {
      b <- bmat (delta, w, dold)
      xnew <- v %*% b %*% xold
      dnew <- as.matrix (dist (xnew))
      snew <- sum (w * (delta - dnew) ^ 2) / 2
      if (verbose) {
        cat (
          formatC (itel, width = 4, format = "d"),
          formatC (
            sold,
            digits = 10,
            width = 13,
            format = "f"
          ),
          formatC (
            snew,
            digits = 10,
            width = 13,
            format = "f"
          ),
          "\n"
        )
      }
      if ((itel == itmax) || (sold - snew) < eps)
        break ()
      itel <- itel + 1
      sold <- snew
      dold <- dnew
      xold <- xnew
    }
    return (list (
      x = xnew,
      d = dnew,
      s = snew,
      itel = itel
    ))
  }

oneMore <- function (g, u) {
  v <- bruteForce (g, u)$x
  nv <- nrow (v)
  s <- matrix (0, 2, 2)
  ev <- rep (0, nv)
  for (i in 1:nv) {
    s[1, 1] <- v[i, 1]
    s[2, 2] <- v[i, 2]
    s[1, 2] <- s[2, 1] <- v[i, 3]
    ee <- eigen (s)
    ev[i] <- min (ee$values)
    if (ev[i] < 0) {
      yy <- ee$vectors[, 2]
      hh <- c (yy[1] ^ 2, yy[2] ^ 2, 2 * yy[1] * yy [2])
      g <- rbind (g,-hh)
      u <- c (u, 0)
    }
  }
  return (list (
    v = v,
    g = g,
    u = u,
    e = ev
  ))
}

makeG <- function (x) {
  n <- nrow (x)
  p <- ncol (x)
  m <- n - p - 1
  k <- qr.Q(qr(cbind(1, x, diag (n))))[,-c(1:(p + 1))]
  g <- matrix (0, n * (n - 1) / 2, m * (m + 1) / 2)
  l <- 1
  if (m == 1) {
    g[, 1] <- lower_triangle (outer (k, k))
  }
  else {
    for (i in 1:m) {
      g[, l] <- lower_triangle (outer(k[, i], k[, i]))
      l <- l + 1
    }
    for (i in 1:(m - 1))
      for (j in (i + 1):m) {
        g[, l] <-
          lower_triangle (outer(k[, i], k[, j]) + outer(k[, j], k[, i]))
        l <- l + 1
      }
  }
  return (g)
}

iStress <-
  function (x,
            delta,
            w = rep (1, length (delta)),
            only = TRUE) {
    m <- length (delta)
    n <- (1 + sqrt (1 + 8 * m)) / 2
    x <- matrix (x, n, length (x) / n)
    d <- lower_triangle (as.matrix (dist (x)))
    g <- makeG (x)
    h <- (d / w) * makeG (x)
    u <- -colSums(w * (delta - d) * h)
    v <- crossprod (h, w * h)
    s <- solve.QP (
      Dmat = v,
      dvec = u,
      Amat = -t(h),
      bvec = -d
    )
    ds <- d  - h %*% s$solution
    is <- sum (w * (delta - ds) ^ 2)
    if (only)
      return (is)
    else
      return (list (istress = is, delta = fill_symmetric (ds)))
  }
