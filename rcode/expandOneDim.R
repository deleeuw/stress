expandRho <- function (delta, w = wdef(nrow(delta)), x, y) {
  n <- nrow(delta)
  s0 <- s1 <- s2 <- s3 <- 0
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j)
        next
      del <- delta[i, j]
      www <- w[i, j]
      dxx <- sum((x[i, ] - x[j, ]) ^ 2)
      dxr <- ifelse (dxx < 1e-15, 1, dxx)
      dsx <- sqrt(dxr)
      dyy <- sum((y[i, ] - y[j, ]) ^ 2)
      dxy <- sum((x[i, ] - x[j, ]) * (y[i, ] - y[j, ]))
      vxy <- dyy - ((dxy) ^ 2) / dxr
      s0 <- s0 + www * del * dsx
      s1 <- s1 + www * (del / dsx) * dxy
      s2 <- s2 + www * (del / dsx) * vxy
      s3 <- s3 + www * dxy * (del / (dsx ^ 3)) * vxy
    }
  }
  return(c(s0, s1, s2 / 2, -s3 / 2))
}

expandEta2 <- function (delta, w = wdef(nrow(delta)), x, y) {
  n <- nrow(delta)
  s0 <- s1 <- s2 <- s3 <- 0
  for (i in 1:n) {
    for (j in 1:n) {
      dxx <- sum((x[i, ] - x[j, ]) ^ 2)
      dsx <- sqrt(dxx)
      dyy <- sum((y[i, ] - y[j, ]) ^ 2)
      dxy <- sum((x[i, ] - x[j, ]) * (y[i, ] - y[j, ]))
      s0 <- s0 + w[i, j] * dxx
      s1 <- s1 + w[i, j] * dxy
      s2 <- s2 + w[i, j] * dyy
    }
  }
  return(c(s0, 2 * s1, s2, s3))
}

expandStress <- function (delta, w = wdef(nrow(delta)), x, y) {
  return (c(1, 0, 0, 0) + expandEta2(delta, w, x, y) - 2 * expandRho(delta, w, x, y))
}

expandTester <- function (delta,
                    w = wdef(nrow(delta)),
                    x,
                    y,
                    left = -1,
                    right = 1,
                    length = 1001,
                    order = 3) {
  w <- w / sum (w)
  delta <- delta / sqrt(sum(w * delta ^ 2))
  x <- apply(x, 2, function (x)
    x - mean(x))
  y <- apply(y, 2, function (x)
    x - mean(x))
  d <- as.matrix(dist(x))
  s <- sum(w * d * delta) / sum(w * d * d)
  d <- s * d
  x <- s * x
  h <- expandStress (delta, w, x, y)
  SEQ <- seq(left, right, length = length)
  sig <- rep(0, length)
  sag <- rep(h[1], length)
  for (i in 1:length) {
    z <- x + SEQ[i] * y
    d <- as.matrix (dist(z))
    eta2 <- sum(w * d ^ 2)
    rho <- sum(w * delta * d)
    sig[i] <- 1 - 2 * rho + eta2
    if (order > 0) {
      sag[i] <- sag[i] + SEQ[i] * h[2]
    }
    if (order > 1) {
      sag[i] <- sag[i] + (SEQ[i] ^ 2) * h[3]
    }
    if (order > 2) {
      sag[i] <- sag[i] + (SEQ[i] ^ 3) * h[4]
    }
  }
  return(cbind(sig, sag))
}

