

csupper <- function (lbd, mbd) {
  n <- length (lbd)
  m <- length (mbd)
  mad <- nad <- matrix (0, m, n)
  sad <- rep(0, n)
  plot(lbd,
       lbd,
       type = "n",
       ylab = "",
       ylim = c(0, 2))
  for (k in 1:m) {
    ly <- mbd[k]
    yy <- ly * z1 + (1 - ly) * z2
    dy <- dist (yy)
    by <- as.matrix(-delta / dy)
    diag (by) <- -rowSums(by)
    for (l in 1:n) {
      xx <- lbd[l] * z1 + (1 - lbd[l]) * z2
      dx <- dist (xx)
      rx <- sum (xx * (by %*% yy))
      nx <- sum (dx ^ 2)
      mad[k, l] <- 1 - 2 * rx + nx
    }
    lines (lbd, mad [k,])
    abline (v = ly)
  }
  lines(lbd,
        apply(mad, 2, min),
        type = "l",
        col = "BLUE",
        lwd = 3)
  for (l in 1:n) {
    xx <- lbd[l] * z1 + (1 - lbd[l]) * z2
    dx <- dist (xx)
    sad[l] <- 1 - 2 * sum (dx * delta) + sum (dx ^ 2)
  }
  lines (lbd,
         sad,
         type = "l",
         col = "RED",
         lwd = 3)
}

aglower <- function (lbd, mbd) {
  n <- length (lbd)
  m <- length (mbd)
  mad <- matrix (0, m, n)
  sad <- rep(0, n)
  plot(lbd,
       lbd,
       type = "n",
       ylab = "",
       ylim = c(0, 2))
  for (k in 1:m) {
    ly <- mbd[k]
    yy <- ly * z1 + (1 - ly) * z2
    dy <- dist (yy)
    ry <- sum (delta * dy)
    by <- as.matrix(-delta / dy)
    diag (by) <- -rowSums(by)
    for (l in 1:n) {
      xx <- lbd[l] * z1 + (1 - lbd[l]) * z2
      dx <- dist (xx)
      sx <- sum (xx * (by %*% xx))
      nx <- sum (dx ^ 2)
      mad[k, l] <- 1 - ry + nx - sx
    }
    lines (lbd, mad [k,])
    abline (v = ly)
  }
  lines(lbd,
        apply(mad, 2, max),
        type = "l",
        col = "BLUE",
        lwd = 3)
  for (l in 1:n) {
    xx <- lbd[l] * z1 + (1 - lbd[l]) * z2
    dx <- dist (xx)
    sad[l] <- 1 - 2 * sum (dx * delta) + sum (dx ^ 2)
  }
  lines (lbd,
         sad,
         type = "l",
         col = "RED",
         lwd = 3)
}
