set.seed(12345)
delta <- 1 - diag (4)
delta <- delta * sqrt (2 / sum (delta ^ 2))
a<-3+sqrt(3)
x<-matrix (c(0, 2*a, a, a, 0, 0, a*sqrt(3), a*sqrt(3)/3), 4, 2)
x <- apply (x, 2, function (x) x - mean (x))
dx<- as.matrix(dist(x))
lb <- sum (delta * dx) / sum (dx ^ 2)
x <- lb * x
z <- matrix (c(0, -2, sqrt(3), -2-sqrt(3)/3, 0, 0, sqrt(3),2 +sqrt(3)), 4, 2)
z <- apply (z, 2, function (x) x - mean (x))

bmat <- function (x) {
  dx <- as.matrix (dist (x))
  b <- - delta / (dx + diag (4))
  diag (b) <- - rowSums (b)
  return (b)
}

stress <- function (x) {
  dx <- as.matrix (dist (x))
  return (sum ((delta - dx) ^ 2) / 2)
}

smacof <- function (xold, itmax = 10000, eps = 1e-15) {
  fold <- stress (xold)
  itel <- 1
  path <- numeric(0)
  repeat {
    b <- bmat (xold)
    xnew <- b %*% xold / 4
    fnew <- stress (xnew)
    path <- c(path, fnew)
    print (noquote (formatC (fnew, di = 15, wi = 18, fo = "f")))
    if (((fold - fnew) < eps) || (itel == itmax)) break
    fold <- fnew
    xold <- xnew
    itel <- itel + 1
  }
  return (list (x = xnew, itel = itel, stress = fnew, path = path))
}