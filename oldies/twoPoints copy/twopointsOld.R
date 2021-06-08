delta <- as.matrix (dist (diag (4)))
delta <- delta / sqrt (sum (delta ^2))
x1 <- matrix (c(0,0,1,0,0,1,1,1), 4, 2, byrow = TRUE)
x2 <- matrix (c(0,0,1,0,.5,.5*sqrt(3)), 3, 2, byrow = TRUE)
x2 <- rbind(x2, colSums(x2) / 3)
x1 <- apply (x1, 2, function (z) z - mean (z))
x2 <- apply (x2, 2, function (z) z - mean (z))
d1 <- as.matrix (dist (x1))
d2 <- as.matrix (dist (x2))
lbd1 <- sum (delta * d1) / sum (d1 ^ 2)
x1 <- x1 * lbd1
d1 <- d1 * lbd1
lbd2 <- sum (delta * d2) / sum(d2 ^ 2)
x2 <- x2 * lbd2
d2 <- d2 * lbd2
bmat <- function (delta, d) {
  b <- - as.matrix (delta) / (as.matrix (d) + diag (4))
  diag (b) <- - rowSums (b)
  return (b)
}
b1 <- bmat (delta, d1) / 4
b2 <- bmat (delta, d2) / 4
y1 <- cbind (x1, b1 %*% x1)
y2 <- cbind (x2, b2 %*% x2)
vv <- function (i, j) {
  a <- matrix (0, 2, 2)
  a[1, 1] <- sum ((x1[i, ]- x1[j,]) ^ 2)
  a[2, 2] <- sum ((x2[i, ]- x2[j,]) ^ 2)
  a[1, 2] <- a[2, 1] <- sum ((x1[i, ]- x1[j,]) * (x2[i, ]- x2[j, ]))
  return (a)
}
stress <- function (x, y) {
  st <- 0
  z <- c(x, y)
  for (i in 1:4)
    for (j in 1:4) {
      dd <- sum (vv (i, j) * outer (z, z))
      st <- st + ((delta[i,j] - sqrt (dd)) ^ 2)
    }
  return (st)
}

pdf("stresser.pdf")

aa <- bb <- seq (-2, 2, length = 100)
z <- matrix (0, 100, 100)
for (i in 1:100) for (j in 1:100)
  z[i, j]<- stress (aa[i], bb[j])
contour(aa, bb, z, levels = seq(0,.50,length=50))


smacof <- function (x, y) {
  xold <- x * x1 + y * x2
  dold <- as.matrix (dist (xold))
  print (sum ((delta - dold) ^ 2))
  repeat {
    b <- bmat (delta, dold)
    xnew <- b %*% xold / 4
    dnew <- as.matrix (dist (xnew))
    print (sum ((delta - dnew) ^ 2))
    if (max (abs (dold - dnew)) < 1e-6) break
    xold <- xnew
    dold <- dnew
  }
}

aa <- seq (-1.2, -.8, length = 100)
bb <- seq (1.4, 1.6, length = 100)
z <- matrix (0, 100, 100)
for (i in 1:100) for (j in 1:100)
  z[i, j]<- stress (aa[i], bb[j])
contour(aa, bb, z, levels = seq(.1106,.1130,length=50))

h<-optim(c(-1.2,1.5),fn = function(x) stress (x[1],x[2]), lower = c(-1.2,1.4), upper = c(-0.8,1.6), method="L-BFGS-B", control = list (trace = 6))
xx<-h$par[1]*x1+h$par[2]*x2
dd <- as.matrix (dist(xx))
bb <- bmat (delta, dd)
yb <- cbind (xx, bb%*%xx / 4)

aa <- seq (-0.1, .1, length = 100)
bb <- seq (0.9, 1.1, length = 100)
z <- matrix (0, 100, 100)
for (i in 1:100) for (j in 1:100)
  z[i, j]<- stress (aa[i], bb[j])
contour(aa, bb, z, levels = seq(0.066,.08,length=100))

bb <- seq (-0.1, .1, length = 100)
aa <- seq (0.9, 1.1, length = 100)
z <- matrix (0, 100, 100)
for (i in 1:100) for (j in 1:100)
  z[i, j]<- stress (aa[i], bb[j])
contour(aa, bb, z, levels = seq(0.025,.030,length=100))

aa <- seq (-1, -.8, length = 100)
bb <- seq (1.4, 1.5, length = 100)
z <- matrix (0, 100, 100)
for (i in 1:100) for (j in 1:100)
  z[i, j]<- stress (aa[i], bb[j])
contour(aa, bb, z, levels = seq(.112,.113,length=50))

dev.off()
