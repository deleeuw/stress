delta <- as.matrix (dist (diag (4)))
delta <- delta * sqrt (2 / sum (delta ^ 2))
x <- matrix (c(0,0,1,0,0,1,1,1), 4, 2, byrow = TRUE)
y <- matrix (c(0,0,1,0,.5,.5*sqrt(3)), 3, 2, byrow = TRUE)
y <- rbind(y, colSums(y) / 3)
x <- apply (x, 2, function (z) z - mean (z))
y <- apply (y, 2, function (z) z - mean (z))
dx <- as.matrix (dist (x))
dy <- as.matrix (dist (y))
sx <- sum (delta * dx) / sum (dx ^ 2)
sy <- sum (delta * dy) / sum (dy ^ 2)
x <- sx * x
y <- sy * y
dx <- sx * dx
dy <- sy * dy
bx <- - delta / (dx + diag (4))
diag (bx) <- - rowSums (bx)
by <- - delta / (dy + diag (4))
diag (by) <- - rowSums (by)
mprint (cbind (x, bx %*% x / 4), 10, 15)
mprint (cbind (y, by %*% y / 4), 10, 15)
sx <- sum ((delta - dx) ^ 2) / 2
sy <- sum ((delta - dy) ^ 2) / 2
####
asum <- 2 * 4 * matrix (c (sum(x ^ 2), sum (x * y), sum (x * y), sum (y ^ 2)), 2, 2)
bsum <- chol (asum)
da <- matrix (0, 4, 4)
z1 <- c(1, 0)
for (i in 1:4) for (j in 1:4) da[i,j] <- sqrt (sum (vv (i, j, x, y) * outer(z1, z1)))
db <- matrix (0, 4, 4)
z2 <- c(0, 1)
for (i in 1:4) for (j in 1:4) db[i,j] <- sqrt (sum (vv (i, j, x, y) * outer(z2, z2)))
du <- matrix (0, 4, 4)
z1 <- bsum[,1]
for (i in 1:4) for (j in 1:4) du[i,j] <- sqrt (sum (vnorm(vv (i, j, x, y), x, y) * outer(z1, z1)))
dv <- matrix (0, 4, 4)
z2 <- bsum[,2]
for (i in 1:4) for (j in 1:4) dv[i,j] <- sqrt (sum (vnorm(vv (i, j, x, y), x, y) * outer(z2, z2)))
b1 <- bmat (z1[1], z1[2], x, y, delta)$b
r1 <- sum (z1 * b1 %*% z1)
s1 <- stress (z1[1], z1[2], x, y, delta)
b2 <- bmat (z2[1], z2[2], x, y, delta)$b
s2 <- stress (z2[1], z2[2], x, y, delta)
mprint (cbind (z1, b1 %*% z1), 10, 15)
mprint (cbind (z2, b2 %*% z2), 10, 15)


for (i in 1:4) for (j in 1:4) {
  z <- matrix (0, 2, 2)
  z[, 1] <- x[i, ] - x[j, ]
  z[, 2] <- y[i, ] - y[j, ]
  e <- eigen (crossprod (z))$values
  print (c(i, j, min (e)))
}

i<-1
j<-2
z <- matrix (0, 2, 2)
z[, 1] <- x[1, ] - x[2, ]
z[, 2] <- y[1, ] - y[2, ]
e0<-eigen(crossprod(z))