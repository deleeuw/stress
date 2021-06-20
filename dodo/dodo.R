library(numDeriv)
set.seed(12345)
delta <- dist(matrix(rnorm(18), 9, 2))
z <- seq(-2*pi,2*pi, length = 10)[1:9]
f <- function (z) {
  x <- cbind(sin(z), cos(z))[1:9,]
  d <- dist (x)
  return (sum ((delta - d) ^ 2) / 2)
}
ng <- grad (f, z)
nh <- hessian(f, z)
x <- cbind(sin(z), cos(z))
d <- dist (x)
b <- -as.matrix(delta / d)
v <- 9 * diag(9) - 1
diag(b) <- -rowSums(b)
dxp <- cbind(cos(z),-sin(z))
dsx <- (v - b) %*% x
fg <- rep(0, 9)
for (s in 1:9) {
  gs <- matrix(0, 9, 2)
  gs[s,] <- dxp[s, ]
  fg[s] <- sum (dsx * gi)
}
ddsx <- matrix(0, 18, 18)
ddsx[1:9,1:9] <- v - b
ddsx[10:18,10:18] <- v - b

for (s in 1:9) {
  for (t in 1:9) {
    
  }
}