library(combinat)

isRotation <- function (x) {
  n <- nrow (x)
  x <- apply (x, 2, function (x) x - mean (x))
  s <- sum (x ^ 2)
  p <- permn (1:n)
  m <- length (p)
  u <- matrix(0, m, m)
  for (i in 1:m) {
    for (j in 1:m) {
    u[i,j] <- s - sum(svd(crossprod(x[p[[i]],], x[p[[j]],]))$d)
    } 
  }
  return(u)
}

xs <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), 4, 2)

u<-isRotation(xs)
for(i in 1:24) {
  print(which(u[i,]==0))
}
