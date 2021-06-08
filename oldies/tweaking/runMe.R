dyn.load("smacof.so")

source("smacofR.R")
source("smacofRC.R")
source("utilsRC.R")

m <- 100
n <- 50
wv <- dv <- rep (1, n * (n - 1) / 2)
wm <- dm <- 1 - diag (n)
set.seed(12345)
user1 <- rep(0, m)
user2 <- rep(0, m)
user3 <- rep(0, m)
for (j in 1:m) {
  xold <- rnorm(2*n)
  t1 <- system.time(h1 <- smacofR(wm, dm, 2, xold = matrix(xold, n, 2), eps = 1e-6,  itmax = 1000), gcFirst = TRUE)
  t2 <- system.time(h2 <- smacofRC(wv, dv, 2, xold = xold, eps = 1e-6, itmax = 1000), gcFirst = TRUE)
  t3 <- system.time(h3 <- smacofRCU(wv, dv, 2, xold = xold, eps = 1e-6, itmax = 1000), gcFirst = TRUE)
  user1[j] <- t1[1]
  user2[j] <- t2[1]
  user3[j] <- t3[1]
}