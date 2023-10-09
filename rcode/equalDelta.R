source("rcode/common/smacof.R")
source("rcode/common/indexing.R")
source("rcode/common/io.R")
       
set.seed(12345)

equalDelta <- function(n, p, m) {
  w <- delta <- wdef(n)
  w <- w / sum(w)
  delta <- delta / sqrt(sum(w * delta ^ 2))
  z <- matrix(0, m, 5)
  y <- list()
  for (i in 1:m) {
    x <- matrix (rnorm(n * p), n, p)
    x <- apply (x, 2, function(x)
      x - mean(x))
    d <- as.matrix(dist(x))
    s <- sum (w * d * delta) / sum (w * d * d)
    x <- x * s
    d <- d * s
    h <-
      smacofR(
        w,
        delta,
        p,
        xold = x,
        verbose = FALSE,
        eps = 1e-15,
        xstop = TRUE,
        itmax = 10000
      )
    z[i,] <- c(i, h$s, h$itel, min (eigen(h$h)$values), max(abs(h$g)))
    y[[i]] <- h$x
  }
  return(list(z, y))
}