strainAdd <-
  function (delta,
            w = rep (1, length (delta)),
            p = 2,
            itmax = 100,
            eps = 1e-6,
            verbose = TRUE) {
    delta <- as.matrix (delta ^ 2)
    n <- nrow(delta)
    
  }

strainWeight <-
  function (delta,
            w = rep (1, length (delta)),
            p = 2,
            itmax = 100,
            eps = 1e-6,
            verbose = TRUE) {
    
  }

alscal <-
  function (delta,
            p,
            x = torgerson (delta, p),
            w = wdef (nrow (x)),
            itmax = 1000,
            eps = 1e-6,
            verbose = TRUE,
            check = TRUE) {
    n <- nrow (x)
    delta <- delta ^ 2
    d <- as.matrix (dist (x)) ^ 2
    sold <- sum (w * (delta - d) ^ 2)
    wsum <- rowSums (w)
    itel <- 1
    snew <- sold
    repeat {
      for (k in 1:n) {
        t4 <- wsum[k]
        for (s in 1:p) {
          u <- x[, s] - x[k, s]
          t0 <- snew
          t1 <- t2 <- t3 <- 0
          for (i in 1:n) {
            t1 <- t1 + 4 * w[i, k] * (d[i, k] - delta[i, k]) * u[i]
            t2 <-
              t2 + 2 * w[i, k] * ((d[i, k] - delta[i, k]) + 2 * u[i] ^ 2)
            t3 <- t3 + 4 * w[i, k] * u[i]
          }
          pp <- polynomial(c(t0,-t1, t2,-t3, t4))
          qq <- deriv (pp)
          ss <- solve (qq)
          ss <- Re (ss[which (abs (Im (ss)) < 1e-10)])
          tt <- predict (pp, ss)
          snew <- min (tt)
          root <- ss[which.min (tt)]
          x[k, s] <- x[k, s] + root
          for (i in (1:n)[-k]) {
            d[i, k] <- d[i, k] - 2 * root * u[i] + root ^ 2
            d[k, i] <- d[i, k]
          }
        }
      }
      if (verbose) {
        cat(
          "itel ",
          formatC(itel, width = 6, format = "d"),
          "sold ",
          formatC(
            sold,
            digits = 6,
            width = 15,
            format = "f"
          ),
          "snew ",
          formatC(
            snew,
            digits = 6,
            width = 15,
            format = "f"
          ),
          "\n"
        )
      }
      if (((sold - snew) < eps) || (itel == itmax)) {
        break
      }
      sold <- snew
      itel <- itel + 1
    }
    return (list (
      x = x,
      sstress = snew,
      itel = itel
    ))
  }

jeffrey <- function(a) {
  h <-
    .C("jeffreyC",
       a = as.double (a),
       minwhere = as.double (0),
       minvalue = as.double (0))
  return (list.remove (h, 1))
}

