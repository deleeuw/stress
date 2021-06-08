dyn.load("jacobi.so")

jacobi <- function (a,
                    itmax = 10,
                    eps = 1e-6,
                    verbose = FALSE) {
  m <- length (a)
  n <- round ((sqrt (1 + 8 * m) - 1) / 2)
  r <- -((1:n) ^ 2) / 2 + (2 * n + 3) * (1:n) / 2 - n
  h <-
    .C(
      "jacobiC",
      as.integer(n),
      eval = as.double (a),
      evec = as.double (rep(0, n * n)),
      as.double (rep(0, n)),
      as.double (rep(0, n)),
      as.integer(itmax),
      as.double(eps)
    )
  return (list (eval = h$eval[r],
                evec = matrix(h$evec, n, n)))
}


triangle2matrix <- function (x) {
  m <- length (x)
  n <- round ((sqrt (1 + 8 * m) - 1) / 2)
  h <-
    .C("trimat", as.integer (n), as.double (x), as.double (rep (0, n * n)))
  return (matrix(h[[3]], n, n))
}

matrix2triangle <- function (x) {
  n <- dim(x)[1]
  m <- n * (n + 1) / 2
  h <-
    .C("mattri", as.integer (n), as.double (x), as.double (rep (0, m)))
  return (h[[3]])
}

trianglePrint <- function (x, width = 6, precision = 4) {
  n <- round ((sqrt (1 + 8 * m) - 1) / 2)
  h <- .C("pritru", as.integer(n), as.integer(width), as.integer(precision), as.double (x))
}

matrixPrint <- function (x, width = 6, precision = 4) {
  n <- nrow (x)
  m <- ncol (x)
  h <- .C("primat", as.integer(n), as.integer(m), as.integer(width), as.integer(precision), as.double (x))
}

jacobiR <-
  function(a,
           eps1 = 1e-10,
           eps2 = 1e-6,
           itmax = 100,
           vectors = TRUE,
           verbose = FALSE) {
    n <- nrow(a)
    k <- diag(n)
    itel <- 1
    mx <- 0
    saa <- sum(a ^ 2)
    repeat {
      for (i in 1:(n - 1))
        for (j in (i + 1):n) {
          aij <- a[i, j]
          bij <- abs(aij)
          if (bij < eps1)
            next()
          mx <- max(mx, bij)
          am <- (a[i, i] - a[j, j]) / 2
          u <- c(aij, -am)
          u <- u / sqrt(sum(u ^ 2))
          c <- sqrt((1 + u[2]) / 2)
          s <- sign(u[1]) * sqrt((1 - u[2]) / 2)
          ss <- s ^ 2
          cc <- c ^ 2
          sc <- s * c
          ai <- a[i, ]
          aj <- a[j, ]
          aii <- a[i, i]
          ajj <- a[j, j]
          a[i, ] <- a[, i] <- c * ai - s * aj
          a[j, ] <- a[, j] <- s * ai + c * aj
          a[i, j] <- a[j, i] <- 0
          a[i, i] <- aii * cc + ajj * ss - 2 * sc * aij
          a[j, j] <- ajj * cc + aii * ss + 2 * sc * aij
          if (vectors) {
            ki <- k[, i]
            kj <- k[, j]
            k[, i] <- c * ki - s * kj
            k[, j] <- s * ki + c * kj
          }
        }
      ff <- sqrt(saa - sum(diag(a) ^ 2))
      if (verbose)
        cat(
          "Iteration ",
          formatC(itel, digits = 4),
          "maxel ",
          formatC(mx, width = 10),
          "loss ",
          formatC(ff, width = 10),
          "\n"
        )
      if ((mx < eps1) || (ff < eps2) || (itel == itmax))
        break()
      itel <- itel + 1
      mx <- 0
    }
    d <- diag(a)
    o <- order(d, decreasing = TRUE)
    if (vectors)
      return(list(values = d[o], vectors = k[, o]))
    else
      return(values = d[o])
  }
