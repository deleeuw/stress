gramy <- function (y, v) {
  r <- length (y)
  s <- sum (y[[1]] * (v %*% y[[1]]))
  y[[1]] <- y[[1]] / sqrt (s)
  for (j in 2:r) {
    for (i in 1:(j - 1)) {
      s <- sum (y[[i]] * (v %*% y[[j]]))
      y[[j]] <- y[[j]] - s * y[[i]]
    }
    s <- sum (y[[j]] * v %*% y[[j]])
    y[[j]] <- y[[j]] / sqrt (s)
  }
  return (y)
}

basis <- function (w, n, p) {
  set.seed(12345)
  m <- n * p - p * (p + 1) / 2
  y <- as.list(rep(0, m))
  v <- -as.matrix(w)
  diag(v) <- -rowSums(v)
  for (j in 1:m) {
    y[[j]] <- matrix(rnorm(n * p), n, p)
    for (s in 1:p) {
      if (s == 1) {
        y[[j]][, s] <- y[[j]][, s] - mean(y[[j]][, s])
      } else {
        ss <- 1:(s - 1)
        y[[j]][ss, s] <- 0
        y[[j]][-ss, s] <- y[[j]][-ss, s] - mean(y[[j]][-ss, s])
      }
    }
  }
  return (gramy (y, v))
}


binary <- function (n) {
  x <- matrix(0, 2 ^ n, n)
  for (j in 1:n) {
    x[, n - j + 1] <-
      rep(c(rep(0, 2 ^ (j - 1)), rep(1, 2 ^ (j - 1))), 2 ^ (n - j))
  }
  row.names(x) <- as.character(0:((2 ^ n) - 1))
  return(x)
}

dissim <- function (w) {
  n <- attr(w, "Size")
  delta <- as.dist(matrix(1, n, n))
  return(delta / sqrt(sum(w * delta ^ 2)))
}

rhoval <- function (alpha, w, delta, base) {
  z <- alpha[1] * base[[1]]
  for (i in 2:length(base)) {
    z <- z + alpha[i] * base [[i]]
  }
  return(sum(w * delta * dist(z)))
}

w <- as.dist(matrix(1, 5, 5))
delta <- dissim(w)
base <- basis(w, 5, 2)
alpha <- rep(1 / sqrt(7), 7)
outer <- 2 * binary(7) - 1
inner <- rbind(diag(7), -diag(7))
outval <-
  apply (outer, 1, function(x)
    rhoval(x, w = w, delta = delta, base = base))
inval <-
  apply (inner, 1, function(x)
    rhoval(x, w = w, delta = delta, base = base))

sma2rho <-
  function (alpha,
            w,
            delta,
            base,
            itmax = 1000,
            eps = 1e-6,
            verbose = FALSE) {
    n <- length(alpha)
    itel <- 1
    alphaold <- alpha
    rhoold <- rhoval (alphaold, w, delta, base)
    repeat {
      z <- alphaold[1] * base[[1]]
      for (i in 2:n) {
        z <- z + alphaold[i] * base [[i]]
      }
      b <- -as.matrix (w * delta / dist (z))
      diag(b) <- -rowSums(b)
      c <- matrix (0, n, n)
      for (i in 1:n) {
        for (j in 1:n) {
          c[i, j] <- sum (base[[i]] * (b %*% base[[j]]))
        }
      }
      alphanew <- drop(c %*% alphaold)
      alphanew <- alphanew / sqrt(sum (alphanew ^ 2))
      rhonew <- rhoval (alphanew, w, delta, base)
      if (verbose) {
        cat(
          "itel ",
          formatC (itel, format = "d"),
          " rhoold ",
          formatC(
            rhoold,
            digits = 6,
            width = 8,
            format = "f"
          ),
          " rhonew ",
          formatC(
            rhonew,
            digits = 6,
            width = 8,
            format = "f"
          ),
          "\n"
        )
      }
      if (((rhonew - rhoold) < eps) || (itel == itmax)) {
        break
      }
      itel <- itel + 1
      rhoold <- rhonew
      alphaold <- alphanew
    }
    return (list(
      alpha = alphanew,
      rho = rhonew,
      itel = itel
    ))
  }

addVertices <- function() {
  
}
