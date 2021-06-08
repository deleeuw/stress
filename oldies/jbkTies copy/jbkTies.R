dyn.load("jbkTies.so")

mySortDouble <- function (x, y, w = rep(1, length(x))) {
  n <- length (x)
  xind <- rep (0, n)
  h <-
    .C(
      "mySortDouble",
      x = as.double (x),
      y = as.double (y),
      w = as.double (w),
      xind = as.integer (xind),
      n = as.integer(n)
    )
  return (h)
}

mySortInteger <- function (x) {
  n <- length (x)
  xind <- rep (0, n)
  h <-
    .C(
      "mySortInteger",
      x = as.integer (x),
      xind = as.integer (xind),
      n = as.integer(n)
    )
  return (h)
}

tieBlock <- function (x) {
  n <- length (x)
  h <-
    .C(
      "tieBlock",
      x = as.double (x),
      iblks = as.integer (rep(0, n)),
      as.integer (n),
      nblk = as.integer (0)
    )
  return (h)
}

makeBlocks <- function  (x, w = rep (1, length(x)), iblks) {
  n <- length (x)
  nblk <- max (iblks)
  xblks <- rep (0, nblk)
  wblks <- rep (0, nblk)
  h <-
    .C(
      "makeBlocks",
      x = as.double (x),
      w = as.double (w),
      xblks = as.double (xblks),
      wblks = as.double (wblks),
      iblks = as.integer (iblks),
      n = as.integer(n),
      nblk = as.integer (nblk)
    )
  return (h)
}

sortBlocks <- function (y, w, xind, iblks) {
  n <- length (y)
  nblk <- max (iblks)
  h <- .C(
    "sortBlocks",
    y = as.double (y),
    w = as.double (w),
    xind = as.integer (xind),
    as.integer(iblks),
    as.integer (n),
    as.integer (nblk)
  )
  return (h)
}

jbkPava <- function (x, w = rep(1, length(x))) {
  h <-
    .C("jbkPava",
       x = as.double(x),
       w = as.double (w),
       n = as.integer(length(x)))
  return (h)
}

monregR <-
  function (x,
            y,
            w = rep (1, length (x)),
            ties = 1) {
    f <- sort(unique(x))
    g <- lapply(f, function (z)
      which(x == z))
    n <- length (x)
    k <- length (f)
    if (ties == 1) {
      h <- lapply (g, function (x)
        y[x])
      f <- lapply (g, function (x)
        w[x])
      m <- rep (0, n)
      for (i in 1:k) {
        ii <- order (h[[i]])
        g[[i]] <- g[[i]][ii]
        h[[i]] <- h[[i]][ii]
        f[[i]] <- f[[i]][ii]
      }
      r <- jbkPava (unlist(h), unlist(f))$x
      s <- r[order (unlist (g))]
    }
    if (ties == 2) {
      h <- lapply (g, function (x)
        y[x])
      f <- lapply (g, function (x)
        w[x])
      s <- sapply(1:k, function(i) sum (h[[i]] * f[[i]]))
      r <- sapply (f, sum)
      m <- s / r
      r <- jbkPava (m, r)$x
      s <- rep (0, n)
      for (i in 1:k)
        s[g[[i]]] <- r[i]
    }
    if (ties == 3) {
      h <- lapply (g, function (x)
        y[x])
      f <- lapply (g, function (x)
        w[x])
      s <- sapply(1:k, function(i) sum (h[[i]] * f[[i]]))
      r <- sapply (f, sum)
      m <- s / r
      r <- jbkPava (m, r)$x
      s <- rep (0, n)
      for (i in 1:k)
        s[g[[i]]] <- y[g[[i]]] + (r[i] - m[i])
    }
    return (s)
  }

monregRC <- function (x,
                      y,
                      w = rep(1, length(x)),
                      ties = 1) {
  a <- mySortDouble(x, y, w)
  b <- tieBlock (a$x)$iblks
  if (ties == 1) {
    c <- sortBlocks (a$y, a$w, a$xind, b)
    d <- jbkPava (c$y, c$w)$x
    e <- mySortInteger (c$xind)$xind
    return (d[e])
  }
  if (ties == 2) {
    c <- makeBlocks (a$y, a$w, b)
    d <- jbkPava (c$xblks, c$wblks)$x
    e <- mySortInteger (a$xind)$xind
    return (d[b][e])
  }
  if (ties == 3) {
    c <- makeBlocks (a$y, a$w, b)
    d <- jbkPava (c$xblks, c$wblks)$x
    u <- a$y - (c$xblks - d)[b]
    e <- mySortInteger (a$xind)$xind
    return (u[e])
  }
}
