dorpolRC <- function (n) {
  h <-
    .C("dorpol", as.integer(n), as.double (rep(0, n * n)))
  return (matrix(h[[2]], n, n))
}

docentRC <- function (x) {
  m <- length (x)
  n <- round((sqrt (1 + 8 * m) - 1) / 2)
  h <-
    .C("docent", as.integer (n), as.double (x), as.double (rep(0, m)))
  return(h$y)
}

primatRC <- function (n,
                      m,
                      x,
                      width = 6,
                      precision = 4) {
  h <-
    .C(
      "primat",
      as.integer (n),
      as.integer (m),
      as.integer (width),
      as.integer (precision),
      as.double (x)
    )
}

pritrlRC <- function (x,
                      width = 6,
                      precision = 4) {
  m <- length (x)
  n <- round ((sqrt (1 + 8 * m) + 1) / 2)
  h <-
    .C("pritrl",
       as.integer (n),
       as.integer (width),
       as.integer (precision),
       as.double (x))
}


pritruRC <- function (x,
                      width = 6,
                      precision = 4) {
  m <- length (x)
  n <- round ((sqrt (1 + 8 * m) - 1) / 2)
  h <-
    .C("pritru",
       as.integer (n),
       as.integer (width),
       as.integer (precision),
       as.double (x))
}


priarrRC <- function (n,
                      m,
                      r,
                      x,
                      width = 6,
                      precision = 4) {
  h <-
    .C(
      "primat",
      as.integer (n),
      as.integer (m),
      as.integer (r),
      as.integer(width),
      as.integer (precision),
      as.double (x)
    )
}

mpowerRC <- function (x, p) {
  m <- length (x)
  n <- round ((sqrt (1 + 8 * m) - 1) / 2)
  h <-
    .C("mpower",
       as.integer(n),
       as.double (x),
       as.double(p),
       xpow = as.double(rep(0, m)))
  return (h$xpow)
}

trimatRC <- function (x) {
  m <- length (x)
  n <- round ((sqrt (1 + 8 * m) - 1) / 2)
  h <-
    .C("trimat", as.integer(n), as.double (x), y = as.double (rep(0, n * n)))
  return(h$y)
}

mattriRC <- function (x) {
  n <- round (sqrt (length (x)))
  m <- n * (n + 1) / 2
  h <-
    .C("mattri", as.integer(n), as.double (x), y = as.double (rep(0, m)))
  return(h$y)
}

mutrmaRC <- function (n, m, a, x) {
  h <-
    .C("mutrma",
       as.integer(n),
       as.integer(m),
       as.double (a),
       as.double (x),
       y = as.double (rep(0, n * m)))
  return (h$y)
}
