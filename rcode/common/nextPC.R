
nextPermutation <- function (x) {
  if (all (x == (length(x):1)))
    return (NULL)
  z <- .C("nextPermutation", as.integer(x), as.integer(length(x)))
  return (z[[1]])
}

nextCombination <- function (x, n) {
  m <- length (x)
  if (all (x == ((n - m) + 1:m)))
    return (NULL)
  z <-
    .C("nextCombination",
       as.integer(n),
       as.integer (m),
       as.integer(x))
  return (z[[3]])
}
