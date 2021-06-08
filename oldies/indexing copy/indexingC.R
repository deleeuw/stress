dyn.load("indexing.so")

fArrayFirstRC <- function (cell, dimension) {
  rank <- length (cell)
  h <-
    .C(
      "fArrayFirstGlue",
      as.integer(cell),
      as.integer(dimension),
      as.integer(rank),
      index = as.integer(0)
    )
  return (h$index)
}


fSupSymIncreasingFirstRC <- function (cell) {
  rank <- length(cell)
  h <-
    .C("fSupSymIncreasingGlue",
       as.integer(cell),
       as.integer(rank),
       index = as.integer(0))
  return (h$index)
}

fArrayFirstInverseRC <- function(index, dimension) {
  rank <- length (dimension)
  h <-
    .C(
      "fArrayFirstInverse",
      as.integer(rank),
      as.integer (dimension),
      as.integer(index),
      cell = as.integer(rep(0, rank))
    )
  return (h$cell)
}

fSupSymIncreasingFirstInverseRC <- function (dimension, rank, index) {
  h <-
    .C(
      "fSupSymIncreasingFirstInverse",
      as.integer(dimension),
      as.integer(rank),
      as.integer(index),
      cell = as.integer(rep(0, rank))
    )
  return (h$cell)
}
