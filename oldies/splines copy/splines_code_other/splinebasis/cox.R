dyn.load("cox.so")

knots <- c(0, 0, 0, .3, .5, .6, 1, 1, 1)

coxR <- function (i, x, order, knots) {
  h <-
    .C(
      "cox",
      as.integer(length(knots)),
      as.integer(i),
      as.integer(order),
      as.double(x),
      as.double(knots)
    )
  return (h)
}
