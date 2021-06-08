dyn.load("mySort.so")

mySort <- function (x) {
  h <-
    .C(
      "mySort",
      values = as.double (x),
      indices = as.integer(1:length(x)),
      as.integer(length(x))
    )
  return (list (values = h$values, indices = h$indices))
}