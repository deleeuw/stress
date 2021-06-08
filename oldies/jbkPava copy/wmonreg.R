dyn.load("wmonreg.so")

wmonreg <- function (x, w = rep(1, length(x))) {
h <- .C("wmonreg", as.double(x), as.double (w), as.integer(length(x)))
return (h[[1]])
}
