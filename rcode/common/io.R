matrixPrint <- function (x,
                    digits = 6,
                    width = 8,
                    format = "f",
                    flag = "+") {
  print (noquote (formatC (
    x,
    digits = digits,
    width = width,
    format = format,
    flag = flag
  )))
}

iterationWrite <- function (labels, values, digits, width, format) {
  for (i in 1:length(labels)) {
    cat (labels[i],
         formatC(
           values[i],
           di = digits[i],
           wi = width[i],
           fo = format[i]
         ),
         " ")
  }
  cat("\n")
}

rotateEllipse <- function (x) {
  z <- (x[1,] + x[2,]) / 2
  x <- x - matrix (z, nrow(x), 2, byrow = TRUE)
  s <- sqrt (sum (x[1,] ^ 2))
  r <- matrix(c(x[1, 1], x[1, 2],-x[1, 2], x[1, 1]), 2, 2) / s
  x <- x %*% r
  e <- as.matrix (dist (x))
  d <- mean (rowSums(e[-(1:2), 1:2]))
  a <- d / 2
  c <- abs (x[1, 1])
  b <- sqrt (a ^ 2 - c ^ 2)
  return (list(
    x = x,
    a = a,
    b = b,
    c = c
  ))
}

plotEllipse <- function (x) {
  r <- rotateEllipse (x)
  f <- seq (0, 2 * pi, length = 100)
  z <- cbind(sin (f), cos (f))
  z[, 1] <- z[, 1] * r$a
  z[, 2] <- z[, 2] * r$b
  plot(z,
       type = "l",
       col = "RED",
       lwd = 2)
  text (r$x, as.character (1:nrow(r$x)))
  abline(h = 0)
  abline(v = 0)
}

draw_ellipse <- function (center,
                          radius,
                          a = diag (2),
                          np = 100,
                          ...) {
  par (pty = "s")
  e <- eigen(a)
  k <- e$vectors
  lbd <- e$values
  seq <- seq(0, 2 * pi, length = np)
  scos <- (radius * sin (seq)) / lbd[1]
  ccos <- (radius * cos (seq)) / lbd[2]
  sico <- k %*% rbind(scos, ccos) + center
  plot (sico[1,], sico[2,], type = "l", ...)
}






