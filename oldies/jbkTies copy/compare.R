h <- function (n, m, ties) {
  mm <- 0
  for (j in 1:m) {
    x <- sample (1:10, n, replace = TRUE)
    y <- rnorm (n)
    h <- monregRC (x, y, ties = ties)
    g <- monregR (x, y, ties = ties)
    mm <- max (mm, max (abs (h - g)))
  }
  print (mm)
}
  
  h(n, m, 1)
  h(n, m, 2)
  h(n, m, 3)
  