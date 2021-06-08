f <- function (n, m, blks, ties) {
  for (j in 1:m) {
    x <- sample (1:blks, n, replace = TRUE)
    y <- rnorm (n)
    h <- monregR (x, y, ties = ties)
  }
}

g <- function (n, m, blks, ties) {
  for (j in 1:m) {
    x <- sample (1:blks, n, replace = TRUE)
    y <- rnorm (n)
    h <- monregRC (x, y, ties = ties)
  }
}

n <- 10000
m <- 100
ties <- c(1, 2, 3)
blks <- c(2, 10, 100, 1000, 10000)
set.seed(12345)
for (i in blks) {
  for (j in ties) {
    cat(
      "blks =",
      formatC (i, digits = 6, format = "d"),
      "ties = ",
      formatC (j, digits = 6, format = "d"),
      "\n"
    )
    print(system.time(f(n, m, i, j)))
    print(system.time(g(n, m, i, j)))
  }
}
