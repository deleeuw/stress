baseplot <- function (x,
                      y,
                      z,
                      wx = TRUE,
                      wy = TRUE,
                      wz = TRUE) {
  par(pty="s")
  plot(
    x,
    xlim = c(-3, 3),
    ylim = c(-3, 3),
    xlab = "",
    ylab = "",
    type  = "n"
  )
  if (wx)
    points(x, col = "RED", cex = 1)
  if (wy)
    points(y, col = "BLUE", cex = 1)
  if (wz)
    points(z, col = "GREEN", cex = 1)
  mx <- apply(x, 2, mean)
  my <- apply(y, 2, mean)
  mz <- apply(z, 2, mean)
  if (wx)
    points(
      matrix(mx, 1, 2),
      col = "RED",
      pch = 5,
      cex = 2,
      lwd = 2
    )
  if (wy)
    points(
      matrix(my, 1, 2),
      col = "BLUE",
      pch = 5,
      cex = 2,
      lwd = 2
    )
  if (wz)
    points(
      matrix(mz, 1, 2),
      col = "GREEN",
      pch = 5,
      cex = 2,
      lwd = 2
    )
  if (wx)
    for (i in 1:10) {
      lines(rbind(x[i,], mx))
    }
  if (wy)
    for (i in 1:5) {
      lines(rbind(y[i,], my))
    }
  if (wz)
    for (i in 1:5) {
      lines(rbind(z[i,], mz))
    }
}
