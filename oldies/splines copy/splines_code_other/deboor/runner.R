source("deboor.R")

set.seed
innerknots <- c(1,2,3)
multiplicities <- c(1,1,1)
order <- 3
degree <- 2
lowend <- 0
highend <- 4

x <- seq (0, 4, length = 1000)
h <- bsplineBasis (x, innerknots, multiplicities, order, lowend, highend)
k <- ncol (h)
pdf ("plotter_111.pdf")
par (mfrow=c(3,3))
#l<-1-t(apply(h,1,cumsum))
for (j in 1:k)
  plot (x, h[, j], type="l", col = "RED", lwd = 3)
dev.off()

ff <- function (n, m) {
	set.seed(12345)
    for (j in 1:m) {
    	x <- runif (n, 0, 4)
        h <- bsplineBasis (x, innerknots, multiplicities, order, lowend, highend)
    }
}


