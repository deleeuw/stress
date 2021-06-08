source ("splineBasis.R")

x <- seq(0, 1, by = .05)
innerknots <- c(.3,.5,.6)

h00 <- bsplineBasis(x, 0, innerknots, lowknot  = 0, highknot = 1, type = 0)
h01 <- bsplineBasis(x, 0, innerknots, lowknot  = 0, highknot = 1, type = 1)
h02 <- bsplineBasis(x, 0, innerknots, lowknot  = 0, highknot = 1, type = 2)

h10 <- bsplineBasis(x, 1, innerknots, lowknot  = 0, highknot = 1, type = 0)
h11 <- bsplineBasis(x, 1, innerknots, lowknot  = 0, highknot = 1, type = 1)
h12 <- bsplineBasis(x, 1, innerknots, lowknot  = 0, highknot = 1, type = 2)


h20 <- bsplineBasis(x, 2, innerknots, lowknot  = 0, highknot = 1, type = 0)
h21 <- bsplineBasis(x, 2, innerknots, lowknot  = 0, highknot = 1, type = 1)
h22 <- bsplineBasis(x, 2, innerknots, lowknot  = 0, highknot = 1, type = 2)
