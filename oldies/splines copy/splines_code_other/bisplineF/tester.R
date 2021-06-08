source ("BISpline.R")

x <- seq(0, 1, by = .05)
innerknots <- c(.3,.5,.6)

h01 <- BISpline(x, 0, innerknots, lowknot  = 0, highknot = 1, type = 1)
h02 <- BISpline(x, 0, innerknots, lowknot  = 0, highknot = 1, type = 2)

h11 <- BISpline(x, 1, innerknots, lowknot  = 0, highknot = 1, type = 1)
h12 <- BISpline(x, 1, innerknots, lowknot  = 0, highknot = 1, type = 2)


h21 <- BISpline(x, 2, innerknots, lowknot  = 0, highknot = 1, type = 1)
h22 <- BISpline(x, 2, innerknots, lowknot  = 0, highknot = 1, type = 2)
