data(morse, package = "smacof")
morse <- as.matrix (morse)
w <- matrix (1, 36, 36) - diag (36)
morse <- 2 * morse / sqrt (sum (w * morse * morse))
lbd <- (0:10000)/1000
sink("morse.out")
hMorse <-runPenalty (w, morse, p = 1, lbd = lbd, eps = 1e-10, cut = 1e-10)
sink()