set.seed(12345)
library(ggm)
dyn.load("hsweep.so")
b <- crossprod(matrix(rnorm(1000000),1000,1000))/1000

rSweep <- function (a, indi) {
for (j in indi) {
    pv <- a[j, j]
    if (pv == 0) next ()
    pr <- a[j, -j]
    pc <- a[-j, j]
    a[j, j] <- -1 / pv
    a[j, -j] <- pr / pv
    a[-j, j] <- pc / pv
    a[-j, -j] <- a[-j, -j] - outer(pc, pr) / pv
    }
return (a)
}

fSweep <- function (b , ind, eps = 1e-10, trianin = FALSE, trianout = FALSE) {
	if (trianin) {
		m <- length (b)
		n <- as.integer ( (sqrt (1 + 8 * m) - 1) / 2)
		a <- b
	} else {
		n <- nrow (b)
		m <- (n * n + n) / 2
		a <- b[upper.tri (b, diag = TRUE)]
	}
	for (k in ind) {
		a <- .Fortran("gsweep",
			diago = as.double (a[cumsum (1:n)]),
			trian = as.double (a),
			k = as.integer (k),
			l = as.integer (0),
			m = as.integer (m),
			n = as.integer (n),
			e = as.double (eps),
			ifault = as.integer (0))$trian
		}
	if (trianout) {
		return (a)
	} else {
		return (tri2sym (a))
	}
}

cSweep <- function (b , ind, eps = 1e-10, trianin = FALSE, trianout = FALSE) {
	if (trianin) {
		m <- length (b)
		n <- as.integer ( (sqrt (1 + 8 * m) - 1) / 2)
		a <- b
	} else {
		n <- nrow (b)
		m <- (n * n + n) / 2
		a <- b[upper.tri (b, diag = TRUE)]
	}
	a <- .C("hsweep",
		diago = as.double (a[cumsum (1:n)]),
		trian = as.double (a),
		n = as.integer (n),
		m = as.integer (m),
		ind = as.integer (ind),
		npiv = as.integer (length (ind)),
		eps = as.double (eps))$trian
	if (trianout) {
		return (a)
	} else {
		return (tri2sym (a))
	}
}

tri2sym <- function (a) {
	m <- length (a)
	n <- as.integer ( (sqrt (1 + 8 * m) - 1) / 2)
	b <- matrix (0, n, n)
	b[upper.tri (b, diag = TRUE)] <- a
	b <- b + t(b)
	diag (b) <- diag (b) / 2
	return (b)
}

sweepTimer <- function (howMany = 5, whichPart = 1:500) {
	print ("element sweeping in pure R")
	gc()
	for (i in 1:howMany) {
		print (system.time(rSweep(b, 1:500)))
	}
	print ("element sweeping in R and fortran from R")
	gc()
	for (i in 1:howMany) {
		print (system.time(fSweep(b, 1:500)))
	}
	print ("element sweeping in R and C from R and fortran from C")
	gc()
	for (i in 1:howMany) {
		print (system.time(cSweep(b, 1:500)))
	}
	print ("matrix sweeping from ggm using solve")
	gc()
	for (i in 1:howMany) {
		print (system.time(swp(b, 1:500)))
	}
	print ("plain old solve")
	gc()
	for (i in 1:howMany) {
		print (system.time(solve(b)))
	}
}
