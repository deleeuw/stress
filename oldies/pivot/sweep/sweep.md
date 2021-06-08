# Partial inverses by Sweeping
Jan de Leeuw  
March 21, 2016  



#Simple Sweep

Suppose \(A\) is a square non-singular matrix. A *simple sweep* on the \(k^{th}\) element replaces \(A\) by the matrix with elements

\begin{align*}
a_{kk}&=\frac{1}{a_{kk}},\\
a_{ki}&=-\frac{a_{ki}}{a_{kk}},\\
a_{ik}&=\frac{a_{ik}}{a_{kk}},\\
a_{ij}&=a_{ij}-\frac{a_{ik}a_{kj}}{a_{kk}}
\end{align*}
for \(j\not= k\) and \(i\not= k\).

In R that becomes

```r
sweep <- function (a, k) {
  p <- a [k, k] 
  a [k, k] <- 1 /p
  a [-k, -k] <- a [-k, -k] - outer (a [-k, k], a [k, -k]) / p
  a [k, -k] <- -a[k, -k] / p
  a [-k, k] <- a[-k, k] / p
  return (a)
}
```
#Example Sweep



```r
set.seed(12345)
mprint (a <- matrix (sample(1:9), 3, 3))
```

```
##      [,1]  [,2]  [,3] 
## [1,]  7.00  9.00  4.00
## [2,]  8.00  3.00  2.00
## [3,]  6.00  1.00  5.00
```

```r
mprint (sweep(a, 1))
```

```
##      [,1]  [,2]  [,3] 
## [1,]  0.14 -1.29 -0.57
## [2,]  1.14 -7.29 -2.57
## [3,]  0.86 -6.71  1.57
```
#Using Sweeping for the Inverse


```r
sweepSubset <- function (a, k, verbose = FALSE) {
  for (i in k) {
    a <- sweep (a, i)
    if (verbose) mprint (a)
  }
  return (a)
}
```

```r
mprint (b <- cSweep (a, 1 : nrow (a)))
```

```
##      [,1]  [,2]  [,3] 
## [1,]  0.05  0.16  0.03
## [2,]  0.16  0.08 -0.09
## [3,]  0.03 -0.09 -0.26
```

```r
mprint (b <- fSweep (a, 1 : nrow (a)))
```

```
##      [,1]  [,2]  [,3] 
## [1,]  0.05  0.16  0.03
## [2,]  0.16  0.08 -0.09
## [3,]  0.03 -0.09 -0.26
```

```r
mprint (b <- rSweep (a, 1 : nrow (a)))
```

```
##      [,1]  [,2]  [,3] 
## [1,]  0.06 -0.20  0.03
## [2,] -0.14  0.05  0.09
## [3,] -0.05  0.23 -0.25
```

```r
mprint (a %*% b)
```

```
##      [,1]  [,2]  [,3] 
## [1,] -1.00  0.00  0.00
## [2,]  0.00 -1.00  0.00
## [3,]  0.00  0.00 -1.00
```
Using Sweeping for Regression
-------

```r
x <- matrix (rnorm (20), 10, 2)
y <- rnorm (10)
mprint (b<-qr.solve (x, y))
```

```
## [1]  0.11 -0.07
```

```r
mprint (sum((y-x%*%b) ^ 2))
```

```
## [1]  6.39
```

```r
h <- crossprod (cbind (x,y))
mprint (sweepSubset (h, 1 : 2, verbose = TRUE))
```

```
##               y    
##    0.05  0.03 -0.11
##   -0.03 13.31 -0.95
## y  0.11 -0.95  6.46
##               y    
##    0.05  0.00 -0.11
##    0.00  0.08  0.07
## y  0.11 -0.07  6.39
##               y    
##    0.05  0.00 -0.11
##    0.00  0.08  0.07
## y  0.11 -0.07  6.39
```
The partial inverse
----
Suppose \(k\) is a subset of \(\{1,2,\cdots,n\}\) and \(\overline{k}\) is the complimentary subset. In our examples we usually take, for convenience, \(k=\{1,2,\cdots,k\}\). The *partial inverse* of \(A\) on \(k\), which generalizes the simple sweep, is
\[
\mathbf{inv}_k(A)=\begin{bmatrix}
A_{kk}^{-1}&-A_{kk}^{-1}A_{k\overline{k}}\\
A_{\overline k k}A_{kk}^{-1}&A_{\overline k\overline k}-A_{\overline k k}A_{kk}^{-1}A_{k\overline{k}}
\end{bmatrix}.
\]

Why is this called a partial inverse ? If
\[
\begin{bmatrix}
A_{kk}&A_{k\overline k}\\
A_{\overline k k}&A_{\overline k\overline k}
\end{bmatrix}
\begin{bmatrix}
x\\y
\end{bmatrix}=
\begin{bmatrix}
u\\v
\end{bmatrix},
\]
then
\[
\begin{bmatrix}
A_{kk}^{-1}&-A_{kk}^{-1}A_{k\overline{k}}\\
A_{\overline k k}A_{kk}^{-1}&A_{\overline k\overline k}-A_{\overline k k}A_{kk}^{-1}A_{k\overline{k}}
\end{bmatrix}
\begin{bmatrix}
u\\y
\end{bmatrix}=
\begin{bmatrix}
x\\v
\end{bmatrix},
\]
And conversely. Of course all this assumes that \(A_{kk}\) is non-singular. The following results are true for an arbitrary index set \(k\) and its complement
\(\overline k\).


1. $A=\mathbf{inv}_k(\mathbf{inv}_k(A))$
2. $\mathbf{inv}_k(A)=(\mathbf{inv}_{\overline k}(A))^{-1}$
3. $\mathbf{inv}_{\overline k}(\mathbf{inv}_k(A))=A^{-1}$
4. $\mathbf{inv}_{\ell}(\mathbf{inv}_k(A))=\mathbf{inv}_{k}(\mathbf{inv}_\ell(A))$
5. $\mathbf{inv}_{k}(A)=\mathbf{inv}_{\overline k}(A^{-1})$
6. $(\mathbf{inv}_{k}(A))^{-1}=\mathbf{inv}_{k}(A^{-1})$

#Timing


```r
set.seed(12345)
b <- crossprod(matrix(rnorm(1000000),1000,1000))/1000
```

```
## [1] "element sweeping in pure R"
##    user  system elapsed 
##  14.543   1.397  15.155 
##    user  system elapsed 
##  14.388   1.362  14.925 
##    user  system elapsed 
##  14.351   1.352  14.872 
##    user  system elapsed 
##  14.404   1.359  14.947 
##    user  system elapsed 
##  14.441   1.393  15.046 
## [1] "element sweeping in R and fortran from R"
##    user  system elapsed 
##   1.491   0.316   1.811 
##    user  system elapsed 
##   1.353   0.294   1.650 
##    user  system elapsed 
##   1.531   0.258   1.793 
##    user  system elapsed 
##   1.522   0.258   1.783 
##    user  system elapsed 
##   1.468   0.260   1.732 
## [1] "element sweeping in R and C from R and fortran from C"
##    user  system elapsed 
##   0.459   0.010   0.471 
##    user  system elapsed 
##   0.447   0.009   0.457 
##    user  system elapsed 
##   0.456   0.010   0.467 
##    user  system elapsed 
##   0.464   0.010   0.475 
##    user  system elapsed 
##   0.456   0.010   0.468 
## [1] "just solve it"
##    user  system elapsed 
##   0.206   0.005   0.071 
##    user  system elapsed 
##   0.199   0.004   0.069 
##    user  system elapsed 
##   0.204   0.005   0.071 
##    user  system elapsed 
##   0.208   0.005   0.084 
##    user  system elapsed 
##   0.200   0.004   0.071
```

#Code

##R Code


```r
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
```

##C code


```cpp
void gsweep_(double *, double *, int *, int *, int *, int *, double *, int *);

void
hsweep(double *s, double *trian, int *n, int *m, int *ind, int *np, double *eps) {
	int ifault = 0, l = 0;
	for (int i = 0; i < *np; i++) {
		gsweep_ (s, trian, ind + i, &l, m, n, eps, &ifault);
	}
```

#NEWS

#References
