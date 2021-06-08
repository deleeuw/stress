Sweeping and the partial inverse
========================================================

Simple Sweep
------
Suppose \(A\) is a square non-singular matrix. A *simple sweep* on the \(k^{th}\) element replaces \(A\) by the matrix with elements
\[
\begin{align*}
a_{kk}&=\frac{1}{a_{kk}},\\
a_{ki}&=-\frac{a_{ki}}{a_{kk}},\\
a_{ik}&=\frac{a_{ik}}{a_{kk}},\\
a_{ij}&=a_{ij}-\frac{a_{ik}a_{kj}}{a_{kk}}
\end{align*}
\]
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

Example Sweep
-----



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

Using Sweeping for the Inverse
-------

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
mprint (b <- sweepSubset (a, 1 : nrow (a), verbose = TRUE))
```

```
##      [,1]  [,2]  [,3] 
## [1,]  0.14 -1.29 -0.57
## [2,]  1.14 -7.29 -2.57
## [3,]  0.86 -6.71  1.57
##      [,1]  [,2]  [,3] 
## [1,] -0.06  0.18 -0.12
## [2,]  0.16 -0.14 -0.35
## [3,] -0.20  0.92  3.94
##      [,1]  [,2]  [,3] 
## [1,] -0.06  0.20 -0.03
## [2,]  0.14 -0.05 -0.09
## [3,]  0.05 -0.23  0.25
##      [,1]  [,2]  [,3] 
## [1,] -0.06  0.20 -0.03
## [2,]  0.14 -0.05 -0.09
## [3,]  0.05 -0.23  0.25
```

```r
mprint (a %*% b)
```

```
##      [,1]  [,2]  [,3] 
## [1,]  1.00 -0.00  0.00
## [2,]  0.00  1.00  0.00
## [3,] -0.00  0.00  1.00
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

### First
\[
A = \mathbf{inv}_k(\mathbf{inv}_k(A))
\]

```r
mprint (a - sweepSubset (sweepSubset (a, 1 : 2), 1 : 2))
```

```
##      [,1]  [,2]  [,3] 
## [1,] -0.00 -0.00 -0.00
## [2,] -0.00 -0.00 -0.00
## [3,] -0.00 -0.00  0.00
```

### Second
\[
\mathbf{inv}_k(A)=(\mathbf{inv}_{\overline k}(A))^{-1}
\]

```r
mprint (sweepSubset (a, 1 : 2) - solve (sweepSubset (a, 3)))
```

```
##      [,1]  [,2]  [,3] 
## [1,]  0.00 -0.00 -0.00
## [2,] -0.00  0.00  0.00
## [3,]  0.00 -0.00 -0.00
```

### Third
\[
\mathbf{inv}_{\overline k}(\mathbf{inv}_k(A))=A^{-1}
\]

```r
mprint (solve (a) - sweepSubset (sweepSubset (a, 1 : 2), 3))
```

```
##      [,1]  [,2]  [,3] 
## [1,] -0.00  0.00  0.00
## [2,]  0.00  0.00 -0.00
## [3,]  0.00 -0.00  0.00
```

### Fourth
\[
\mathbf{inv}_{\ell}(\mathbf{inv}_k(A))=\mathbf{inv}_{k}(\mathbf{inv}_\ell(A))
\]

```r
mprint (sweepSubset (sweepSubset (a, 1 : 2), 3) - sweepSubset (sweepSubset (a, 3), 1 : 2))
```

```
##      [,1]  [,2]  [,3] 
## [1,]  0.00 -0.00 -0.00
## [2,] -0.00  0.00  0.00
## [3,] -0.00  0.00  0.00
```

### Fifth
\[
\mathbf{inv}_{k}(A)=\mathbf{inv}_{\overline k}(A^{-1})
\]

```r
mprint (sweepSubset (a, 1 : 2) - sweepSubset (solve (a), 3))
```

```
##      [,1]  [,2]  [,3] 
## [1,]  0.00 -0.00 -0.00
## [2,] -0.00  0.00  0.00
## [3,]  0.00 -0.00  0.00
```

### Sixth
\[
(\mathbf{inv}_{k}(A))^{-1}=\mathbf{inv}_{k}(A^{-1})
\]

```r
mprint (solve (sweepSubset (a, 1 : 2)) - sweepSubset (solve (a), 1 : 2))
```

```
##      [,1]  [,2]  [,3] 
## [1,]  0.00  0.00  0.00
## [2,]  0.00  0.00  0.00
## [3,] -0.00 -0.00 -0.00
```

