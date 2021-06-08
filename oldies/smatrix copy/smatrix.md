---
title: "R/C Routines for Symmetric Matrices in Lower Triangular Column-Major Storage"
author: "Jan de Leeuw"
date: "Version June 27, 2018"
output:
  html_document:
    keep_md: yes
    number_sections: yes
    toc: yes
  word_document:
    toc: yes
  pdf_document:
    keep_tex: yes
    number_sections: yes
    toc: yes
fontsize: 12pt
graphics: yes
bibliography: smatrix.bib
abstract: dodo
---
<style type="text/css">

body{ /* Normal  */
   font-size: 18px;
}
td {  /* Table  */
   font-size: 18px;
}
h1 { /* Header 1 */
 font-size: 28px;
 color: DarkBlue;
}
h2 { /* Header 2 */
 font-size: 22px;
 color: DarkBlue;
}
h3 { /* Header 3 */
 font-size: 18px;
 color: DarkBlue;
}
code.r{ /* Code block */
  font-size: 18px;
}
pre { /* Code block */
  font-size: 18px
}
</style>

<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  TeX: { equationNumbers: { autoNumber: "AMS" } }
});
</script>




**Note:** This is a working paper which will be expanded/updated frequently. All suggestions for improvement are welcome. The directory [gifi.stat.ucla.edu/smatrix](http://gifi.stat.ucla.edu/smatrix) has a pdf version, the bib file, the R and C code, and the complete Rmd file.

# Introduction

Suppose $A$ is a symmetric matrix, say

```
##     2     3     4     5     6 
##     3     4     5     6     7 
##     4     5     6     7     8 
##     5     6     7     8     9 
##     6     7     8     9    10
```
We store the elements on the matrix as a vector, using the lower triangle (including the diagonal) with the first column first, then the second column, and so on. In R this is simply

```r
a[outer(1:5,1:5,">=")]
```

```
##  [1]  2  3  4  5  6  4  5  6  7  6  7  8  8  9 10
```
This defines a symmetric matrix $A$ in LTCM (lower triangle column-major) storage mode. Our job in this paper is to rewrite some of the basic linear algebra algorithms for symmetric matrices so that they apply to LTCM matrices. Clearly this implies we save memory, but since accessing elements will be more complicated we sacrifice some speed. Much of this is available in LAPACK/LINPACK and probably in many other libraries, but I am just interested in rewriting some of the code in my own dialect. 

# Conventions

Our R functions are wrappers for C functions. Since we want to be able to use our C functions outside the R environment as well, they should use a few of the R resources as possible. Thus we use the .C() interface, and we try to avoid memory allocation and I/O
in the C functions as much as possible. We leave that to the caller, which is either R or a main routine in C. By the .C() conventions all C functions return void, and pass their arguments by reference (have pointers as arguments). 

To find element $(i,j)$ in a matrix we use C functions `MINDEX, TINDEX, SINDEX` to compute the location of the element in the 
linear storage. `MINDEX` is for full matrices, `TINDEX` for LTCM matrices, and `SINDEX` for strictly lower triangular matrices. 
In our routines row and column indices start at one, not at zero. Thus we also use `VINDEX` to access elements of a vector,
so that in C `a[VINDEX(i)] = a[i-1]`. To complete the set we also have `UINDEX`, which is used for a sequence of LTCM matrices
in linear storage. Because these routines are used heavily, we define them to be static inline in the header file `smatrix.h`.

Note that our conventions to start indices at one, and to use column-major storage, are both chosen to conform to R (and thus, ultimately, to FORTRAN). 

# Routines

## cholesky

The routine `scholesky()` decomposes a symmetric positive semi-definite matrix $A$ as $A=LL'$, with $L$ lower-triangular.
It constructs the columns of $L$ one by one, in order, and overwrites $A$ with $L$. Along the way it computes the determinant
as the product of the pivots. If it encounters a zero pivot, which can only happen if the matrix is singular, it just jumps over it.

```r
x<-matrix(c(1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,-1,-1),4,4)
a<-tcrossprod(x[,c(1,4)])[outer(1:4,1:4,">=")]
mprint(triangle2matrix(a), w = 3, d = 0)
```

```
##      [,1] [,2] [,3] [,4]
## [1,]   2    2    0    0 
## [2,]   2    2    0    0 
## [3,]   0    0    2    2 
## [4,]   0    0    2    2
```

```r
h<-symmetricCholesky(a, eps = 1e-10)
g<-.C("pritru", as.integer(4), as.integer(8), as.integer(4), h[[2]])
```

```
##    1.4142 
##    1.4142    0.0000 
##    0.0000    0.0000    1.4142 
##    0.0000    0.0000    1.4142    0.0000
```
Cholesky with pivots is implemented in `pcholesky()`.

```r
h<-pivotCholesky(a, eps = 1e-10)
g<-.C("pritru", as.integer(4), as.integer(8), as.integer(4), h[[2]])
```

```
##    1.4142 
##    0.0000    1.4142 
##    0.0000    1.4142    0.0000 
##    1.4142    0.0000    0.0000    0.0000
```
7.8886090522\times 10^{-31} 
4 
4, 1, 2, 3


```r
b<-(x %*% (c(1,-1,1,1) * t(x)))[outer(1:4,1:4,">=")]
#mprint(triangle2matrix(b), w = 3, d = 0)
#h<-pivotCholesky(b, eps = 1e-10)
#h

#mprint(triangle2matrix(h[[2]]),  w = 4, d = 2, f ="+")
```
≈
## sweep

The `ssweep()` implements a sequence of the symmetric version of the Beaton sweep (@beaton_64, @goodnight_79, @deleeuw_E_16j).
Given our storage mode, it is obviously important to choose sweeps that preserve symmetry.

After a full complement of $n$ symmetric sweeps a non-singular symmetric matrix $A$ is transformed to $-A^{-1}$. Also, the determinant and rank are returned.

If we encounter a zero diagonal element along the way, then the sweep is skipped. The resulting returned matrix will contain a partial inverse of the swept rows and columns. See @lange_10 (p 95-99) or @deleeuw_E_15. For positive semi-definite $A$ the principal sub-matrix of non-swept elements is still positive semi-definite, so a zero diagonal elements means the whole row and column in the sub-matrix are zero. Pivoting would make the situation more clear, but I have not implemented it yet.


```r
x<-matrix(c(1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,-1,-1),4,4)
a<-tcrossprod(x[,1:2])[outer(1:4,1:4,">=")]
mprint(triangle2matrix(a), w = 3, d = 0)
```

```
##      [,1] [,2] [,3] [,4]
## [1,]   2    0    0    2 
## [2,]   0    2    2    0 
## [3,]   0    2    2    0 
## [4,]   2    0    0    2
```

```r
h<-symmetricSweep(a, eps = 1e-10)
mprint(triangle2matrix(h[[2]]),  w = 4, d = 2, f ="+")
```

```
##      [,1]  [,2]  [,3]  [,4] 
## [1,] -0.50 +0.00 +0.00 +1.00
## [2,] +0.00 -0.50 +1.00 +0.00
## [3,] +0.00 +1.00 +0.00 +0.00
## [4,] +1.00 +0.00 +0.00 +0.00
```
The determinant is 4 and the rank is 2.


```r
a<-(outer(x[,1],x[,1])-outer(x[,2],x[,2]))[outer(1:4,1:4,">=")]
mprint(triangle2matrix(a), w = 3, d = 0)
```

```
##      [,1] [,2] [,3] [,4]
## [1,]   0    2    2    0 
## [2,]   2    0    0    2 
## [3,]   2    0    0    2 
## [4,]   0    2    2    0
```

```r
h<-symmetricSweep(a, eps = 1e-10)
mprint(triangle2matrix(h[[2]]),  w = 4, d = 2, f ="+")
```

```
##      [,1]  [,2]  [,3]  [,4] 
## [1,] +0.00 +2.00 +2.00 +0.00
## [2,] +2.00 +0.00 +0.00 +2.00
## [3,] +2.00 +0.00 +0.00 +2.00
## [4,] +0.00 +2.00 +2.00 +0.00
```

```r
x<-outer(1:4,0:3,"^")
a<-crossprod(x, diag(c(1,-1,-1,-1))*t(x))[outer(1:4,1:4,">=")]
mprint(triangle2matrix(a), w = 6, d = 0)
```

```
##      [,1]   [,2]   [,3]   [,4]  
## [1,]      1      1      1      1
## [2,]      1     -4     -8    -16
## [3,]      1     -8    -81   -243
## [4,]      1    -16   -243  -4096
```

```r
h<-symmetricSweep(a, eps = 1e-10)
mprint(triangle2matrix(h[[2]]),  w = 10, d = 8, f ="+")
```

```
##      [,1]        [,2]        [,3]        [,4]       
## [1,] -0.79026225 -0.22202617 +0.01234654 -0.00005812
## [2,] -0.22202617 +0.25101548 -0.02971761 +0.00072830
## [3,] +0.01234654 -0.02971761 +0.01834001 -0.00096894
## [4,] -0.00005812 +0.00072830 -0.00096894 +0.00029877
```

```r
triangle2matrix(a)%*%triangle2matrix(h[[2]])
```

```
##                  [,1]             [,2]             [,3]             [,4]
## [1,] -1.000000000e+00  6.938893904e-18  3.469446952e-18  0.000000000e+00
## [2,]  0.000000000e+00 -1.000000000e+00  2.775557562e-17 -4.336808690e-19
## [3,]  0.000000000e+00  0.000000000e+00 -1.000000000e+00  5.204170428e-18
## [4,] -4.440892099e-16 -8.881784197e-16  7.771561172e-16 -1.000000000e+00
```
The determinant is -1.101199\times 10^{6} and the rank is 4.


```r
h<-pivotSweep(a, eps = 1e-10)
mprint(triangle2matrix(h[[2]]), w = 4, d = 2, f ="+")
```

```
##      [,1]  [,2]   [,3]    [,4]    
## [1,] +1.00 +1.00  +1.00   +1.00   
## [2,] +1.00 -4.00  -8.00   -16.00  
## [3,] +1.00 -8.00  -81.00  -243.00 
## [4,] +1.00 -16.00 -243.00 -4096.00
```
The determinant is 1 and the rank is 0.

## invtri

## jacobi

If $A_k$ are $m$ symmetric matrices, then compute $K$ with $K'K=KK'=I$ such that the sum of squares of the diagonal elements of the $H_k=K'A_kK$ is as large as possible (and the sum of squares of the off-diagonal elements is as small as possible). This is
basically the routine discussed in more detail in @deleeuw_E_18g.

This is the special case with $m=1$ if smjacobi, i.e. the classical Jacobi method to compute the eigen decomposition of a single symmetric matrix, previously implemented in this form by @deleeuw_E_17o. 

## smpinverse

If $A=K\Lambda K'$ is the eigen decomposition of $A$ then this routine computes $K\Phi K'$, where $\phi_i=f(\lambda_i)$, where $f$ can be  (type argument equal to one, two, or three)
$$
f_1(\lambda)=\begin{cases}\lambda^{-1}&\text{ if }\lambda\not= 0,\\0&\text{ otherwise}.\end{cases}
$$
$$
f_2(\lambda)=\begin{cases}\sqrt{\lambda}&\text{ if }\lambda\geq 0,\\0&\text{ otherwise}.\end{cases}
$$
$$
f_3(\lambda)=\begin{cases}1/\sqrt{\lambda}&\text{ if }\lambda\geq 0,\\0&\text{ otherwise}.\end{cases}
$$
# Code

## R code

### smatrix.R


```r
dyn.load("smatrix.so")

symmetricCholesky <- function (a, eps = 1e-10) {
  m <- length (a)
  n <- (sqrt (1 + 8 * m) - 1) / 2
  h <- .C(
    "scholesky",
    n = as.integer(n),
    a = as.double (a),
    singular = as.integer(0),
    indefinite = as.integer(0),
    det = as.double (0),
    eps = as.double(eps)
  )
  return (h)
}

pivotCholesky <- function (a, eps = 1e-10) {
  m <- length (a)
  n <- (sqrt (1 + 8 * m) - 1) / 2
  h <- .C(
    "pcholesky",
    n = as.integer(n),
    a = as.double (a),
    det = as.double (0),
    rank = as.integer(0),
    order = as.integer(1:n),
    indefinite = as.integer(0),
    eps = as.double(eps)
  )
  return (h)
}

symmetricSweep <- function (a, eps = 1e-10) {
  m <- length (a)
  n <- (sqrt (1 + 8 * m) - 1) / 2
  h <- .C(
    "ssweep",
    n = as.integer(n),
    a = as.double (a),
    det = as.double (0),
    rank = as.integer(0),
    singular = as.integer(0),
    eps = as.double (eps)
  )
  return (h)
}

pivotSweep <- function (a, eps = 1e-10) {
  m <- length (a)
  n <- (sqrt (1 + 8 * m) - 1) / 2
  h <- .C(
    "psweep",
    n = as.integer(n),
    a = as.double (a),
    det = as.double (0),
    rank = as.integer(0),
    order = as.integer(1:n),
    singular = as.integer(0),
    eps = as.double (eps)
  )
  return (h)
}


triangle2matrix <- function (x) {
  m <- length (x)
  n <- round ((sqrt (1 + 8 * m) - 1) / 2)
  h <-
    .C("trimat", as.integer (n), as.double (x), as.double (rep (0, n * n)))
  return (matrix(h[[3]], n, n))
}

matrix2triangle <- function (x) {
  n <- dim(x)[1]
  m <- n * (n + 1) / 2
  h <-
    .C("mattri", as.integer (n), as.double (x), as.double (rep (0, m)))
  return (h[[3]])
}

triangle2triangle <- function (x) {
  m <- length (x)
  n <- round ((sqrt (1 + 8 * m) - 1) / 2)
  h <-
    .C("tritri", as.integer (n), as.double (x), as.double (rep (0, n * n)))
  return (matrix(h[[3]], n, n))
}

matrix2print <- function (x, w = 15, p = 10) {
  n <- nrow (x)
  m <- ncol (x)
  h <-
    .C(
      "primat",
      as.integer(n),
      as.integer(m),
      as.integer(w),
      as.integer(p),
      as.double (x)
    )
}
```

## C Code

### smatrix.h


```c
#ifndef SMATRIX_H
#define SMATRIX_H

#define USING_R

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef USING_R
#include <R.h>
#endif

static inline int VINDEX(const int);
static inline int MINDEX(const int, const int, const int);
static inline int SINDEX(const int, const int, const int);
static inline int TINDEX(const int, const int, const int);
static inline int AINDEX(const int, const int, const int, const int, const int);
static inline int UINDEX(const int, const int, const int, const int, const int);

static inline double SQUARE(const double);
static inline double THIRD(const double);
static inline double FOURTH(const double);
static inline double FIFTH(const double);

static inline double MAX(const double, const double);
static inline double MIN(const double, const double);
static inline int IMIN(const int, const int);
static inline int IMAX(const int, const int);

// index function for vector

static inline int VINDEX(const int i) { return i - 1; }

// index function for rectangular matrix

static inline int MINDEX(const int i, const int j, const int n) {
  return (i - 1) + (j - 1) * n;
}

// index function for symmetric matrix

static inline int AINDEX(const int i, const int j, const int k, const int n,
                         const int m) {
  return (i - 1) + (j - 1) * n + (k - 1) * n * m;
}

// index function for lower triangle without diagonal

static inline int SINDEX(const int i, const int j, const int n) {
  return ((j - 1) * n) - (j * (j - 1) / 2) + (i - j) - 1;
}

// index function for lower triangle with diagonal

static inline int TINDEX(const int i, const int j, const int n) {
  return ((j - 1) * n) - ((j - 1) * (j - 2) / 2) + (i - (j - 1)) - 1;
}

// index function for multiple lower triangle with diagonal

static inline int UINDEX(const int i, const int j, const int k, const int n,
                         const int m) {
  return ((k - 1) * n * (n + 1) / 2) + ((j - 1) * n) - ((j - 1) * (j - 2) / 2) +
         (i - (j - 1)) - 1;
}

static inline double SQUARE(const double x) { return x * x; }
static inline double THIRD(const double x) { return x * x * x; }
static inline double FOURTH(const double x) { return x * x * x * x; }
static inline double FIFTH(const double x) { return x * x * x * x * x; }

static inline double MAX(const double x, const double y) {
  return (x > y) ? x : y;
}

static inline double MIN(const double x, const double y) {
  return (x < y) ? x : y;
}

static inline int IMAX(const int x, const int y) { return (x > y) ? x : y; }

static inline int IMIN(const int x, const int y) { return (x < y) ? x : y; }

#endif /* SMATRIX_H */
```

### smatrix.c


```c

#include "smatrix.h"
#include "sroutines.h"

#include "cholesky.c"
#include "invtri.c"
#include "jacobi.c"
#include "sweep.c"
#include "utility.c"
```



# References
