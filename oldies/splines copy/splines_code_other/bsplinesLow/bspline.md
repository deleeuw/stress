Low Order B-Splines
========================================================

We start with a sequence of knots to divide the real line into intervals. A *spline* of order `k` is a polynomial of degree `k` in each interval, but with the  property that at the knots the spline has continuous derivatives of degree `k-1`. A B-spline of order `k` has the additional property that it is non-zero in only `k+1` consecutive intervals. B-splines take values between zero and one, and their integral over the real line is one.

B-splines of arbitrary degree are normally computed using convenient and stable recursion formulas, but we give code for some explicit computation of low-order B-splines below. First we define some constants.


```r
set.seed (12345)
knots <- sort (rnorm (20))
knots
```

```
##  [1] -1.8180 -0.9193 -0.8864 -0.7505 -0.4535 -0.3316 -0.2842 -0.2762
##  [9] -0.1162 -0.1093  0.2987  0.3706  0.5202  0.5855  0.6059  0.6301
## [17]  0.7095  0.8169  1.1207  1.8173
```

```r
grid <- seq (-1.5, 1.5, length = 1001)
```



```
## Loading fonts...
## Loading fonts finished
```


Degree zero B-splines
-------
A B-spline of degree zero is a step function which is zero everywhere, except on one of the intervals, where it is one.

```r
ZbSpline <- function (x, knots, k = 1) {
  ZbSplineSingle <- function (x, knots, k = 1) {
	k0 <- knots [k]
	k1 <- knots [k + 1]
	if ((x >= k0) && (x < k1)) {
		return (1)
	}
	return (0)
	}
	return (sapply (x, function (z) ZbSplineSingle (z, knots, k)))
}
```


















