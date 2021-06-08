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
