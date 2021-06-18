#include "../include/smacof.h"

// Euclidean distances

void smacofDist(const double *x, const int n, const int p, double *dist) {
  for (int j = 1; j <= (n - 1); j++) {
    for (int i = j + 1; i <= n; i++) {
      double d = 0.0;
      for (int s = 1; s <= p; s++) {
        d += SQUARE(x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]);
      }
      dist[SINDEX(i, j, n)] = sqrt(d);
    }
  }
  return;
}

// double center symmetrtc hollow matrix in strict lower triangular storage
// diagonal elements of the doubly centered matrix are left implicit

void smacofDoubleCenter(double *x, const int n) {
  double *rsums = (double *)calloc((size_t)n, sizeof(double));
  double s = 0.0, t = 0.0;
  for (int i = 1; i <= n; i++) {
    s = 0.0;
    for (int j = 1; j <= n; j++) {
      if (i > j) {
        s += x[SINDEX(i, j, n)];
        t += x[SINDEX(i, j, n)];
      }
      if (i < j) {
        s += x[SINDEX(j, i, n)];
        t += x[SINDEX(j, i, n)];
      }
    }
    rsums[i] = s / ((double)n);
  }
  t /= ((double)SQUARE(n));
  for (int j = 1; j <= (n - 1); j++) {
    for (int i = j + 1; i <= n; i++) {
      x[SINDEX(i, j, n)] = x[SINDEX(i, j, n)] - rsums[i] - rsums[j] + t;
    }
  }
  free(rsums);
  return;
}

// multiply doubly centered symmetric matrix stored in strictly lower
// triangular form (without diagonal) with a general matrix

void smacofDCMultX(const double *dcmat, double *x, const int n, const int p) {
  double t = 0.0, u = 0.0;
  double *v = (double *)calloc((size_t)n, sizeof(double));
  for (int s = 1; s <= p; s++) {
    for (int i = 1; i <= n; i++) {
      t = 0.0;
      for (int j = 1; j <= n; j++) {
        u = x[MINDEX(j, s, n)] - x[MINDEX(i, s, n)];
        if (i > j) {
          t += dcmat[SINDEX(i, j, n)] * u;
        }
        if (i < j) {
          t += dcmat[SINDEX(j, i, n)] * u;
        }
      }
      v[VINDEX(i)] = t;
    }
    for (int i = 1; i <= n; i++) {
      x[MINDEX(i, s, n)] = v[VINDEX(i)];
    }
  }
  free(v);
  return;
}
