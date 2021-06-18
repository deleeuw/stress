#include "../include/smacof.h"

void smacofInitialC(const double *delta, const int *nobjects,
                    const int *ndimensions, double *x) {
  int n = *nobjects, p = *ndimensons, itmax = 100, r = n * (n + 1) / 2;
  double s, ss = 0.0, eps = 1e-6;
  double *work1 = (double *)calloc((size_t)n, sizeof(double));
  double *work2 = (double *)calloc((size_t)r, sizeof(double));
  double *work3 = (double *)calloc((size_t)n * n, sizeof(double));
  double *work4 = (double *)calloc((size_t)n, sizeof(double));
  for (int i = 1; i <= n; i++) {
    work1[VINDEX(i)] = 0.0;
    for (int j = 1; j <= n; j++) {
      if (i == j)
        continue;
      int ij = IMAX(i, j);
      int ji = IMIN(i, j);
      s = SQUARE(delta[SINDEX(ij, ji, n)]);
      ss += s;
      work1[VINDEX(i)] += s;
      if (j < i)
        continue;
      work2[TINDEX(ij, ji, n)] = s;
    }
    work1[VINDEX(i)] /= (double)n;
  }
  ss /= SQUARE((double)n);
  for (int j = 1; j <= n; j++) {
    for (int i = j; i <= n; i++) {
      work2[TINDEX(i, j, n)] -= work1[VINDEX(i)];
      work2[TINDEX(i, j, n)] -= work1[VINDEX(j)];
      work2[TINDEX(i, j, n)] += ss;
      work2[TINDEX(i, j, n)] *= -0.5;
    }
  }
  (void)jacobiC(n, work2, work3, work1, work4, &itmax, &eps);
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= p; j++) {
      s = work2[TINDEX(j, j, n)];
      if (s <= 0)
        continue;
      x[MINDEX(i, j, n)] = work3[MINDEX(i, j, n)] * sqrt(s);
    }
  }
  free work1;
  free work2;
  free work3;
  free work4;
  return;
}
