
void invtri(const int *nn, double *x, const double *eps, bool *ierr) {
  int n = *nn;
  double s = 0.0, t = 0.0;
  for (int k = 1; k <= n; k++) {
    t = x[TINDEX(k, k, n)];
    if (fabs(t) < *eps) {
      *ierr = true;
      return;
    }
    x[TINDEX(k, k, n)] = 1.0 / t;
    for (int i = k + 1; i <= n; i++) {
      s = 0.0;
      for (int j = k; j < i; j++) {
        s += x[TINDEX(i, j, n)] * x[TINDEX(j, k, n)];
      }
      x[TINDEX(i, k, n)] = -s / x[TINDEX(i, i, n)];
    }
  }
  return;
}