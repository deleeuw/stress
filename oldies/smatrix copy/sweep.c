
void ssweep(const int *nn, double *a, double *ppiv, int *npiv, bool *singular,
            const double *eps) {
  int n = *nn;
  double dg = 0.0;
  *singular = false;
  *ppiv = 1.0;
  *npiv = 0;
  for (int k = 1; k <= n; k++) {
    dg = a[TINDEX(k, k, n)];
    if (fabs(dg) < *eps) {
      *singular = true;
      continue;
    }
    (*ppiv) *= dg;
    (*npiv)++;
    a[TINDEX(k, k, n)] = -1 / dg;
    for (int j = 1; j <= k - 1; j++) {
      for (int i = j; i <= k - 1; i++) {
        a[TINDEX(i, j, n)] -= a[TINDEX(k, i, n)] * a[TINDEX(k, j, n)] / dg;
      }
    }
    for (int j = k + 1; j <= n; j++) {
      for (int i = j; i <= n; i++) {
        a[TINDEX(i, j, n)] -= a[TINDEX(i, k, n)] * a[TINDEX(j, k, n)] / dg;
      }
    }
    for (int j = 1; j <= k - 1; j++) {
      for (int i = k + 1; i <= n; i++) {
        a[TINDEX(i, j, n)] -= a[TINDEX(i, k, n)] * a[TINDEX(k, j, n)] / dg;
      }
    }
    for (int j = 1; j <= k - 1; j++) {
      a[TINDEX(k, j, n)] /= dg;
    }
    for (int j = k + 1; j <= n; j++) {
      a[TINDEX(j, k, n)] /= dg;
    }
  }
  return;
}

void psweep(const int *nn, double *a, double *ppiv, int *npiv, int *order,
            bool *singular, const double *eps) {
  int n = *nn, q = 0.0, itp = 0;
  double dg = 0.0, t = 0.0, tmp = 0.0, tp1 = 0.0, tp2 = 0.0, tp3 = 0.0;
  *singular = false;
  *ppiv = 1.0;
  *npiv = 0;
  for (int i = 1; i <= n; i++) {
    order[VINDEX(i)] = i;
  }
  for (int k = 1; k <= n; k++) {
    t = 0.0;
    q = 0;
    for (int i = k; i <= n; i++) {
      dg = a[TINDEX(i, i, n)];
      if (fabs(dg) > t) {
        t = dg;
        q = i;
      }
    }
    if (t < *eps) {
      *singular = true;
      return;
    }
    itp = order[VINDEX(k)];
    order[VINDEX(k)] = order[VINDEX(q)];
    order[VINDEX(q)] = itp;
    tp1 = a[TINDEX(k, k, n)];
    tp2 = a[TINDEX(q, q, n)];
    tp3 = a[TINDEX(q, k, n)];
    for (int i = 1; i <= k - 1; i++) {
      tmp = a[TINDEX(k, i, n)];
      a[TINDEX(k, i, n)] = a[TINDEX(q, i, n)];
      a[TINDEX(q, i, n)] = tmp;
    }
    for (int i = k + 1; i <= q - 1; i++) {
      tmp = a[TINDEX(i, k, n)];
      a[TINDEX(i, k, n)] = a[TINDEX(q, i, n)];
      a[TINDEX(q, i, n)] = tmp;
    }
    for (int i = q + 1; i <= n; i++) {
      tmp = a[TINDEX(i, q, n)];
      a[TINDEX(i, q, n)] = a[TINDEX(i, k, n)];
      a[TINDEX(i, k, n)] = tmp;
    }
    a[TINDEX(q, q, n)] = tp1;
    a[TINDEX(k, k, n)] = tp2;
    a[TINDEX(q, k, n)] = tp3;
    (*ppiv) *= t;
    (*npiv)++;
    a[TINDEX(k, k, n)] = -1 / t;
    for (int j = 1; j <= k - 1; j++) {
      for (int i = j; i <= k - 1; i++) {
        a[TINDEX(i, j, n)] -= a[TINDEX(k, i, n)] * a[TINDEX(k, j, n)] / t;
      }
    }
    for (int j = k + 1; j <= n; j++) {
      for (int i = j; i <= n; i++) {
        a[TINDEX(i, j, n)] -= a[TINDEX(i, k, n)] * a[TINDEX(j, k, n)] / t;
      }
    }
    for (int j = 1; j <= k - 1; j++) {
      for (int i = k + 1; i <= n; i++) {
        a[TINDEX(i, j, n)] -= a[TINDEX(i, k, n)] * a[TINDEX(k, j, n)] / t;
      }
    }
    for (int j = 1; j <= k - 1; j++) {
      a[TINDEX(k, j, n)] /= t;
    }
    for (int j = k + 1; j <= n; j++) {
      a[TINDEX(j, k, n)] /= t;
    }
  }
  return;
}
