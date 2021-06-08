
void trimat(const int *n, const double *x, double *y) {
  int nn = *n;
  for (int i = 1; i <= nn; i++) {
    for (int j = 1; j <= nn; j++) {
      y[MINDEX(i, j, nn)] =
          (i >= j) ? x[TINDEX(i, j, nn)] : x[TINDEX(j, i, nn)];
    }
  }
  return;
}

void tritri(const int *n, const double *x, double *y) {
  int nn = *n;
  for (int i = 1; i <= nn; i++) {
    for (int j = 1; j <= nn; j++) {
      y[MINDEX(i, j, nn)] = (i >= j) ? x[TINDEX(i, j, nn)] : 0.0;
    }
  }
  return;
}

void mattri(const int *n, const double *x, double *y) {
  int nn = *n;
  for (int j = 1; j <= nn; j++) {
    for (int i = j; i <= nn; i++) {
      y[TINDEX(i, j, nn)] = x[MINDEX(i, j, nn)];
    }
  }
  return;
}

void primat(const int *n, const int *m, const int *w, const int *p,
            const double *x) {
  for (int i = 1; i <= *n; i++) {
    for (int j = 1; j <= *m; j++) {
#ifdef USING_R
      Rprintf(" %*.*f ", *w, *p, x[MINDEX(i, j, *n)]);
#else
      printf(" %*.*f ", *w, *p, x[MINDEX(i, j, *n)]);
#endif
    }
#ifdef USING_R
    Rprintf("\n");
#else
    printf("\n");
#endif
  }
#ifdef USING_R
  Rprintf("\n\n");
#else
  printf("\n\n");
#endif
  return;
}

void pritru(const int *n, const int *w, const int *p, const double *x) {
  for (int i = 1; i <= *n; i++) {
    for (int j = 1; j <= i; j++) {
#ifdef USING_R
      Rprintf(" %*.*f ", *w, *p, x[TINDEX(i, j, *n)]);
#else
      printf(" %*.*f ", *w, *p, x[TINDEX(i, j, *n)]);
#endif
    }
#ifdef USING_R
    Rprintf("\n");
#else
    printf("\n");
#endif
  }
#ifdef USING_R
  Rprintf("\n\n");
#else
  printf("\n\n");
#endif
  return;
}

void prisru(const int *n, const int *m, const int *w, const int *p,
            const double *x) {
  for (int k = 1; k <= *m; k++) {
    for (int i = 1; i <= *n; i++) {
      for (int j = 1; j <= i; j++) {
#ifdef USING_R
        Rprintf(" %*.*f ", *w, *p, x[UINDEX(i, j, k, *n, *m)]);
#else
        printf(" %*.*f ", *w, *p, x[UINDEX(i, j, k, *n, *m)]);
#endif
      }
#ifdef USING_R
      Rprintf("\n");
#else
      printf("\n");
#endif
    }
#ifdef USING_R
    Rprintf("\n\n");
#else
    printf("\n\n");
#endif
  }
  return;
}
