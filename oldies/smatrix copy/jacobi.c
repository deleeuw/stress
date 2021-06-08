
void smjacobi(const int *nn, const int *mm, double *a, double *evec,
              double *fini, double *f, int *itel, const int *itmax,
              const double *eps, const bool *verbose, const bool *vectors) {
  int n = *nn, m = *mm;
  double d = 0.0, cost = 0.0, sint = 0.0, u = 0.0, p = 0.0, q = 0.0, r = 0.0,
         piil = 0.0, pijl = 0.0, lbd = 0.0, dd = 0.0, pp = 0.0, dp = 0.0;
  double fold = 0.0, fnew = 0.0, oldi = 0.0, oldj = 0.0;
  *itel = 1;
  if (*vectors) {
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
        evec[MINDEX(i, j, n)] = (i == j) ? 1.0 : 0.0;
      }
    }
  }
  for (int k = 1; k <= m; k++) {
    for (int i = 1; i <= n; i++) {
      fold += SQUARE(a[UINDEX(i, i, k, n, m)]);
    }
  }
  *fini = fold;
  while (true) {
    for (int j = 1; j <= n - 1; j++) {
      for (int i = j + 1; i <= n; i++) {
        dd = 0.0, pp = 0.0, dp = 0.0;
        for (int k = 1; k <= m; k++) {
          p = a[UINDEX(i, j, k, n, m)];
          q = a[UINDEX(i, i, k, n, m)];
          r = a[UINDEX(j, j, k, n, m)];
          d = (q - r) / 2.0;
          dd += SQUARE(d);
          pp += SQUARE(p);
          dp += p * d;
        }
        lbd = ((dd + pp) - sqrt(SQUARE(dd - pp) + 4.0 * SQUARE(dp))) / 2.0;
        u = dp / sqrt(SQUARE(dp) + SQUARE(lbd - pp));
        if ((fabs(dp) < 1e-15) && (pp <= dd)) {
          continue;
        }
        cost = sqrt((1 + u) / 2);
        sint = -sqrt((1 - u) / 2);
        if (*vectors) {
          for (int l = 1; l <= n; l++) {
            piil = evec[MINDEX(l, i, n)];
            pijl = evec[MINDEX(l, j, n)];
            evec[MINDEX(l, i, n)] = cost * piil - sint * pijl;
            evec[MINDEX(l, j, n)] = sint * piil + cost * pijl;
          }
        }
        for (int k = 1; k <= m; k++) {
          p = a[UINDEX(i, j, k, n, m)];
          q = a[UINDEX(i, i, k, n, m)];
          r = a[UINDEX(j, j, k, n, m)];
          for (int l = 1; l <= n; l++) {
            if ((l == i) || (l == j))
              continue;
            int il = IMIN(i, l);
            int li = IMAX(i, l);
            int jl = IMIN(j, l);
            int lj = IMAX(j, l);
            oldi = a[UINDEX(li, il, k, n, m)];
            oldj = a[UINDEX(lj, jl, k, n, m)];
            a[UINDEX(li, il, k, n, m)] = cost * oldi - sint * oldj;
            a[UINDEX(lj, jl, k, n, m)] = sint * oldi + cost * oldj;
          }
          a[UINDEX(i, i, k, n, m)] =
              SQUARE(cost) * q + SQUARE(sint) * r - 2 * cost * sint * p;
          a[UINDEX(j, j, k, n, m)] =
              SQUARE(sint) * q + SQUARE(cost) * r + 2 * cost * sint * p;
          a[UINDEX(i, j, k, n, m)] =
              cost * sint * (q - r) + (SQUARE(cost) - SQUARE(sint)) * p;
        }
      }
    }
    fnew = 0.0;
    for (int k = 1; k <= m; k++) {
      for (int i = 1; i <= n; i++) {
        fnew += SQUARE(a[UINDEX(i, i, k, n, m)]);
      }
    }
    if (*verbose == true) {
#ifdef USING_R
      Rprintf("itel %4d fold %15.10f fnew %15.10f\n", *itel, fold, fnew);
#else
      printf("itel %4d fold %15.10f fnew %15.10f\n", *itel, fold, fnew);
#endif
    }
    if (((fnew - fold) < *eps) || (*itel == *itmax))
      break;
    fold = fnew;
    (*itel)++;
  }
  *f = fnew;
  return;
}

void sjacobi(const int *nn, double *a, double *evec, double *fini, double *f,
             int *itel, const int *itmax, const double *eps,
             const bool *verbose, const bool *vectors) {
  int n = *nn;
  double d = 0.0, cost = 0.0, sint = 0.0, u = 0.0, p = 0.0, q = 0.0, r = 0.0,
         piil = 0.0, pijl = 0.0, lbd = 0.0, dd = 0.0, pp = 0.0, dp = 0.0;
  double fold = 0.0, fnew = 0.0, oldi = 0.0, oldj = 0.0;
  *itel = 1;
  if (*vectors) {
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
        evec[MINDEX(i, j, n)] = (i == j) ? 1.0 : 0.0;
      }
    }
  }
  for (int i = 1; i <= n; i++) {
    fold += SQUARE(a[TINDEX(i, i, n)]);
  }
  *fini = fold;
  while (true) {
    for (int j = 1; j <= n - 1; j++) {
      for (int i = j + 1; i <= n; i++) {
        dd = 0.0, pp = 0.0, dp = 0.0;
        p = a[TINDEX(i, j, n)];
        q = a[TINDEX(i, i, n)];
        r = a[TINDEX(j, j, n)];
        d = (q - r) / 2.0;
        dd = SQUARE(d);
        pp = SQUARE(p);
        dp = p * d;
        lbd = ((dd + pp) - sqrt(SQUARE(dd - pp) + 4.0 * SQUARE(dp))) / 2.0;
        u = dp / sqrt(SQUARE(dp) + SQUARE(lbd - pp));
        if ((fabs(dp) < 1e-15) && (pp <= dd)) {
          continue;
        }
        cost = sqrt((1 + u) / 2);
        sint = -sqrt((1 - u) / 2);
        if (*vectors) {
          for (int l = 1; l <= n; l++) {
            piil = evec[MINDEX(l, i, n)];
            pijl = evec[MINDEX(l, j, n)];
            evec[MINDEX(l, i, n)] = cost * piil - sint * pijl;
            evec[MINDEX(l, j, n)] = sint * piil + cost * pijl;
          }
        }
        p = a[TINDEX(i, j, n)];
        q = a[TINDEX(i, i, n)];
        r = a[TINDEX(j, j, n)];
        for (int l = 1; l <= n; l++) {
          if ((l == i) || (l == j))
            continue;
          int il = IMIN(i, l);
          int li = IMAX(i, l);
          int jl = IMIN(j, l);
          int lj = IMAX(j, l);
          oldi = a[TINDEX(li, il, n)];
          oldj = a[TINDEX(lj, jl, n)];
          a[TINDEX(li, il, n)] = cost * oldi - sint * oldj;
          a[TINDEX(lj, jl, n)] = sint * oldi + cost * oldj;
        }
        a[TINDEX(i, i, n)] =
            SQUARE(cost) * q + SQUARE(sint) * r - 2 * cost * sint * p;
        a[TINDEX(j, j, n)] =
            SQUARE(sint) * q + SQUARE(cost) * r + 2 * cost * sint * p;
        a[TINDEX(i, j, n)] =
            cost * sint * (q - r) + (SQUARE(cost) - SQUARE(sint)) * p;
      }
    }
    fnew = 0.0;
    for (int i = 1; i <= n; i++) {
      fnew += SQUARE(a[TINDEX(i, i, n)]);
    }
    if (*verbose == true) {
#ifdef USING_R
      Rprintf("itel %4d fold %15.10f fnew %15.10f\n", *itel, fold, fnew);
#else
      printf("itel %4d fold %15.10f fnew %15.10f\n", *itel, fold, fnew);
#endif
    }
    if (((fnew - fold) < *eps) || (*itel == *itmax))
      break;
    fold = fnew;
    (*itel)++;
  }
  *f = fnew;
  return;
}

void smpinverse(int *nn, double *a, int *type, int *itmax, double *eps) {
  double fini = 0.0, f = 0.0, s = 0.0, dg = 0.0;
  bool verbose = false, vectors = true;
  int n = *nn, itel = 1;
#ifdef USING_R
  double *evec = (double *)Calloc((size_t)n * n, double);
  double *diag = (double *)Calloc((size_t)n, double);
#else
  double *evec = (double *)calloc((size_t)n * n, sizeof(double));
  double *diag = (double *)calloc((size_t)n, sizeof(double));
#endif
  (void)sjacobi(nn, a, evec, &fini, &f, &itel, itmax, eps, &verbose, &vectors);
  for (int i = 1; i <= n; i++) {
    dg = a[TINDEX(i, i, n)];
    if (*type == 0) {
      diag[VINDEX(i)] = (fabs(dg) < *eps) ? 0.0 : 1.0 / dg;
    }
    if (*type == 1) {
      diag[VINDEX(i)] = (dg < *eps) ? 0.0 : sqrt(dg);
    }
    if (*type == 2) {
      diag[VINDEX(i)] = (dg < *eps) ? 0.0 : 1.0 / sqrt(dg);
    }
  }
  for (int j = 1; j <= n; j++) {
    for (int i = j; i <= n; i++) {
      s = 0.0;
      for (int k = 1; k <= n; k++) {
        s += evec[MINDEX(i, k, n)] * evec[MINDEX(j, k, n)] * diag[VINDEX(k)];
      }
      a[TINDEX(i, j, n)] = s;
    }
  }
#ifdef USING_R
  (void)Free(evec);
  (void)Free(diag);
#else
  (void)free(evec);
  (void)free(diag);
#endif
  return;
}
