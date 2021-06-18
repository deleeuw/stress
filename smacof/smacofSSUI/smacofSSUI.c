#include "../include/smacof.h"

int main() {
  double delta[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double x[8] = {1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
  int n = 4, p = 2, itmax = 100;
  bool verbose = true, speedup = false;
  double eps = 1e-15;
  (void)smacofSSUI(delta, x, &n, &p, &speedup, &itmax, &eps, &verbose);
  printf("\n\n");
  (void)primat(4, 2, 15, 10, x);
  double *dist = (double *)calloc((size_t)6, sizeof(double));
  (void)smacofDist(x, n, p, dist);
  (void)primat(1, 6, 15, 10, delta);
  (void)primat(1, 6, 15, 10, dist);
  (void)smacofGradientU(dist, delta, n, p, x);
  (void)primat(4, 2, 15, 10, x);
  free(dist);
  return EXIT_SUCCESS;
}

void smacofSSUI(double *delta, double *x, const int *nobjects,
                const int *ndimensions, const bool *speedup, const int *itmax,
                const double *eps, const bool *verbose) {
  int n = *nobjects, p = *ndimensions, acc = *speedup;
  int m = n * (n - 1) / 2, np = n * p, itel = 1;
  double sold = 0.0, smid = 0.0, snew = 0.0;
  // normalize delta
  (void)smacofNormDeltaU(delta, m);
  // compute initial distances
  double *dist = (double *)calloc((size_t)m, sizeof(double));
  (void)smacofDist(x, n, p, dist);
  // scale initial configuration
  (void)smacofScaleXU(x, dist, delta, m, np);
  // compute initial stress
  sold = smacofLossU(dist, delta, m);
  while (true > false) {
    (void)smacofGuttmanU(dist, delta, n, p, acc, x);
    (void)smacofDist(x, n, p, dist);
    smid = smacofLossU(dist, delta, m);
    (void)smacofIntervalU(dist, delta, m);
    snew = smacofLossU(dist, delta, m);
    if (verbose) {
      printf("itel = %4d sold = %15.10f smid = %15.10f snew = %15.10f\n", itel,
             sold, smid, snew);
    }
    if ((itel == *itmax) || ((sold - snew) < *eps)) {
      break;
    }
    itel += 1;
    sold = snew;
  }
  free(dist);
  return;
}

void smacofNormDeltaU(double *delta, const int m) {
  double s = 0.0, r = 0.0;
  for (int k = 0; k < m; k++) {
    s += SQUARE(delta[k]);
  }
  r = sqrt(((double)m) / s);
  for (int k = 0; k < m; k++) {
    delta[k] *= r;
  }
  return;
}

void smacofScaleXU(double *x, double *dist, const double *delta, const int m,
                   const int np) {
  double sd1 = 0.0, sd2 = 0.0;
  for (int k = 0; k < m; k++) {
    sd1 += dist[k] * delta[k];
    sd2 += SQUARE(dist[k]);
  }
  double lbd = sd1 / sd2;
  for (int k = 0; k < m; k++) {
    dist[k] *= lbd;
  }
  for (int k = 0; k < np; k++) {
    x[k] *= lbd;
  }
}

double smacofLossU(const double *dist, const double *delta, const int m) {
  double stress = 0.0;
  for (int k = 0; k < m; k++) {
    stress += SQUARE(delta[k] - dist[k]);
  }
  return stress;
}

void smacofGuttmanU(const double *dist, const double *delta, const int n,
                    const int p, const bool speedup, double *x) {
  int k;
  for (int s = 1; s <= p; s++) {
    double *y = calloc((size_t)n, sizeof(double));
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
        if (j == i) {
          continue;
        }
        if (i > j) {
          k = SINDEX(i, j, n);
        }
        if (j > i) {
          k = SINDEX(j, i, n);
        }
        if (dist[k] < 1e-15) {
          continue;
        }
        y[VINDEX(i)] +=
            delta[k] * (x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]) / dist[k];
      }
      if (speedup) {
        y[VINDEX(i)] = 2 * y[VINDEX(i)] / ((double)n) - x[MINDEX(i, s, n)];
      } else {
        y[VINDEX(i)] /= ((double)n);
      }
    }
    for (int i = 1; i <= n; i++) {
      x[MINDEX(i, s, n)] = y[VINDEX(i)];
    }
    free(y);
  }
}

void smacofGradientU(double *dist, double *delta, const int n, const int p,
                     double *x) {
  for (int s = 1; s <= p; s++) {
    double *y = (double *)calloc((size_t)n, sizeof(double));
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
        int k = 0;
        if (j == i) {
          continue;
        }
        if (i > j) {
          k = SINDEX(i, j, n);
        }
        if (j > i) {
          k = SINDEX(j, i, n);
        }
        if (dist[k] < 1e-15) {
          continue;
        }
        y[VINDEX(i)] += (1 - delta[k] / dist[k]) *
                        (x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]);
      }
    }
    for (int i = 1; i <= n; i++) {
      x[MINDEX(i, s, n)] = y[VINDEX(i)];
    }
    free(y);
  }
  return;
}

void smacofIntervalU(const double *dist, double *delta, int m) {
  double acum = 0.0, bcum = 0.0, ccum = 0.0, dcum = 0.0;
  double s = 0.0, deltamin = INFINITY;
  for (int k = 0; k < m; k++) {
    acum += delta[k];
    bcum += dist[k];
    deltamin = MIN(deltamin, delta[k]);
  }
  acum /= (double)m;
  bcum /= (double)m;
  for (int k = 0; k < m; k++) {
    ccum += (delta[k] - acum) * (dist[k] - bcum);
    dcum += SQUARE(delta[k] - acum);
  }
  double a0 = ccum / dcum;
  double b0 = -(a0 * acum - bcum);
  if ((a0 > 1e-15) && ((a0 * deltamin + b0) > 1e-15)) {
    s = 0.0;
    for (int k = 0; k < m; k++) {
      delta[k] = a0 * delta[k] + b0;
      s += SQUARE(delta[k]);
    }
  } else {
    ccum = 0.0;
    dcum = 0.0;
    for (int k = 0; k < m; k++) {
      ccum += (delta[k] - deltamin) * dist[k];
      dcum += SQUARE(delta[k] - deltamin);
    }
    double a1 = ccum / dcum;
    double b1 = -a1 * deltamin;
    if (a1 > 1e-15) {
      s = 0.0;
      for (int k = 0; k < m; k++) {
        delta[k] = a0 * delta[k] + b0;
        s += SQUARE(delta[k]);
      }
    } else {
      return;
    }
  }
  double r = sqrt(((double)m) / s);
  for (int k = 0; k < m; k++) {
    delta[k] *= r;
  }
  return;
}
