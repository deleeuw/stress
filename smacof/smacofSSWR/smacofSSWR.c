#include "../include/smacof.h"

int main() {
  double delta[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double w[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  double x[8] = {1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
  int n = 4, p = 2, itmax = 100;
  bool verbose = true, speedup = false;
  double eps = 1e-10;
  (void)smacofSSWR(delta, w, x, &n, &p, &speedup, &itmax, &eps, &verbose);
  printf("\n\n");
  (void)primat(4, 2, 15, 10, x);
  double *dist = (double *)calloc((size_t) 6, sizeof(double));
  (void)smacofDist(x, n, p, dist);
  (void)smacofGradientW(dist, w, delta, n, p, x); 
  (void)primat(4, 2, 15, 10, x);
  free(dist);
  return EXIT_SUCCESS;
}

// Single Matrix, Symmetric, Weighted Stress, Ratio Transform

void smacofSSWR(double *delta, const double *w, double *x, const int *nobjects,
                const int *ndimensions, const bool *speedup, const int *itmax,
                const double *eps, const bool *verbose) {
  int n = *nobjects, p = *ndimensions, acc = *speedup;
  int m = n * (n - 1) / 2, np = n * p, itel = 1;
  double sold = 0.0, snew = 0.0;
  // normalize delta
  (void)smacofNormDeltaW(w, delta, m);
  // Compute the MP inverse of V and test for irreducibility
  double *vinv = (double *)calloc((size_t)m, sizeof(double));
  (void)smacofInvertVW(w, vinv, n, m);
  // Compute initial distances
  double *dist = (double *)calloc((size_t)m, sizeof(double));
  (void)smacofDist(x, n, p, dist);
  // Scale initial configuration and distances
  (void)smacofScaleXW(x, dist, w, delta, m, np);
  // Compute initial stress
  sold = smacofLossW(dist, w, delta, m);
  while (true) {
    (void)smacofGuttmanW(dist, w, delta, vinv, n, p, acc, x);
    (void)smacofDist(x, n, p, dist);
    snew = smacofLossW(dist, w, delta, m);
    if (verbose) {
      printf("itel = %4d sold = %15.10f snew = %15.10f\n", itel, sold, snew);
    }
    if ((itel == *itmax) || ((sold - snew) < *eps)) {
      break;
    }
    itel += 1;
    sold = snew;
  }
  free(vinv);
  free(dist);
  return;
}

void smacofNormDeltaW(const double *w, double *delta, const int m) {
  double s = 0.0;
  for (int k = 0; k < m; k++) {
    s += w[k] * SQUARE(delta[k]);
  }
  double r = sqrt(((double) m) / s);
  for (int k = 0; k < m; k++) {
    delta[k] *= r;
  }
  return;
}

void smacofInvertVW(const double *w, double *vinv, const int n, const int m) {
  double ni = 1.0 / (double)n, s = 0.0;
  for (int k = 0; k < m; k++) {
    vinv[k] = -w[k] + ni;
  }
  double *vsum = (double *)calloc((size_t)n, sizeof(double));
  for (int i = 1; i <= n; i++) {
    s = 0.0;
    for (int j = i + 1; j <= n; j++) {
      s += w[SINDEX(j, i, n)];
    }
    for (int j = 1; j < i; j++) {
      s += w[SINDEX(i, j, n)];
    }
    vsum[VINDEX(i)] = s + ni;
  }
  double *vkrw = (double *)calloc((size_t)n, sizeof(double));
  for (int k = 1; k <= n; k++) {
    for (int j = 1; j <= n; j++) {
      if (j == k) {
        vkrw[VINDEX(k)] = 0.0;
      }
      if (j < k) {
        vkrw[VINDEX(j)] = vinv[SINDEX(k, j, n)];
      }
      if (j > k) {
        vkrw[VINDEX(j)] = vinv[SINDEX(j, k, n)];
      }
    }
    s = vsum[VINDEX(k)];
    if (s < 1e-15) {
      printf("Error: W is not irreducible.\n");
      exit(EXIT_FAILURE);
    }
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= i; j++) {
        if (i == j) {
          vsum[VINDEX(i)] -= SQUARE(vkrw[VINDEX(i)]) / s;
        } else {
          vinv[SINDEX(i, j, n)] -= vkrw[VINDEX(i)] * vkrw[VINDEX(j)] / s;
        }
      }
    }
    for (int j = 1; j < k; j++) {
      vinv[SINDEX(k, j, n)] /= s;
    }
    for (int i = k + 1; i <= n; i++) {
      vinv[SINDEX(i, k, n)] /= s;
    }
    vsum[VINDEX(k)] = -1.0 / s;
  }
  for (int k = 0; k < m; k++) {
    vinv[k] = -vinv[k] - ni;
  }
  free(vsum);
  free(vkrw);
  return;
}

void smacofScaleXW(double *x, double *dist, const double *w,
                   const double *delta, const int m, const int np) {
  double sd1 = 0.0, sd2 = 0.0;
  for (int k = 0; k < m; k++) {
    sd1 += w[k] * dist[k] * delta[k];
    sd2 += w[k] * SQUARE(dist[k]);
  }
  double lbd = sd1 / sd2;
  for (int k = 0; k < m; k++) {
    dist[k] *= lbd;
  }
  for (int k = 0; k < np; k++) {
    x[k] *= lbd;
  }
}

double smacofLossW(const double *dist, const double *w, const double *delta,
                   const int m) {
  double stress = 0.0;
  for (int k = 0; k < m; k++) {
    stress += w[k] * SQUARE(delta[k] - dist[k]);
  }
  return stress;
}

void smacofGuttmanW(const double *dist, const double *w, const double *delta,
                    const double *vinv, const int n, const int p,
                    const bool speedup, double *x) {
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
        y[VINDEX(i)] += w[k] * delta[k] *
                        (x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]) / dist[k];
      }
    }
    double *z = calloc((size_t)n, sizeof(double));
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
        z[VINDEX(i)] += vinv[k] * (y[VINDEX(j)] - y[VINDEX(i)]);
      }
      if (speedup) {
        z[VINDEX(i)] = 2 * z[VINDEX(i)] - x[MINDEX(i, s, n)];
      }
    }
    for (int i = 1; i <= n; i++) {
      x[MINDEX(i, s, n)] = z[VINDEX(i)];
    }
    free(y);
    free(z);
  }
}

void smacofGradientW(double *dist, double *w, double *delta, const int n,
                       const int p, double *x) { 
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
        y[VINDEX(i)] += w[k] * (1 - delta[k] / dist[k]) *
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

