
#include "smacof.h"

void smacofDistC(const double *x, const int *n, const int *p, double *d) {
    int k = 1, nn = *n, pp = *p;
    for (int j = 1; j <= nn - 1; j++) {
        for (int i = j + 1; i <= nn; i++) {
            double dij = 0.0;
            for (int s = 1; s <= pp; s++) {
                dij += SQUARE(x[MINDEX(i, s, nn)] - x[MINDEX(j, s, nn)]);
            }
            d[VINDEX(k)] = sqrt(dij);
            k++;
        }
    }
}

void smacofLossC(const double *dist, const double *w, const double *delta,
                 const int *m, double *stress) {
    int mm = *m;
    *stress = 0.0;
    for (int k = 1; k <= mm; k++) {
        *stress += w[VINDEX(k)] * SQUARE(delta[VINDEX(k)] - dist[VINDEX(k)]);
    }
    *stress /= 2.0;
    return;
}

void smacofBmatC(const double *dist, const double *w, const double *delta,
                 const int *n, double *b) {
    int nn = *n;
    for (int i = 1; i <= nn; i++) {
        b[TINDEX(i, i, nn)] = 0.0;
    }
    for (int j = 1; j <= nn - 1; j++) {
        for (int i = j + 1; i <= nn; i++) {
            int k = SINDEX(i, j, nn);
            double dinv = (dist[k] == 0.0) ? 0.0 : 1.0 / dist[k];
            double elem = w[k] * delta[k] * dinv;
            b[TINDEX(i, j, nn)] = -elem;
            b[TINDEX(i, i, nn)] += elem;
            b[TINDEX(j, j, nn)] += elem;
        }
    }
    return;
}

void smacofVmatC(const double *w, const int *n, double *v) {
    int nn = *n;
    for (int i = 1; i <= nn; i++) {
        v[TINDEX(i, i, nn)] = 0.0;
    }
    for (int j = 1; j <= nn - 1; j++) {
        for (int i = j + 1; i <= nn; i++) {
            double elem = w[SINDEX(i, j, nn)];
            v[TINDEX(i, j, nn)] = -elem;
            v[TINDEX(i, i, nn)] += elem;
            v[TINDEX(j, j, nn)] += elem;
        }
    }
    return;
}

void smacofHmatC(const double *x, const double *bmat, const double *vmat,
                 const double *w, const double *delta, const double *dist,
                 const int *n, const int *p, double *work, double *h) {
    int nn = *n, pp = *p, np = nn * pp;
    for (int s = 1; s <= pp; s++) {
        for (int t = 1; t <= s; t++) {
            for (int j = 1; j <= nn; j++) {
                for (int i = j; i <= nn; i++) {
                    work[TINDEX(i, j, nn)] = 0.0;
                }
            }
            for (int j = 1; j <= nn - 1; j++) {
                for (int i = j + 1; i <= nn; i++) {
                    double f1 = (x[MINDEX(i, s, nn)] - x[MINDEX(j, s, nn)]);
                    double f2 = (x[MINDEX(i, t, nn)] - x[MINDEX(j, t, nn)]);
                    double f3 = THIRD(dist[SINDEX(i, j, nn)]);
                    double f4 = delta[SINDEX(i, j, nn)] * w[SINDEX(i, j, nn)];
                    f3 = (f3 < 1e-10) ? 0 : 1 / f3;
                    double elem = f1 * f2 * f4 * f3;
                    work[TINDEX(i, j, nn)] = -elem;
                    work[TINDEX(i, i, nn)] += elem;
                    work[TINDEX(j, j, nn)] += elem;
                }
            }
            if (s == t) {
                for (int j = 1; j <= nn; j++) {
                    for (int i = j; i <= nn; i++) {
                        h[TINDEX((s - 1) * nn + i, (s - 1) * nn + j, nn * pp)] =
                            bmat[TINDEX(i, j, nn)] - work[TINDEX(i, j, nn)];
                    }
                }
            }
            if (s != t) {
                for (int i = 1; i <= nn; i++) {
                    for (int j = 1; j <= nn; j++) {
                        int ij = IMIN(i, j);
                        int ji = IMAX(i, j);
                        int sj = (s - 1) * nn + j;
                        int si = (t - 1) * nn + i;
                        h[TINDEX(sj, si, np)] = -work[TINDEX(ji, ij, nn)];
                    }
                }
            }
        }
    }
    return;
}

void smacofGuttmanC(const double *x, const double *bmat, const double *vinv,
                    const int *n, const int *p, double *work, double *y) {
    (void)mutrma(n, p, bmat, x, work);
    (void)mutrma(n, p, vinv, work, y);
    return;
}

void smacofGradientC(const double *x, const double *bmat, const double *vmat,
                     const int *n, const int *p, double *work, double *y) {
    int nn = *n, pp = *p;
    (void)mutrma(n, p, vmat, x, y);
    (void)mutrma(n, p, bmat, x, work);
    for (int i = 1; i <= nn; i++) {
        for (int s = 1; s <= pp; s++) {
            y[MINDEX(i, s, nn)] -= work[MINDEX(i, s, nn)];
        }
    }
    return;
}

void smacofHessianC(const double *x, const double *bmat, const double *vmat,
                    const double *w, const double *delta, const double *dist,
                    const int *n, const int *p, double *work, double *h) {
    int nn = *n, pp = *p, np = nn * pp;
    (void)smacofHmatC(x, bmat, vmat, w, delta, dist, n, p, work, h);
    for (int j = 1; j <= np; j++) {
        for (int i = j; i <= np; i++) {
            h[TINDEX(i, j, np)] = -h[TINDEX(i, j, np)];
        }
    }
    for (int s = 1; s <= pp; s++) {
        for (int j = 1; j <= nn; j++) {
            for (int i = j; i <= nn; i++) {
                h[TINDEX((s - 1) * nn + i, (s - 1) * nn + j, np)] +=
                    vmat[TINDEX(i, j, nn)];
            }
        }
    }
    return;
}

void smacofNewtonC() { return; }

void smacofInitialC(const double *delta, const int *n, const int *p,
                    double *work1, double *work2, double *work3, double *work4,
                    double *x) {
    int nn = *n, pp = *p, itmax = 100;
    double s, ss = 0.0, eps = 1e-6;
    for (int i = 1; i <= nn; i++) {
        work1[VINDEX(i)] = 0.0;
        for (int j = 1; j <= nn; j++) {
            if (i == j) continue;
            int ij = IMAX(i, j);
            int ji = IMIN(i, j);
            s = SQUARE(delta[SINDEX(ij, ji, nn)]);
            ss += s;
            work1[VINDEX(i)] += s;
            if (j < i) continue;
            work2[TINDEX(ij, ji, nn)] = s;
        }
        work1[VINDEX(i)] /= (double)nn;
    }
    ss /= SQUARE((double)nn);
    for (int j = 1; j <= nn; j++) {
        for (int i = j; i <= nn; i++) {
            work2[TINDEX(i, j, nn)] -= work1[VINDEX(i)];
            work2[TINDEX(i, j, nn)] -= work1[VINDEX(j)];
            work2[TINDEX(i, j, nn)] += ss;
            work2[TINDEX(i, j, nn)] *= -0.5;
        }
    }
    (void)jacobiC(n, work2, work3, work1, work4, &itmax, &eps);
    for (int i = 1; i <= nn; i++) {
        for (int j = 1; j <= pp; j++) {
            s = work2[TINDEX(j, j, nn)];
            if (s <= 0) continue;
            x[MINDEX(i, j, nn)] = work3[MINDEX(i, j, nn)] * sqrt(s);
        }
    }
    return;
}

void smacofUpdateC(const double *xold, const double *w, const double *delta,
                   const double *vinv, const int *n, const int *p, double *dnew,
                   double *bmat, double *work, double *snew, double *xnew) {
    int nn = *n, pp = *p, mm = nn * (nn - 1) / 2;
    (void)mutrma(n, p, bmat, xold, work);
    (void)mutrma(n, p, vinv, work, xnew);
    for (int j = 1; j <= nn - 1; j++) {
        for (int i = j + 1; i <= nn; i++) {
            double dij = 0.0;
            for (int s = 1; s <= pp; s++) {
                dij += SQUARE(xnew[MINDEX(i, s, nn)] - xnew[MINDEX(j, s, nn)]);
            }
            dnew[SINDEX(i, j, nn)] = sqrt(dij);
        }
    }
    for (int i = 1; i <= nn; i++) {
        bmat[TINDEX(i, i, nn)] = 0.0;
    }
    for (int j = 1; j <= nn - 1; j++) {
        for (int i = j + 1; i <= nn; i++) {
            int k = SINDEX(i, j, nn);
            double dinv = (dnew[k] == 0.0) ? 0.0 : 1.0 / dnew[k];
            double elem = w[k] * delta[k] * dinv;
            bmat[TINDEX(i, j, nn)] = -elem;
            bmat[TINDEX(i, i, nn)] += elem;
            bmat[TINDEX(j, j, nn)] += elem;
        }
    }
    *snew = 0.0;
    for (int k = 1; k <= mm; k++) {
        *snew += w[VINDEX(k)] * SQUARE(delta[VINDEX(k)] - dnew[VINDEX(k)]);
    }
    *snew /= 2.0;
}
