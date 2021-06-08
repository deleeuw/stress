#include <stdio.h>

double bns(const int, const int, const int, const double, const double *);
double bms(const int, const int, const int, const double, const double *);
double ims(const int, const int, const int, const double, const double *);
double jms(const int, const int, const int, const double, const double *);

void splinebasis(const int *, const int *, const int *, const double *,
                 const double *, const int *, double *);

inline int MINDEX(const int, const int, const int);
inline int VINDEX(const int);

inline int MINDEX(const int i, const int j, const int nrow) {
    return (j - 1) * nrow + (i - 1);
}

inline int VINDEX(const int i) { return i - 1; }

void splinebasis(const int *d, const int *n, const int *m, const double *x,
                 const double *knots, const int *type, double *basis) {
    int mm = *m, dd = *d, nn = *n;
    int k = mm - dd - 1, i, j;
    for (i = 1; i <= nn; i++) {
        if (x[VINDEX(i)] == knots[VINDEX(mm)]) {
            basis[MINDEX(i, k, nn)] = 1.0;
            for (j = 1; j <= (k - 1); j++) {
                basis[MINDEX(i, j, nn)] = 0.0;
            }
        } else {
            for (j = 1; j <= k; j++) {
                if (*type == 0) {
                    basis[MINDEX(i, j, nn)] =
                        bms(mm, j, dd + 1, x[VINDEX(i)], knots);
                }
                if (*type == 1) {
                    basis[MINDEX(i, j, nn)] =
                        bns(mm, j, dd + 1, x[VINDEX(i)], knots);
                }
                if (*type == 2) {
                    basis[MINDEX(i, j, nn)] =
                        jms(mm, j, dd + 1, x[VINDEX(i)], knots);
                }
            }
        }
    }
}

double bns(const int nknots, const int nspline, const int updegree,
           const double x, const double * knots) {
    double y, y1, y2, temp1, temp2, k0 = knots[VINDEX(nspline + 1)],
                                    k1 = knots[VINDEX(nspline)],
                                    k2 = knots[VINDEX(nspline + updegree)],
                                    k3 = knots[VINDEX(nspline + updegree - 1)];
    if (updegree == 1) {
        if ((x >= k1) && (x < k0))
            y = 1.0;
        else
            y = 0.0;
    } else {
        temp1 = 0.0;
        if ((k3 - k1) > 0) temp1 = (x - k1) / (k3 - k1);
        temp2 = 0.0;
        if ((k2 - k0) > 0) temp2 = (k2 - x) / (k2 - k0);
        y1 = bns(nknots, nspline, updegree - 1, x, knots);
        y2 = bns(nknots, nspline + 1, updegree - 1, x, knots);
        y = temp1 * y1 + temp2 * y2;
    }
    return y;
}

double bms(const int nknots, const int nspline, const int updegree,
           const double x, const double * knots) {
    double y, y1, y2, temp1, temp2,
           k0 = knots[VINDEX(nspline + 1)], k1 = knots[VINDEX(nspline)],
           k2 = knots[VINDEX(nspline + updegree)], kd = 1.0 / (k0 - k1);
    if (updegree == 1) {
        if ((x >= k1) && (x < k0))
            y = kd;
        else
            y = 0.0;
    } else {
        temp1 = 0.0;
        temp2 = 0.0;
        if ((k2 - k1) > 0) {
            temp1 = (x - k1) / (k2 - k1);
            temp2 = (k2 - x) / (k2 - k1);
        }
        y1 = bms(nknots, nspline, updegree - 1, x, knots);
        y2 = bms(nknots, nspline + 1, updegree - 1, x, knots);
        y = temp1 * y1 + temp2 * y2;
    }
    return y;
}


double jms(const int nknots, const int nspline, const int updegree,
           const double x, const double * knots) {
    double sum = 0.0;
    if (x >= knots[VINDEX(nspline + updegree)]) return 1.0;
    for (int j = nspline; j < nspline + updegree; j++) {
        sum += bns(nknots, j, updegree + 1, x, knots);
    }
    return sum;
}

double cox (const int nknots, const int nspline, const int order,
            const double x, const double * knots, double * v) {
    int k;
    for (int j = 2; j <= nknots; j++) {
        if (x < knots[VINDEX(j)]) && (x >= knots[VINDEX(j - 1)]) k = j;
    }
    printf ("%10.4f %10.4f %10.4f\n", knots[VINDEX(j - 1)], x, knots[VINDEX(j)]);
}
