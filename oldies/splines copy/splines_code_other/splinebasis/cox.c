#include <stdio.h>

void cox (const int *, const int *, const int *, const int *, const double *, const double *, double *, double *);
int jintf (const double, const int, const double *);

inline int MINDEX(const int, const int, const int);
inline int VINDEX(const int);

inline int MINDEX(const int i, const int j, const int nrow) {
    return (j - 1) * nrow + (i - 1);
}

inline int VINDEX(const int i) { return i - 1; }

void cox (const int *nknots, const int *i, const int *order, const int *norm, const double *x, const double * knots, double *v, double *vint) {
    int jint = jintf (*x, *nknots, knots), k = *order, n = *nknots, jintmk = jint - k, jj1 = n - jint;
    double dk = (double) *order, xvalue = *x;
    if (jint <= n) {
        vint[VINDEX(1)] = (*norm == 1) ? 1.0 / dk : (knots[VINDEX(jint)] - knots[VINDEX(jintmk)]) / dk;
    }
    double e1 = xvalue - knots[VINDEX(jint)], e2 = knots[VINDEX(jint + 1)] - xvalue;
    v[VINDEX(1)] = 1.0 / (knots[VINDEX(jint + 1)] - knots[VINDEX(jint)]);
    vint[VINDEX(1)] = e1 * v[VINDEX(1)];
    if (k == 1) {
        if (*norm == 1) return;
        v[VINDEX(1)] = 1.0;
        vint[VINDEX(1)] = e1;
        return;
    }
}

int jintf (const double x, const int nknots, const double *knots) {
    double xleft = knots[VINDEX(1)], xrite = knots[VINDEX(nknots)];
    if (x < xleft) {
        return 1;
    }
    if (x > xrite) {
        return nknots;
    }
    if (x < xrite) {
        for (int j = 1; j < nknots; j++) {
            if ((x < knots[VINDEX(j + 1)]) && (x >= knots[VINDEX(j)])) {
                return j;
            }
        }
    }
    else {
        if (knots[VINDEX(nknots - 1)] < xrite) {
            return nknots;
        }
        else {
            for (int j = nknots; j >= 1; j--) {
                if (knots[VINDEX(j)] < xrite) {
                    return j;
                }
            }

        }
    }
    return 0;
}
