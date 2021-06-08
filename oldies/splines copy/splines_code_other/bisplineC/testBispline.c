#include <stdio.h>
#include <stdlib.h>
#include "bispline.h"

double x[] = {0.4, 1.1, 1.2, 2.5, 3.3, 5.5};
double knots[] = {0, 0, 0, 1, 2, 3, 4, 5, 6, 6, 6};
int degree = 2;

extern void bisplineC(const int *, const double *, const int *, const int *,
                      const double *, int *, double *, double *, int *);
extern void bisplineBasisC(const int *, const double *, const int *,
                           const int *, const double *, const int *, double *,
                           double *);

extern inline int INDEX(const int);
extern inline int IMIN(const int, const int);
extern inline int IMAX(const int, const int);
extern inline int MINDEX(const int, const int, const int);

void bisplineC(const int *nknots, const double *knots, const int *order,
               const int *norm, const double *xvalue, int *jint, double *v,
               double *vint, int *ifail) {
    int n = *nknots, k = *order, rm = *norm;
    double x = *xvalue, xleft = knots[INDEX(1)], xrite = knots[INDEX(n)];
    if (n < (k + 1)) {
        *ifail = 1;
        return;
    }
    if (n < 1) {
        *ifail = 2;
        return;
    }
    if ((rm != 1) && (rm != 2)) {
        *ifail = 5;
        return;
    }
    for (int i = 1; i <= n; i++) {
        v[INDEX(i)] = 0.0;
        vint[INDEX(i)] = 0.0;
    }
    if (x < xleft) {
        *jint = 1;
        return;
    }
    if (x >= xrite) {
        int i1 = n, i2 = n - 1;
        for (int i = 1; i <= i2; i++) {
            int irmi = n - i;
            if (knots[INDEX(irmi)] >= x) {
                i1 = irmi - 1;
            }
        }
        *jint = i1;
        if (*jint <= 0) {
            *ifail = 4;
            return;
        }
    }
    int il, ir;
    if (x < xrite) {
        il = 1;
        ir = n;
        while ((ir - il) > 1) {
            int middle = (int)(ir + il) / 2;
            if (x < knots[INDEX(middle)]) {
                ir = middle;
            } else {
                il = middle;
            }
        }
        *jint = il;
    }
    int jintmk = *jint - k;
    ir = IMIN(n, *jint + k) - 1;
    if (*jint == n) {
        il = IMIN(ir, IMAX(1, jintmk));
    }
    if (*jint < n) {
        il = IMIN(ir, IMAX(1, jintmk + 1));
    }
    for (int i = il; i <= ir; i++) {
        if (knots[INDEX(i)] > knots[INDEX(i + 1)]) {
            *ifail = 3;
            return;
        }
    }
    int i = il, irp1 = ir + 1;
    while (TRUE) {
        int ipk = i + k;
        if (ipk > irp1) {
            break;
        }
        if (knots[INDEX(i)] >= knots[INDEX(ipk)]) {
            *ifail = 4;
            return;
        }
        i++;
    }
    int jj1 = n - *jint;
    double dk = (double)k;
    if (jj1 <= 0) {
        if (rm == 1) {
            vint[INDEX(1)] = 1.0 / dk;
        }
        if (rm == 2) {
            vint[INDEX(1)] = (knots[INDEX(*jint)] - knots[INDEX(jintmk)]) / dk;
        }
        return;
    }
    double e1 = x - knots[INDEX(*jint)], e2 = knots[INDEX(*jint + 1)] - x;
    v[INDEX(1)] = 1.0 / (knots[INDEX(*jint + 1)] - knots[INDEX(*jint)]);
    vint[INDEX(1)] = e1 * v[INDEX(1)];
    if (k == 1) {
        if (rm == 1) {
            return;
        }
        v[INDEX(1)] = 1.0;
        vint[INDEX(1)] = e1;
        return;
    }
    int nk = IMIN(k, jj1);
    if (nk != 1) {
        for (int j = 2; j <= nk; j++) {
            int jintpj = *jint + j;
            if (!((j == k) && (rm == 2))) {
                v[INDEX(j)] = e1 * v[INDEX(j - 1)] /
                              (knots[INDEX(jintpj)] - knots[INDEX(*jint)]);
                vint[INDEX(j)] = e1 * v[INDEX(j)];
            } else {
                v[INDEX(j)] = e1 * v[INDEX(j - 1)];
                vint[INDEX(j)] = e1 * v[INDEX(j)] /
                                 (knots[INDEX(jintpj)] - knots[INDEX(*jint)]);
            }
        }
    }
    if (*jint != 1) {
        int mj = IMIN(*jint, k) - 1;
        for (int j = 1; j <= mj; j++) {
            int jintmj = *jint - j;
            double e3 = x - knots[INDEX(jintmj)];
            if (!(((j + 1) == k) && (rm == 2))) {
                v[INDEX(1)] = e2 * v[INDEX(1)] /
                              (knots[INDEX(*jint + 1)] - knots[INDEX(jintmj)]);
                vint[INDEX(1)] = vint[INDEX(1)] + e3 * v[INDEX(1)];
            } else {
                v[INDEX(1)] = e2 * v[INDEX(1)];
                vint[INDEX(1)] = vint[INDEX(1)] +
                                 e3 * v[INDEX(1)] / (knots[INDEX(*jint + 1)] -
                                                     knots[INDEX(jintmj)]);
            }
            int nkj = IMIN(k - j, jj1);
            if (nkj <= 1) {
                continue;
            }
            for (int l = 2; l <= nkj; l++) {
                int jintpl = *jint + l;
                if (!(((j + l) == k) && (rm == 2))) {
                    v[INDEX(l)] = (e3 * v[INDEX(l - 1)] +
                                   (knots[INDEX(jintpl)] - x) * v[INDEX(l)]) /
                                  (knots[INDEX(jintpl)] - knots[INDEX(jintmj)]);
                    vint[INDEX(l)] = vint[INDEX(l)] + e3 * v[INDEX(l)];
                } else {
                    v[INDEX(l)] = e3 * v[INDEX(l - 1)] +
                                  (knots[INDEX(jintpl)] - x) * v[INDEX(l)];
                    vint[INDEX(l)] = vint[INDEX(l)] +
                                     e3 * v[INDEX(l)] / (knots[INDEX(jintpl)] -
                                                         knots[INDEX(jintmj)]);
                }
            }
        }
    }
    for (int j = 1; j <= k; j++) {
        int jintpj = *jint + j;
        int jtpjmk = jintpj - k;
        if (jintpj > n) {
            return;
        }
        if (jtpjmk >= 1) {
            if (rm == 1) {
                vint[INDEX(j)] = vint[INDEX(j)] / dk;
            }
            if (rm == 2) {
                vint[INDEX(j)] = (knots[INDEX(jintpj)] - knots[INDEX(jtpjmk)]) *
                                 vint[INDEX(j)] / dk;
            }
        } else {
            v[INDEX(j)] = 0.0;
            vint[INDEX(j)] = 0.0;
        }
    }
    return;
}

void bisplineBasisC(const int *nknots, const double *knots, const int *order,
                    const int *nvalues, const double *x, const int *norm,
                    double *bbasis, double *ibasis) {
    int n = *nknots, m = *nvalues, k = *order;
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            bbasis[MINDEX(i, j, m)] = 0.0;
        }
    }
    for (int i = 1; i <= m; i++) {
        int jint = 0, ifail = 0;
        double *v = (double *)calloc(n, sizeof(double));
        double *vint = (double *)calloc(n, sizeof(double));
        double xvalue = x[INDEX(i)];
        (void)bisplineC(nknots, knots, order, norm, &xvalue, &jint, v, vint,
                        &ifail);
        for (int j = 1; j <= k; j++) {
            int j1 = jint - k + j;
            bbasis[MINDEX(i, j1, m)] = v[INDEX(j)];
        }
        free(v);
        free(vint);
    }
}

int main(void) {
    int nknots = 11, order = degree + 1, nvalues = 6, norm = 2;

    double *bbasis = (double *)calloc(nknots * nvalues, sizeof(double));
    double *ibasis = (double *)calloc(nknots * nvalues, sizeof(double));

    (void)bisplineBasisC(&nknots, knots, &order, &nvalues, x, &norm, bbasis,
                         ibasis);

    for (int i = 1; i <= 6; i++) {
        for (int j = 1; j <= 11; j++) {
            printf("%10.6f", bbasis[MINDEX(i, j, 6)]);
        }
        printf("\n");
    }

    free(bbasis);
    free(ibasis);
}