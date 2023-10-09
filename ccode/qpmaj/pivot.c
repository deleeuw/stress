
#include <math.h>
#include <stdio.h>

#define MINDEX(i, j, n) (((j)-1) * (n) + (i)-1)

void pivotOneC (double *, int *, int *, int *, int *, int *, double *eps);
void pivotC (double *, int *, int *, int *, int *, int *, int *, int *, int *, double *eps);
void mprint (double *, int *, int *);

void pivotOneC (double *a, int *kk, int *nn, int *mm, int *type, int *refuse, double *eps) {
    int i, j, ik, kj, ij, kp, s00, s01, s10, n = *nn, m = *mm, k = *kk, tp = *type;
    double pv, ee = *eps;
    *refuse = 0;
    if (tp == 1) {
        s00 = -1;
        s01 = 1;
        s10 = 1;
    }
    if (tp == 2) {
        s00 = -1;
        s01 = -1;
        s10 = -1;
    }
   if (tp == 3) {
        s00 = 1;
        s01 = -1;
        s10 = 1;
    }
   if (tp == 4) {
        s00 = 1;
        s01 = 1;
        s10 = -1;
    }
    kp = MINDEX(k, k, n);
    pv = a[kp];
    if (fabs (pv) < ee) {
        *refuse = 1; 
        return;
    }
    for (i = 1; i <= n; i++) {
        if (i == k) continue;
        ik = MINDEX (i, k , n);
        for (j = 1; j <= m; j++) {
            if (j == k) continue;
            kj = MINDEX (k, j, n);
            ij = MINDEX (i, j, n);
            a[ij] = a[ij] - a[ik] * a[kj] / pv;
        }
    }
    for (i = 1; i <= n; i++) {
        if (i == k) continue;
        ik = MINDEX (i, k, n);
        a[ik] = s10 * a[ik] / pv;
    }
    for (j = 1; j <= m; j++) {
        if (j == k) continue;
        kj = MINDEX (k, j, n);
        a[kj] = s01 * a[kj] / pv;
    }
    a[kp] = s00 / pv;
}


void pivotC (double *a, int *ind, int *jnd, int *done, int *nn, int *mm, int *pp, int *type, int *skip, double *eps) {
    int ipiv, kpiv, lpiv, refuse;
    int n = *nn, m = *mm, p = *pp;
    double fmax, fpiv;
    for (int l = 0; l < p; l++) {
        jnd[l] = 0;
        done[l] = 0;
        skip[l] = 0;
    }
    for (int l = 1; l <= p; l++) {
        fmax = -1.0;
        for (int k = 1; k <= p; k++) {
            kpiv = ind[k - 1];
            if (done[k - 1] == 1) continue;
            fpiv = fabs (a[MINDEX (kpiv, kpiv, n)]);
            if (fpiv > fmax) {
                ipiv = k;
                lpiv = kpiv;
                fmax = fpiv;
            }
        }
        done[ipiv - 1] = 1;
        jnd [l - 1] = lpiv;
        pivotOneC(a, &lpiv, &n, &m, type, &refuse, eps);
        skip[l - 1] = refuse;
    }
}


void mprint (double * a, int *nn, int *mm) {
    int i, j, ij, n = *nn, m = *mm;
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= m; j++) {
            ij = MINDEX (i, j, n);
            printf (" %5.3f ", a[ij]);
        }
    printf ("\n");
    }
}