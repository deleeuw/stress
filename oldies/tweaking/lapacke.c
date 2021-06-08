#include "smacof.h"

void dposv(const int *n, const int *m, double *a, double *b) {
    lapack_int nn = (lapack_int)*n, mm = (lapack_int)*m;
    (void)LAPACKE_dposv(LAPACK_COL_MAJOR, 'U', nn, mm, a, nn, b, nn);
    return;
}

void dsyevd(const int *n, double *a, double *x) {
    lapack_int nn = (lapack_int)*n;
    (void)LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'V', 'U', nn, a, nn, x);
    return;
}

void dgeqrf(const int *n, const int *m, double *a, double *tau) {
    lapack_int nn = (lapack_int)*n, mm = (lapack_int)*m;
    (void)LAPACKE_dgeqrf(LAPACK_COL_MAJOR, nn, mm, a, nn, tau);
    return;
}

void dorgqr(const int *n, const int *m, double *a, double *tau) {
    lapack_int nn = (lapack_int)*n, mm = (lapack_int)*m;
    (void)LAPACKE_dorgqr(LAPACK_COL_MAJOR, nn, mm, mm, a, nn, tau);
    return;
}

void dortho(const int *n, const int *m, double *a) {
    lapack_int nn = (lapack_int)*n, mm = (lapack_int)*m;
    double *tau = calloc((size_t)mm, sizeof(double));
    (void)LAPACKE_dgeqrf(LAPACK_COL_MAJOR, nn, mm, a, nn, tau);
    (void)LAPACKE_dorgqr(LAPACK_COL_MAJOR, nn, mm, mm, a, nn, tau);
    free(tau);
    return;
}
