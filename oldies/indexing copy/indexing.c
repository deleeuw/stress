
#include "indexing.h"

int binCoef(const int n, const int m) {
    int *work = (int *)calloc((size_t)m + 1, sizeof(int));
    work[0] = 1;
    for (int i = 1; i <= n; i++) {
        for (int j = IMIN(i, m); j > 0; j--) {
            work[j] = work[j] + work[j - 1];
        }
    }
    int choose = work[m];
    free(work);
    return (choose);
}

int int_cmp(const void *x, const void *y) {
    const int *ix = (const int *)x;
    const int *iy = (const int *)y;
    return (*ix - *iy);
}

int fArrayFirst(const int *cell, const int *dimension, const int rank) {
    int aux = 1, result = 1;
    for (int i = 1; i <= rank; i++) {
        result += (cell[VINDEX(i)] - 1) * aux;
        aux *= dimension[VINDEX(i)];
    }
    return (result);
}

int fSupSymIncreasing(int *cell, const int rank) {
    int result = 1;
    (void)qsort(cell, (size_t)rank, sizeof(int), int_cmp);
    for (int i = 1; i <= rank; i++) {
        int k = i + (cell[VINDEX(i)] - 1) - 1;
        result += binCoef(k, i);
    }
    return (result);
}

void fArrayFirstInverse(const int *rank, const int *dimension, const int *index,
                        int *cell) {
    int m = *rank, l = *index, b = 1;
    for (int j = 1; j < m; j++) {
        b *= dimension[VINDEX(j)];
    }
    for (int j = m - 1; j > 0; j--) {
        int k = (l - 1) / b;
        cell[j] = k + 1;
        l -= k * b;
        b /= dimension[j - 1];
    }
    cell[VINDEX(1)] = l;
    return;
}

void fSupSymIncreasingFirstInverse(const int *dimension, const int *rank,
                                   const int *index, int *cell) {
    int n = *dimension, m = *rank, l = *index, v = l - 1;
    for (int k = m; k >= 1; k--) {
        for (int j = 0; j < n; j++) {
            int sj = binCoef(k + j - 1, k), sk = binCoef(k + j, k);
            if (v < sk) {
                cell[VINDEX(k)] = j + 1;
                v -= sj;
                break;
            }
        }
    }
    return;
}

void binCoefGlue(const int *n, const int *m, int *choose) {
    *choose = binCoef(*n, *m);
}

void fArrayFirstGlue(const int *cell, const int *dimension, const int *rank,
                     int *result) {
    *result = fArrayFirst(cell, dimension, *rank);
}

void fSupSymIncreasingGlue(int *cell, const int *rank, int *result) {
    *result = fSupSymIncreasing(cell, *rank);
}