#ifndef SMACOF_H
#define SMACOF_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

static inline int VINDEX(const int);
static inline int MINDEX(const int, const int, const int);
static inline int SINDEX(const int, const int, const int);

static inline double SQUARE(const double);

static inline double MAX(const double, const double);
static inline double MIN(const double, const double);
static inline int IMIN(const int, const int);
static inline int IMAX(const int, const int);

typedef enum transform {
    ratio = 1,
    interval = 2,
    polynomial = 3,
    splinical = 4
} transform;

typedef enum ties { primary = 1, secondary = 2, tertiary = 3 } ties;

typedef struct couple {
    double value;
    int index;
} couple;

typedef struct block {
    double value;
    int size;
    int rank;
} block;

void primat(const int, const int, const int, const int, const double *);
void pritrl(const int, const int, const int, const double *);
void smacofDist(const double *, const int, const int, double *);
void smacofDoubleCenter(double *, const int);
void smacofDCMultX(const double *, double *, const int, const int);
int myComp(const void *, const void *);
void myBlockSort(const double *, const int, int *, couple *, block *);

double smacofLossW(const double *, const double *, const double *, const int);
double smacofLossU(const double *, const double *, const int);

void smacofNormDeltaW(const double *, double *, const int);
void smacofNormDeltaU(double *, const int);

void smacofScaleXW(double *, double *, const double *, const double *,
                   const int, const int);
void smacofScaleXU(double *, double *, const double *, const int, const int);

void smacofInvertVW(const double *, double *, const int, const int);

void smacofGuttmanW(const double *, const double *, const double *,
                    const double *, const int, const int, const bool, double *);
void smacofGuttmanU(const double *, const double *, const int, const int,
                    const bool, double *);
void smacofGradientW(double *dist, double *w, double *delta, const int n,
                     const int p, double *x);
void smacofGradientU(double *dist, double *delta, const int n, const int p,
                     double *x);
void smacofIntervalU(const double *, double *, int);

void smacofIntervalW(const double *, const double *, double *, int);

void smacofMissingU(const double *, double *, int);

void smacofSSWR(double *, const double *, double *, const int *, const int *,
                const bool *, const int *, const double *, const bool *);
void smacofSSUR(double *, double *, const int *, const int *, const bool *,
                const int *, const double *, const bool *);
void smacofSSUI(double *, double *x, const int *, const int *, const bool *,
                const int *, const double *, const bool *);
void smacofSSWI(double *, const double *, double *, const int *, const int *,
                const bool *, const int *, const double *, const bool *);

// VINDEX takes 1,...,n to 0,...,n-1

static inline int VINDEX(const int i) { return i - 1; }

// MINDEX retrieves element (i,j) from an n x m matrix in
// column-major-order strorage

static inline int MINDEX(const int i, const int j, const int n) {
    return (i - 1) + (j - 1) * n;
}

//
// SINDEX retrieves element (i, j) from a strict lower triangular matrix
// of order n. Thus always i > j. Loop
//  for (int j = 1; j < (n - 1); j++) {
//    for (int i = (j + 1); i<= n; i++) {
//    }
//..}
//

static inline int SINDEX(const int i, const int j, const int n) {
    return ((j - 1) * n) - (j * (j - 1) / 2) + (i - j) - 1;
}

static inline double SQUARE(const double x) { return x * x; }

static inline double MAX(const double x, const double y) {
    return (x > y) ? x : y;
}

static inline double MIN(const double x, const double y) {
    return (x < y) ? x : y;
}

static inline int IMAX(const int x, const int y) { return (x > y) ? x : y; }

static inline int IMIN(const int x, const int y) { return (x < y) ? x : y; }

#endif /* SMACOF_H */
