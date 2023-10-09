// This implements the formulas in
//
// D.J. Jeffrey
// Formulae, Algorithms, and Quartic Extrema
// Mathematics Magazine, 1997, 70(5), 341-348
//
// to find the minimum of a quartic polynomial.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SQ(x) ((x) * (x))
#define CU(x) ((x) * (x) * (x))
#define QU(x) ((x) * (x) * (x) * (x))

void smacofJeffrey(double *a, double *minwhere, double *minvalue) {
    if (a[4] <= 0.0) {
        printf("Quartic term must be positive. Exiting...\n");
        exit (EXIT_FAILURE);
    }
    double b0 = 0.0, b1 = 0.0, b2 = 0.0;
    b2 += a[2] / (3.0 * a[4]);
    b2 -= (SQ(a[3])) / (8.0 * SQ(a[4]));
    b1 += a[1] / (2.0 * a[4]);
    b1 += CU(a[3]) / (16.0 * CU(a[4]));
    b1 -= (a[2] * a[3]) / (4.0 * SQ(a[4]));
    b0 += (a[2] * SQ(a[3])) / (16.0 * SQ(a[4]));
    b0 -= (a[1] * a[3]) / (4.0 * a[4]);
    b0 -= (3.0 * QU(a[3])) / (256.0 * CU(a[4]));
    if (fabs(b1) > 1e-15) {
        double s = SQ(b1) + CU(b2) + sqrt(QU(b1) + 2 * SQ(b1) * CU(b2));
        double t = pow(s, (double) 1 / 3);
        double k = t + (SQ(b2) / t) + b2;
        double infp = -0.75 * (k - b2) * (k - 3.0 * b2);
        *minwhere = -b1 / k - 0.25 * a[3] / a[4];
        *minvalue = infp * a[4] + a[0] + b0;
    } else {
        double infp = -2.25 * SQ(fmin(0.0, b2));;
        *minwhere = sqrt(-fmin(0.0, 1.5 * b2));
        *minvalue = infp * a[4] + a[0] + b0;
    }
    return;
}
