#include <math.h>

#define MAX(a, b) (((a) > (b)) ? (a) : (b))

void cleanup(double *a, int *n, int *m, int *ind, double *eps) {
    int i, j, l, ii, jj, nn = *n, mm = *m;
    double s;
    for (i = 0; i < (nn - 1); i++) {
        ii = i * mm;
        for (j = (i + 1); j < nn; j++) {
            s = 0.0;
            jj = j * mm;
            if (ind[j] == 0) continue;
            for (l = 0; l < mm; l++) {
                s = MAX(s, fabs(*(a + ii + l) - *(a + jj + l)));
            }
            if (s < *eps) {
              ind[j] = 0;
            }
        }
    }
}
