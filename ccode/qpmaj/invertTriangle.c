#include <stdio.h>
#include <stdlib.h>

void invertLowerTriangle(const int, double *);

static inline int TINDEX(const int, const int, const int);

// lower triangle in compact storage. make sure i >= j.

static inline int TINDEX(const int i, const int j, const int n) {
  return (-(j * (j - (2 * n + 3))) / 2 - n + (i - j));
}

// invert a lower triangular matrix (with nonzero diagonal) in
// compact column-major storage. Or, equivalently, a upper triangular
// matrix in compact row-major st0rage.

void invertLowerTriangle(const int n, double *x) {
  for (int k = 1; k <= n; k++) {
    x[TINDEX(k, k, n)] = 1.0 / x[TINDEX(k, k, n)];
    for (int i = k + 1; i <= n; i++) {
      double s = 0.0;
      for (int j = 1; j <= i - 1; j++) {
        if (j >= k) {
          s += x[TINDEX(i, j, n)] * x[TINDEX(j, k, n)];
        }
      }
      x[TINDEX(i, k, n)] = -s / x[TINDEX(k, k, n)];
    }
  }
}

double x[15] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
int n = 5;

int main() {
  (void)invertLowerTriangle(n, x);
  for (int i = 0; i < 15; i++) {
    printf(" %10.4f ", x[i]);
  }
  printf("\n");
  return (EXIT_SUCCESS);
}
