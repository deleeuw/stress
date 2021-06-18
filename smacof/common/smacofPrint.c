#include "../include/smacof.h"

// print a general matrix

void smacofPrintMatrix(const int n, const int m, const int w, const int p,
                       const double *x) {
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= m; j++) {
      printf(" %+*.*f ", w, p, x[MINDEX(i, j, n)]);
    }
    printf("\n");
  }
  printf("\n\n");
  return;
}

// print strict lower triangle

void smacofPrintTriangle(const int n, const int w, const int p,
                         const double *x) {
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= i; j++) {
      if (i > j) {
        printf(" %+*.*f ", w, p, x[SINDEX(i, j, n)]);
      }
    }
    printf("\n");
  }
  printf("\n\n");
  return;
}
