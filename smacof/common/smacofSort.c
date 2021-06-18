
#include "../include/smacof.h"

int myComp(const void *px, const void *py) {
  double x = ((struct couple *)px)->value;
  double y = ((struct couple *)py)->value;
  return (int)copysign(1.0, x - y);
}

void myBlockSort(const double *x, const int n, int *nblock, couple *xi, block *yi) {
  for (int i = 0; i < n; i++) {
    xi[i].value = x[i];
    xi[i].index = i;
  }
  (void)qsort(xi, (size_t)n, (size_t)sizeof(couple), myComp);
  int k = 0, r = 1;
  for (int i = 0; i < n - 1; i++) {
    yi[i].size = 1;
  }
  for (int i = 0; i < n - 1; i++) {
    yi[VINDEX(r)].value = xi[i].value;
    yi[VINDEX(r)].rank = r;
    if (xi[i].value == xi[i + 1].value) {
      yi[VINDEX(r)].size += 1;
    } else {
      r += 1;
    }
  }
  *nblock = r;
  return;
}

/*
double x[10] = {3.0, 1.0, 1.0, 5.0, 1.0, 5.0, 1.0, 2.0, 5.0, 2.0};
int n = 10;

int main() {
  couple *xi = (couple *)calloc((size_t)n, (size_t)sizeof(couple));
  block *yi = (block *)calloc((size_t)n, (size_t)sizeof(block));
  int m = 0;
  (void)myBlockSort(x, n, &m, xi, yi);
  yi = (block *)realloc((block *)yi, (size_t) m * sizeof(block));
  for (int i = 0; i < m; i++) {
    printf("value %4.2f size %2d, rank %2d\n", yi[i].value, yi[i].size,
           yi[i].rank);
  }
  free(yi);
  free(xi);
  return (EXIT_SUCCESS);
}
*/
