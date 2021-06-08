#include "array.h"

/*
  a general matrix of shape (4,3)
*/

int main(void) {
  ARRAY a;
  int length = 12, rank = 2, k = 0;
  int *shape = (int *)calloc((size_t)rank, sizeof(int));
  int *indices = (int *)calloc((size_t)rank, sizeof(int));
  *(shape + 0) = 4;
  *(shape + 1) = 3;
  a.content = (double *)calloc((size_t)length, sizeof(double));
  a.length = &length;
  a.shape = shape;
  a.rank = &rank;
  a.decode = GMDECODE;
  a.encode = GMENCODE;
  printf("========\n");
  printf("GMDECODE\n");
  printf("========\n");
  for (int j = 1; j <= *(a.shape + 1); j++) {
    *(indices + 1) = j;
    for (int i = 1; i <= *(a.shape + 0); i++) {
      *(indices + 0) = i;
      (void)a.decode(indices, &k, a.shape, a.rank);
      printf("indices %4d %4d location %4d\n", i, j, k);
    }
  }
  printf("========\n");
  printf("GMENCODE\n");
  printf("========\n");
  for (int i = 0; i < length; i++) {
    (void)a.encode(&i, indices, a.shape, a.rank);
    printf("location %4d indices %4d %4d\n", i, *(indices + 0), *(indices + 1));
  }
  (void)free(shape);
  (void)free(a.content);
  (void)free(indices);
  return EXIT_SUCCESS;
}
