#include "array.h"

/*
  a general array of shape (2,4,3)
*/

int main(void) {
  ARRAY a;
  int length = 24, rank = 3, l = 0;
  int *shape = (int *)calloc((size_t)rank, sizeof(int));
  int *indices = (int *)calloc((size_t)rank, sizeof(int));
  *(shape + 0) = 2;
  *(shape + 1) = 4;
  *(shape + 2) = 3;
  a.content = (double *)calloc((size_t)length, sizeof(double));
  a.length = &length;
  a.shape = shape;
  a.rank = &rank;
  a.decode = GADECODE;
  a.encode = GAENCODE;
  printf("========\n");
  printf("GADECODE\n");
  printf("========\n");
  for (int k = 1; k <= *(a.shape + 2); k++) {
    *(indices + 2) = k;
    for (int j = 1; j <= *(a.shape + 1); j++) {
      *(indices + 1) = j;
      for (int i = 1; i <= *(a.shape + 0); i++) {
        *(indices + 0) = i;
        (void)a.decode(indices, &l, a.shape, a.rank);
        printf("indices %4d %4d %4d location %4d\n", i, j, k, l);
      }
    }
  }
  printf("========\n");
  printf("GAENCODE\n");
  printf("========\n");
  for (int i = 0; i < length; i++) {
    (void)a.encode(&i, indices, a.shape, a.rank);
    printf("location %4d indices %4d %4d %4d\n", i, *(indices + 0),
           *(indices + 1), *(indices + 2));
  }
  (void)free(shape);
  (void)free(a.content);
  (void)free(indices);
}
