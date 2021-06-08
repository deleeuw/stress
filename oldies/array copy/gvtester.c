#include "array.h"

/*
  a vector of shape (10)
*/

int main(void) {
  ARRAY a;
  int length = 10, rank = 1, j = 0;
  int *shape = (int *)calloc((size_t)rank, sizeof(int));
  *shape = 10;
  a.content = (double *)calloc((size_t)length, sizeof(double));
  a.length = &length;
  a.shape = shape;
  a.rank = &rank;
  a.decode = GVDECODE;
  a.encode = GVENCODE;
  printf("========\n");
  printf("GVDECODE\n");
  printf("========\n");
  for (int i = 1; i <= length; i++) {
    (void)a.decode(&i, &j, a.shape, a.rank);
    printf("index %4d location %4d\n", i, j);
  }
  printf("========\n");
  printf("GVENCODE\n");
  printf("========\n");
  for (int i = 0; i < length; i++) {
    (void)a.encode(&i, &j, a.shape, a.rank);
    printf("location %4d index %4d\n", i, j);
  }
  (void)free(shape);
  (void)free(a.content);
  return EXIT_SUCCESS;
}
