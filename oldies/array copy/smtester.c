#include "array.h"

/*
  a symmetric matrix of shape (4,4)
*/

int main(void) {
  ARRAY a;
  int length = 10, rank = 2, shape = 4, k = 0;
  int *indices = (int *)calloc((size_t)rank, sizeof(int));
  a.content = (CTYPE *)calloc((size_t)length, sizeof(CTYPE));
  a.length = &length;
  a.shape = &shape;
  a.rank = &rank;
  a.decode = SMDECODE;
  a.encode = SMENCODE;
  printf("========\n");
  printf("SMDECODE\n");
  printf("========\n");
  for (int j = 1; j <= shape; j++) {
    *(indices + 1) = j;
    for (int i = 1; i <= shape; i++) {
      *(indices + 0) = i;
      (void)a.decode(indices, &k, a.shape, a.rank);
      printf("indices %4d %4d location %4d\n", i, j, k);
    }
  }
  printf("========\n");
  printf("SMENCODE\n");
  printf("========\n");
  for (int i = 0; i < length; i++) {
    (void)a.encode(&i, indices, a.shape, a.rank);
    printf("location %4d indices %4d %4d\n", i, *(indices + 0), *(indices + 1));
  }
  (void)free(a.content);
  (void)free(indices);
  return EXIT_SUCCESS;
}
