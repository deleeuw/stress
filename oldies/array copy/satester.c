#include "array.h"

/*
  a symmetric array of shape (3,3,3)
*/

int main(void) {
  ARRAY a;
  int length = 10, rank = 3, shape = 3, l = 0;
  int *indices = (int *)calloc((size_t)rank, sizeof(int));
  a.content = (double *)calloc((size_t)length, sizeof(double));
  a.length = &length;
  a.shape = &shape;
  a.rank = &rank;
  a.decode = SADECODE;
  a.encode = SAENCODE;
  printf("========\n");
  printf("SADECODE\n");
  printf("========\n");
  for (int k = 1; k <= shape; k++) {
    *(indices + 2) = k;
    for (int j = 1; j <= shape; j++) {
      *(indices + 1) = j;
      for (int i = 1; i <= shape; i++) {
        *(indices + 0) = i;
        if ((i <= j) && (j <= k)) {
          (void)a.decode(indices, &l, a.shape, a.rank);
          printf("indices %4d %4d %4d location %4d\n", i, j, k, l);
        }
      }
    }
  }
  printf("========\n");
  printf("SAENCODE\n");
  printf("========\n");
  for (int i = 0; i < length; i++) {
    (void)a.encode(&i, indices, a.shape, a.rank);
    printf("location %4d indices %4d %4d %4d\n", i, *(indices + 0),
           *(indices + 1), *(indices + 2));
  }
  (void)free(a.content);
  (void)free(indices);
  return EXIT_SUCCESS;
}
