#include "array.h"

int binCoef(const int n, const int m) {
  int result;
  int *work = (int *)calloc((size_t)(m + 1), sizeof(int));
  work[0] = 1;
  for (int i = 1; i <= n; i++) {
    for (int j = (i > m) ? m : i; j > 0; j--) {
      work[j] += work[j - 1];
    }
  }
  result = work[m];
  free(work);
  return (result);
}

void isort(int *cell, const int *n) {
  (void)qsort(cell, (size_t)*n, sizeof(int), icmp);
  return;
}

int icmp(const void *x, const void *y) {
  const int *ix = (const int *)x;
  const int *iy = (const int *)y;
  return (*ix - *iy);
}

void GVENCODE(const int *location, int *indices, const int *shape,
              const int *rank) {
  *indices = *location + 1;
  return;
}

void GVDECODE(const int *indices, int *location, const int *shape,
              const int *rank) {
  *location = *indices - 1;
  return;
}

void GMENCODE(const int *location, int *indices, const int *shape,
              const int *rank) {
  int nrow = *(shape + 0), icol, irow, iloc = *location;
  irow = (iloc % nrow) + 1;
  icol = (iloc / nrow) + 1;
  *(indices + 0) = irow;
  *(indices + 1) = icol;
  return;
}

void GMDECODE(const int *indices, int *location, const int *shape,
              const int *rank) {
  int nrow = *(shape + 0), irow = *(indices + 0), icol = *(indices + 1),
      iloc = (irow - 1) + (icol - 1) * nrow;
  *location = iloc;
  return;
}

void SMDECODE(const int *indices, int *location, const int *shape,
              const int *rank) {
  int nrow = *(shape + 0), irow = *(indices + 0), icol = *(indices + 1),
      iloc = 0;
  if (irow < icol) {
    iloc = irow;
    irow = icol;
    icol = iloc;
  }
  int *maxcol = (int *)calloc((size_t)nrow + 1, sizeof(int));
  for (int i = 0; i <= nrow; i++) {
    *(maxcol + i) = i * nrow - i * (i - 1) / 2;
  }
  iloc = *(maxcol + (icol - 1)) + (irow - icol + 1);
  *location = iloc - 1;
  free(maxcol);
  return;
}

void SMENCODE(const int *location, int *indices, const int *shape,
              const int *rank) {
  int nrow = *(shape + 0), irow, icol, iloc = *location + 1;
  int *maxcol = (int *)calloc((size_t)nrow + 1, sizeof(int));
  for (int i = 0; i <= nrow; i++) {
    *(maxcol + i) = i * nrow - i * (i - 1) / 2;
  }
  for (int i = 1; i <= nrow; i++) {
    if (iloc <= *(maxcol + i)) {
      icol = i;
      irow = iloc - *(maxcol + (i - 1)) + (i - 1);
      break;
    }
  }
  *(indices + 0) = irow;
  *(indices + 1) = icol;
  free(maxcol);
  return;
}

void GADECODE(const int *indices, int *location, const int *shape,
              const int *rank) {
  int aux = 1, res = 1, r = *rank;
  for (int i = 0; i < r; i++) {
    res += (*(indices + i) - 1) * aux;
    aux *= *(shape + i);
  }
  *location = res - 1;
  return;
}

void GAENCODE(const int *location, int *indices, const int *shape,
              const int *rank) {
  int r = *rank, iloc = *location + 1, aux = 1, k = 1;
  for (int j = 1; j < r; j++) {
    aux *= *(shape + j - 1);
  }
  for (int j = r - 1; j > 0; j--) {
    k = (iloc - 1) / aux;
    *(indices + j) = k + 1;
    iloc -= k * aux;
    aux /= *(shape + j - 1);
  }
  *(indices + 0) = iloc;
  return;
}

void SADECODE(const int *indices, int *location, const int *shape,
              const int *rank) {
  int iloc = 1, k = 1, r = *rank;
  int *ind = (int *)calloc((size_t)r, sizeof(int));
  for (int i = 0; i < r; i++) {
    *(ind + i) = *(indices + i);
  }
  (void)isort(ind, rank);
  for (int i = 1; i <= r; i++) {
    k = i + (*(ind + i - 1) - 1) - 1;
    iloc += binCoef(k, i);
  }
  *location = iloc - 1;
  (void)free(ind);
  return;
}

void SAENCODE(const int *location, int *indices, const int *shape,
              const int *rank) {
  int n = *shape, r = *rank, iloc = *location + 1, v = iloc - 1;
  for (int k = r; k >= 1; k--) {
    for (int j = 0; j < n; j++) {
      int sj = binCoef(k + j - 1, k), sk = binCoef(k + j, k);
      if (v < sk) {
        indices[k - 1] = j + 1;
        v -= sj;
        break;
      }
    }
  }
  return;
}
