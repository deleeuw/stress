#include "smatrix.h"
#include "sroutines.h"

/*
 * The routine scholesky() replaces a positive definite matrix A of order n
 * in lower-triangular column-major storage by its Cholesky factor L. Thus
 * A = LL'. The determinent of A is returned as well.
 *
 * If a pivot is zero the routine jumps over the pivot and continues,
 * and sets the singularity warning. Zero of a pivot means less than a small
 * positive eps.
 *
 * If a pivot is negative the routine stops, and sets the indefinite warning.
 * Negativity of a pivot means less than -eps.
 */

void scholesky(const int *nn, double *a, bool *singular, bool *indefinite,
               int *npiv, double *ppiv, const double *eps) {
  int n = *nn;
  double t = 0.0, s = 0.0;
  *npiv = 0;
  *ppiv = 1.0;
  for (int k = 1; k <= n; k++) {
    t = a[TINDEX(k, k, n)];
    if (t < -*eps) {
      *indefinite = true;
      return;
    }
    if (t < *eps) {
      *singular = true;
      continue;
    }
    (*ppiv) *= t;
    (*npiv)++;
    s = sqrt(t);
    for (int i = k; i <= n; i++) {
      a[TINDEX(i, k, n)] = a[TINDEX(i, k, n)] / s;
    }
    for (int i = k + 1; i <= n; i++) {
      for (int j = k + 1; j <= i; j++) {
        a[TINDEX(i, j, n)] -= a[TINDEX(i, k, n)] * a[TINDEX(j, k, n)];
      }
    }
  }
  return;
}

/*
 * The routine pcholesky() replaces a pivoted version of a positive
 * semi-definite matrix A of order n in lower-triangular column-major storage by
 * its Cholesky factor L. Thus we find a permutation P of rows and columns of A
 * such that P'AP=LL'.
 *
 * Pivoting is on the maximum diagonal element. If this is negative the matrix
 * was not positive semidefinite and we return the indefinite flag.  Negatie for
 * a pivot means less than minus a small positive eps. The approximate rank is
 * returned as the number of nonzero pivots, and the dterminant is returned as
 * the product of pivots. On exit the vector order has the pivot permutation.
 */

void pcholesky(const int *nn, double *a, double *ppiv, int *npiv, int *order,
               bool *indefinite, const double *eps) {
  int n = *nn, q = 0, itp = 0;
  double t = 0.0, at = 0.0, tmp = 0.0, tp1 = 0.0, tp2 = 0.0, tp3 = 0.0, s = 0.0;
  *ppiv = 1.0;
  *npiv = 0;
  for (int i = 1; i <= n; i++) {
    order[VINDEX(i)] = i;
  }
  for (int k = 1; k <= n; k++) {
    t = 0.0;
    q = 0;
    for (int i = k; i <= n; i++) {
      at = a[TINDEX(i, i, n)];
      if (at >= t) {
        t = at;
        q = i;
      }
    }
    if (t < -*eps) {
      *indefinite = true;
      return;
    }
    (*npiv)++;
    (*ppiv) *= t;
    itp = order[VINDEX(k)];
    order[VINDEX(k)] = order[VINDEX(q)];
    order[VINDEX(q)] = itp;
    tp1 = a[TINDEX(k, k, n)];
    tp2 = a[TINDEX(q, q, n)];
    tp3 = a[TINDEX(q, k, n)];
    for (int i = 1; i <= k - 1; i++) {
      tmp = a[TINDEX(k, i, n)];
      a[TINDEX(k, i, n)] = a[TINDEX(q, i, n)];
      a[TINDEX(q, i, n)] = tmp;
    }
    for (int i = k + 1; i <= q - 1; i++) {
      tmp = a[TINDEX(i, k, n)];
      a[TINDEX(i, k, n)] = a[TINDEX(q, i, n)];
      a[TINDEX(q, i, n)] = tmp;
    }
    for (int i = q + 1; i <= n; i++) {
      tmp = a[TINDEX(i, q, n)];
      a[TINDEX(i, q, n)] = a[TINDEX(i, k, n)];
      a[TINDEX(i, k, n)] = tmp;
    }
    a[TINDEX(q, q, n)] = tp1;
    a[TINDEX(k, k, n)] = tp2;
    a[TINDEX(q, k, n)] = tp3;
    s = sqrt(t);
    for (int i = k; i <= n; i++) {
      a[TINDEX(i, k, n)] /= s;
    }
    for (int i = k + 1; i <= n; i++) {
      for (int j = k + 1; j <= i; j++) {
        a[TINDEX(i, j, n)] -= a[TINDEX(i, k, n)] * a[TINDEX(j, k, n)];
      }
    }
  }
  return;
}

void bcholesky(const int *nn, double *a, double *d, bool *singular,
               double *ppiv, int *npiv, const double *eps) {
  int n = *nn, wp = 8, pp = 4, r = 0, p = 0, s = 0, k = 1;
  double t = 0.0, mu = (1.0 + sqrt(17.0)) / 8.0, lbd = 0.0, lba = 0.0,
         sig = 0.0, sia = 0.0;
  *npiv = 0;
  *ppiv = 1.0;
  *singular = false;
  while (k <= n) {
    r = 0;
    lbd = 0.0;
    for (int i = k + 1; i <= n; i++) {
      lba = fabs(a[TINDEX(i, k, n)]);
      if (lba > lbd) {
        lbd = lba;
        r = i;
      }
      if (a[TINDEX(k, k, n)] > (mu * lbd)) {
        d[VINDEX(k)] = sqrt(a[TINDEX(k, k, n)]);
        s = 1;
      } else {
        p = 0;
        sig = 0.0;
        for (int i = 1; i < r; i++) {
          sia = fabs(a[TINDEX(r, i, n)]);
          if (sia > sig) {
            sig = sia;
            p = i;
          }
        }
        for (int i = r + 1; i <= n; i++) {
          sia = fabs(a[TINDEX(i, r, n)]);
          if (sia > sig) {
            sig = sia;
            p = i;
          }
        }
        if ((sig * fabs(a[TINDEX(k, k, n)])) >= (mu * SQUARE(lbd))) {
          d[VINDEX(k)] = sqrt(a[TINDEX(k, k, n)]);
          s = 1;
        } else if (fabs(a[TINDEX(r, r, n)]) >= (mu * sig)) {
          /* permute r and k */
          d[VINDEX(k)] = sqrt(a[TINDEX(k, k, n)]);
          s = 1;
        } else {
          /* permute (r,p) to (k,k+1) */
          /* choose 2 x 2 pivot */
          s = 2;
        }
      }
      for (int i = k; i <= n; i++) {
        a[TINDEX(i, k, n)] = a[TINDEX(i, k, n)] / t;
      }
      for (int i = k + 1; i <= n; i++) {
        for (int j = k + 1; j <= i; j++) {
          a[TINDEX(i, j, n)] -= t * a[TINDEX(i, k, n)] * a[TINDEX(j, k, n)];
        }
      }
      (void)pritru(nn, &wp, &pp, a);
    }
  }
  return;
}
