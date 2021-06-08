#include <stdlib.h>
#include <stdio.h>


extern void mySortDouble (double *, double *, double *, int *, const int *);
extern void mySortInteger (int *, int *, const int *);
extern void tieBlock (double *, int *, const int *, int *);
extern void makeBlocks (double *, double *, double *, double *, int *, const int *, const int *);
extern void sortBlocks (double *, double *, int *, const int *, const int *, const int *);
extern void jbkPava (double *, double *, const int *);
void primary (double *, double *, int *, int *, const int *, const int *);
void secondary (double *, double *, int *, const int *, const int *);
void tertiary (double *, double *, int *, const int *, const int *);
void monreg (double *, double *, double *, const int *, const int *);

void monreg (double *x, double *y, double *w, const int *n, const int *ties) {
    int *iblks = (int *) calloc((size_t) * n, sizeof(int));
    int *xind = (int *) calloc((size_t) * n, sizeof(int));
    int *iind = (int *) calloc((size_t) * n, sizeof(int));
    int nblk = 0;
    /*
      * use predictor to sort predictor, weights, and outcome
      */
    (void) mySortDouble (x, y, w, xind, n);
    /*
     * compute predictor tieblocks
     */
    (void) tieBlock (x, iblks, n, &nblk);
    /*
     * do the monotone regression
     */
     for (int i = 0; i <*n; i++)
        printf (" %10.4f\n", y[i]);
    printf("++++++++\n");
    if (*ties == 1)
        (void) primary (y, w, xind, iblks, &nblk, n);
    if (*ties == 2)
        (void) secondary (y, w, iblks, &nblk, n);
    if (*ties == 3)
        (void) tertiary (y, w, iblks, &nblk, n);
    /*
     * use the inverse permutation to put dependent variable back in order
     */
    for (int i = 0; i <*n; i++)
        printf (" %10.4f\n", y[i]);
    printf("++++++++\n");
    (void) mySortInteger (xind, iind, n);
    for (int i = 0; i < *n; i++) {
        w[i] = y[iind[i]];
    }
    free (iblks);
    free (xind);
    free (iind);
    return;
}

void primary (double *y, double *w, int *xind, int *iblks, const int *n, const int *nblk) {
    (void) sortBlocks (y, w, xind, iblks, n, nblk);
    (void) jbkPava (y, w, n);
    return;
}

void secondary (double *y, double *w, int *iblks, const int *nblk, const int *n) {
    double *yblocks = (double *) calloc((size_t) * nblk, (size_t) sizeof(double));
    double *wblocks = (double *) calloc((size_t) * nblk, (size_t) sizeof(double));
    /*
     * compute block weights and averages
     */
    (void) makeBlocks (y, w, yblocks, wblocks, iblks, n, nblk);
    /*
     * do the monotone regression on the blocks
     */
    (void) jbkPava (yblocks, wblocks, nblk);
    /*
     * put the block values make in the result vector
     */
    for (int i = 0; i < *n; i++)
        y[i] = yblocks[iblks[i] - 1];
    free (yblocks);
    free (wblocks);
    return;
}

void tertiary (double *y, double *w, int *iblks, const int *nblk, const int *n) {
    double *yblocks = (double *) calloc((size_t) * nblk, (size_t) sizeof(double));
    double *wblocks = (double *) calloc((size_t) * nblk, (size_t) sizeof(double));
    double *zblocks = (double *) calloc((size_t) * nblk, (size_t) sizeof(double));
    /*
     * compute block weights and averages
     */
    (void) makeBlocks (y, w, yblocks, wblocks, iblks, n, nblk);
    for (int i = 0; i < *nblk; i++) {
        zblocks[i] = yblocks[i];
    }
    /*
     * do the monotone regression on the blocks
     */
    (void) jbkPava (yblocks, wblocks, nblk);
    /*
     * put the block values make in the result vector
     */
    for (int i = 0; i < *n; i++)
        y[i] = y[i] + zblocks[iblks[i] - 1] - yblocks[iblks[i] - 1];
    free (yblocks);
    free (wblocks);
    free (zblocks);
    return;
}

