
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define DEBUG false

struct block {
    double value;
    double weight;
    int size;
    int previous;
    int next;
};

struct quadruple {
    double value;
    double result;
    double weight;
    int index;
};

struct triple {
    double value;
    double weight;
    int index;
};

struct pair {
    int value;
    int index;
};

int myCompDouble(const void *, const void *);
int myCompInteger(const void *, const void *);
void mySortDouble(double *, double *, double *, int *, const int *);
void mySortInteger(int *, int *, const int *);
void mySortInBlock(double *, double *, int *, int *);
void tieBlock(double *, int *, const int *, int *);
void makeBlocks(double *, double *, double *, double *, int *, const int *,
                const int *);
void sortBlocks(double *, double *, int *, const int *, const int *,
                const int *);
void jbkPava(double *, double *, const int *);
void primary(double *, double *, int *, int *, const int *, const int *);
void secondary(double *, double *, int *, const int *, const int *);
void tertiary(double *, double *, int *, const int *, const int *);
void monreg(double *, double *, double *, const int *, const int *);

int myCompDouble(const void *px, const void *py) {
    double x = ((struct quadruple *)px)->value;
    double y = ((struct quadruple *)py)->value;
    return (int)copysign(1.0, x - y);
}

int myCompInteger(const void *px, const void *py) {
    int x = ((struct pair *)px)->value;
    int y = ((struct pair *)py)->value;
    return (int)copysign(1.0, x - y);
}

void mySortInBlock(double *x, double *w, int *xind, int *n) {
    int nn = *n;
    struct triple *xi =
        (struct triple *)calloc((size_t)nn, (size_t)sizeof(struct triple));
    for (int i = 0; i < nn; i++) {
        xi[i].value = x[i];
        xi[i].weight = w[i];
        xi[i].index = xind[i];
    }
    (void)qsort(xi, (size_t)nn, (size_t)sizeof(struct triple), myCompDouble);
    for (int i = 0; i < nn; i++) {
        x[i] = xi[i].value;
        w[i] = xi[i].weight;
        xind[i] = xi[i].index;
    }
    free(xi);
    return;
}

void mySortDouble(double *x, double *y, double *w, int *xind, const int *n) {
    int nn = *n;
    struct quadruple *xi = (struct quadruple *)calloc(
        (size_t)nn, (size_t)sizeof(struct quadruple));
    for (int i = 0; i < nn; i++) {
        xi[i].value = x[i];
        xi[i].result = y[i];
        xi[i].weight = w[i];
        xi[i].index = i + 1;
    }
    (void)qsort(xi, (size_t)nn, (size_t)sizeof(struct quadruple), myCompDouble);
    for (int i = 0; i < nn; i++) {
        x[i] = xi[i].value;
        y[i] = xi[i].result;
        w[i] = xi[i].weight;
        xind[i] = xi[i].index;
    }
    free(xi);
    return;
}

void mySortInteger(int *x, int *k, const int *n) {
    int nn = *n;
    struct pair *xi =
        (struct pair *)calloc((size_t)nn, (size_t)sizeof(struct pair));
    for (int i = 0; i < nn; i++) {
        xi[i].value = x[i];
        xi[i].index = i + 1;
    }
    (void)qsort(xi, (size_t)nn, (size_t)sizeof(struct pair), myCompInteger);
    for (int i = 0; i < nn; i++) {
        x[i] = xi[i].value;
        k[i] = xi[i].index;
    }
    free(xi);
    return;
}

void tieBlock(double *x, int *iblks, const int *n, int *nblk) {
    iblks[0] = 1;
    for (int i = 1; i < *n; i++) {
        if (x[i - 1] == x[i]) {
            iblks[i] = iblks[i - 1];
        } else {
            iblks[i] = iblks[i - 1] + 1;
        }
    }
    *nblk = iblks[*n - 1];
    return;
}

void makeBlocks(double *x, double *w, double *xblks, double *wblks, int *iblks,
                const int *n, const int *nblk) {
    for (int i = 0; i < *nblk; i++) {
        xblks[i] = 0.0;
        wblks[i] = 0.0;
    }
    for (int i = 0; i < *n; i++) {
        xblks[iblks[i] - 1] += w[i] * x[i];
        wblks[iblks[i] - 1] += w[i];
    }
    for (int i = 0; i < *nblk; i++) {
        xblks[i] = xblks[i] / wblks[i];
    }
    return;
}

void sortBlocks(double *y, double *w, int *xind, const int *iblks, const int *n,
                const int *nblk) {
    int *nblks = (int *)calloc((size_t)*nblk, sizeof(int));
    for (int i = 0; i < *n; i++) {
        nblks[iblks[i] - 1]++;
    }
    int k = 0;
    for (int i = 0; i < *nblk; i++) {
        int nn = nblks[i];
        (void)mySortInBlock(y + k, w + k, xind + k, &nn);
        k += nn;
    }
    free(nblks);
    return;
}

void jbkPava(double *x, double *w, const int *n) {
    struct block *blocks = calloc((size_t)*n, sizeof(struct block));
    for (int i = 0; i < *n; i++) {
        blocks[i].value = x[i];
        blocks[i].weight = w[i];
        blocks[i].size = 1;
        blocks[i].previous = i - 1;  // index first element previous block
        blocks[i].next = i + 1;      // index first element next block
    }
    int active = 0;
    do {
        bool upsatisfied = false;
        int next = blocks[active].next;
        if (next == *n)
            upsatisfied = true;
        else if (blocks[next].value > blocks[active].value)
            upsatisfied = true;
        if (!upsatisfied) {
            double ww = blocks[active].weight + blocks[next].weight;
            int nextnext = blocks[next].next;
            blocks[active].value =
                (blocks[active].weight * blocks[active].value +
                 blocks[next].weight * blocks[next].value) /
                ww;
            blocks[active].weight = ww;
            blocks[active].size += blocks[next].size;
            blocks[active].next = nextnext;
            if (nextnext < *n) blocks[nextnext].previous = active;
            blocks[next].size = 0;
        }
        bool downsatisfied = false;
        int previous = blocks[active].previous;
        if (previous == -1)
            downsatisfied = true;
        else if (blocks[previous].value < blocks[active].value)
            downsatisfied = true;
        if (!downsatisfied) {
            double ww = blocks[active].weight + blocks[previous].weight;
            int previousprevious = blocks[previous].previous;
            blocks[active].value =
                (blocks[active].weight * blocks[active].value +
                 blocks[previous].weight * blocks[previous].value) /
                ww;
            blocks[active].weight = ww;
            blocks[active].size += blocks[previous].size;
            blocks[active].previous = previousprevious;
            if (previousprevious > -1) blocks[previousprevious].next = active;
            blocks[previous].size = 0;
        }
        if ((blocks[active].next == *n) && downsatisfied) break;
        if (upsatisfied && downsatisfied) active = next;
    } while (true);
    int k = 0;
    for (int i = 0; i < *n; i++) {
        int blksize = blocks[i].size;
        if (blksize > 0.0) {
            for (int j = 0; j < blksize; j++) {
                x[k] = blocks[i].value;
                k++;
            }
        }
    }
    free(blocks);
}
