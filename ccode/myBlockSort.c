
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int myComp(const void *, const void *);
void myBlockSort(const double *, const int, couple *, block *);

typedef struct couple {
    double value;
    int index;
} couple;

typedef struct block {
    double value;
    int size;
    int rank;
} block;

int myComp(const void *px, const void *py) {
    double x = ((struct couple *)px)->value;
    double y = ((struct couple *)py)->value;
    return (int)copysign(1.0, x - y);
}

void myBlockSort(const double *x, const int n, couple *xi, block *yi) {
    for (int i = 0; i < n; i++) {
        xi[i].value = x[i];
        xi[i].index = i;
    }
    (void)qsort(xi, (size_t)n, (size_t)sizeof(couple), myComp);
    int k = 0, r = 0, *m = 0;
    for (int i = 0; i < n - 1; i++) {
        block[r].value = x[i].value;
        block[r].rank = r;
        if (xi[i].value == xi[i+1].value) {
            block[r].size += 1;            
        } else {
            r += 1;
        }
    }
    for (int i = 0; i < n; i++) {
        printf("value %4.2f size %2d, rank %2d\n", 
            block[i].value, block[i].size, block[i].rank);
    }
    return;
}

double x[10] = {3.0, 1.0, 1.0, 5.0. 1.0. 5.0, 1.0, 2.0, 5.0, 2.0};
int n = 10;

int main () {
        couple *xi =
        (struct couple *)calloc((size_t)n, (size_t)sizeof(couple));
        block *yi =
        (struct couple *)calloc((size_t)n, (size_t)sizeof(block));

    (void) myBlockSort (x, n, xi, yi);
    free(y1)
}
