/*

smacofBlockSort() is a routine to transform a matrix of dissimilarities (in
lower-triangular column-major storage) into a ordinal multidimensional scaling
structure (OMDS structure). An OMDS structure is an array of tie-blocks, where
each tie-block corresponds with a unique dissimilarity value. Tie-blocks have a
value, a size, a vector of weights, and a vector of indices. They are strictly 
ordered by increasing value. The routine is written for the smacof project, 
but it can be used as a preprocessor for any monotone regression problem.

*/

#include "../include/smacof.h"


int sortComp(const void *px, const void *py) {
    double x = ((struct triple *)px)->value;
    double y = ((struct triple *)py)->value;
    return (int)copysign(1.0, x - y);
}

void smacofBlockSort(const double *x, const double *w, const int n, int nblock,
                 block *yi) {
    triple *xi = (triple *)calloc((size_t)n, sizeof(triple));
    yi = (struct block *)calloc((size_t)n, (size_t)sizeof(block));
    for (int i = 0; i < n; i++) {
        xi[i].value = x[i];
        xi[i].weight = w[i];
        xi[i].index = i;
    }
    (void)qsort(xi, (size_t)n, (size_t)sizeof(triple), sortComp);
    int counter = 0;
    while (counter < n) {
        double value = xi[counter].value;
        int size = 0;
        for (int j = counter; j < n; j++) {
            if (xi[j].value == value) {
                size += 1;
            } else {
                break;
            }
        }
        yi[nblock].size = size;
        yi[nblock].value = value;
        yi[nblock].indices = (int *)calloc((size_t)size, sizeof(int));
        yi[nblock].weights = (double *)calloc((size_t)size, sizeof(double));
        for (int i = 0; i < size; i++) {
            yi[nblock].indices[i] = xi[counter + i].index;
            yi[nblock].weights[i] = xi[counter + i].weight;
        }
        counter += size;
        printf("nblock %4d value %4.1f size %4d counter %4d", nblock,
               yi[nblock].value, yi[nblock].size, counter);
        printf("\nindices");
        for (int i = 0; i < size; i++) {
            printf("%4d", yi[nblock].indices[i] + 1);
        }
        printf("\nweights");
        for (int i = 0; i < size; i++) {
            printf("%4.1f", yi[nblock].weights[i]);
        }
        printf("\n");
        if (counter == n) {
            break;
        } else {
            nblock++;
        }
    }
    yi = (block *)realloc(yi, (size_t)nblock * sizeof(block));
    free(xi);
    return;
}

/*
double x1[10] = {3.0, 1.0, 1.0, 5.0, 1.0, 5.0, 1.0, 2.0, 5.0, 2.0};
double x2[10] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
double x3[10] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
double ww[10] = {1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0};
int n = 10;

int main() {
    block *yi = NULL;
    int nblock = 1;
    (void)smacofBlockSort(x1, ww, n, nblock, yi);
    printf("**********************************\n");
    (void)smacofBlockSort(x2, ww, n, nblock, yi);
    printf("**********************************\n");
    (void)smacofBlockSort(x3, ww, n, nblock, yi);
    free(yi);
}
*/
