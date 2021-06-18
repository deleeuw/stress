#include <stdio.h>
#include <stdlib.h>

int main () {
	double b[4] = {1.0, 2.0, 3.0, 4.0};
	double *a = b;
	printf("%10.4f %10.4f\n", a[3], b[3]);
	double *c = (double *)calloc((size_t) 4, sizeof(double));
	c = b;
    printf("%10.4f %10.4f\n", b[3], c[3]);
    free(c);
	printf("%10.4f %10.4f\n", a[3], b[3]);
}