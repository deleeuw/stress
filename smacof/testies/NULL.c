#include <stdio.h>
#include <stdlib.h>

void hownow(double *x, double *y);

int main() {
	double *y;
	*y = 1.00;
//	double *x = NULL;
    double *x = y;
	(void) hownow (x, y);
	printf ("%4.2f \n", *y);
	return EXIT_SUCCESS;	
}

void hownow(double *x, double *y) {
 if (x == NULL) {
 	*y = 0.0;
 } else {
 	*y *= 2.0;
 }
 return;
}