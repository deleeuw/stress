#include <stdio.h>
#include <stdlib.h>

void for_sub_(int *, double *, double *, int *);

int 
main()
  {
    int i=1, ny=3;
    double x=3.14159;
    double y[]={1.1, 2.2, 3.3};
    printf("C calling Fortran subroutine, passing\n");
    printf("i=%d, x=%g\n", i,x);
    printf("y[]=%g,%g,%g\n",y[0],y[1],y[2]);
    for_sub_(&i, &x, y, &ny);
    return EXIT_SUCCESS;
  }
