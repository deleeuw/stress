#include <gsl/gsl_bspline.h>

void BSPLINE(double *x, int *order, int *nbreak, double *brpnts,
             double *results) {
    gsl_bspline_workspace *my_workspace =
        gsl_bspline_alloc((size_t)*order, (size_t)*nbreak);
    size_t ncoefs = gsl_bspline_ncoeffs(my_workspace);
    gsl_vector *values = gsl_vector_calloc(ncoefs);
    gsl_vector *breaks = gsl_vector_calloc((size_t)*nbreak);
    for (int i = 0; i < *nbreak; i++)
        gsl_vector_set(breaks, (size_t)i, brpnts[i]);
    (void)gsl_bspline_knots(breaks, my_workspace);
    (void)gsl_bspline_eval(*x, values, my_workspace);
    for (int i = 0; i < ncoefs; i++) results[i] = (values->data)[i];
    gsl_bspline_free(my_workspace);
    gsl_vector_free(values);
    gsl_vector_free(breaks);
}
