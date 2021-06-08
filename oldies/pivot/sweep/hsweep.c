void gsweep_(double *, double *, int *, int *, int *, int *, double *, int *);

void
hsweep(double *s, double *trian, int *n, int *m, int *ind, int *np, double *eps) {
	int ifault = 0, l = 0;
	for (int i = 0; i < *np; i++) {
		gsweep_ (s, trian, ind + i, &l, m, n, eps, &ifault);
	}
}