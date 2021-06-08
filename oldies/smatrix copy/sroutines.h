#ifndef SROUTINES_H
#define SROUTINES_H

extern void scholesky(const int *, double *, bool *, bool *, double *,
                      const double *);
extern void pcholesky(const int *, double *, double *, int *, int *, bool *,
                      const double *);

/* defined in jacobi.c */

extern void smjacobi(const int *, const int *, double *, double *, double *,
                     double *, int *, const int *, const double *, const bool *,
                     const bool *);
extern void sjacobi(const int *, double *, double *, double *, double *, int *,
                    const int *, const double *, const bool *, const bool *);

/* defined in sweep.c */

extern void ssweep(const int *, double *, double *, int *, bool *,
                   const double *);
extern void psweep(const int *, double *, double *, int *, int *, bool *,
                   const double *);

/* defined in invtri.c */

extern void invtri(const int *, double *, const double *, bool *);

/* defined in utility.c */

extern void primat(const int *, const int *, const int *, const int *,
                   const double *);
extern void pritru(const int *, const int *, const int *, const double *);
extern void prisru(const int *, const int *, const int *, const int *,
                   const double *);
extern void trimat(const int *, const double *, double *);
extern void mattri(const int *, const double *, double *);

#endif /* SROUTINES_H */