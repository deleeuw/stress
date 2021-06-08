#ifndef ARRAY_H
#define ARRAY_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define CTYPE double

typedef void (*ENCODE)(
    const int *, int *, const int *,
    const int *); // encode gives R array indices from C memory location
typedef void (*DECODE)(
    const int *, int *, const int *,
    const int *); // decode gives C memory location from R array indices

typedef struct {
  CTYPE *content;
  int *length;
  int *rank;
  int *shape;
  ENCODE encode;
  DECODE decode;
} ARRAY;


typedef ARRAY GVECTR;
typedef ARRAY GMATRX;
typedef ARRAY SMATRX;
typedef ARRAY GARRAY;
typedef ARRAY SARRAY;

/* protoypes */

extern int binCoef(const int, const int);

extern void isort(int *, const int *);

extern int icmp(const void *, const void *);

extern inline int IMIN(const int, const int);

extern void GVENCODE(const int *, int *, const int *, const int *);

extern void GVDECODE(const int *, int *, const int *, const int *);

extern void GMDECODE(const int *, int *, const int *, const int *);

extern void GMENCODE(const int *, int *, const int *, const int *);

extern void GADECODE(const int *, int *, const int *, const int *);

extern void GAENCODE(const int *, int *, const int *, const int *);

extern void SMDECODE(const int *, int *, const int *, const int *);

extern void SMENCODE(const int *, int *, const int *, const int *);

extern void SADECODE(const int *, int *, const int *, const int *);

extern void SAENCODE(const int *, int *, const int *, const int *);

#endif /* ARRAY_H */
