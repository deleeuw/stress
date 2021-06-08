#define TRUE 1

void bisplineC(const int *, const double *, const int *, const int *,
              const double *, int *, double *, double *, int *);
void bisplineBasisC(const int *, const double *, const int *,
                  const int *, const double *, const int *,
                  double *, double *);

inline int INDEX(const int);
inline int IMIN(const int, const int);
inline int IMAX(const int, const int);
inline int MINDEX(const int, const int, const int);

inline int INDEX(const int i) { return i - 1; }

inline int MINDEX(const int i, const int j, const int n) {
    return (i - 1) + (j - 1) * n;
}

inline int IMIN(const int a, const int b) {
    if (a > b) return b;
    return a;
}

inline int IMAX(const int a, const int b) {
    if (a < b) return b;
    return a;
}
