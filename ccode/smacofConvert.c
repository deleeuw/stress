#include "../include/smacof.h"

void smacofEncode(const int *ip, const int *jp, const int *np, int *kp) {
    int i = *ip, j = *jp, n = *np;
    *kp = i + (j - 1) * n - j * (j + 1) / 2;
    return;
}

void smacofDecode(const int *kp, const int *np, int *ip, int *jp) {
    int j = 1, m = 1, k = *kp, n = *np;
    while (k >= ((j * n) - m + 1)) {
        j += 1;
        m += j;
    }
    *ip = k - (j - 1) * n + m;
    *jp = j;
    return;
}








/*
int main(void) {
    int i = 0, j = 0, k = 0, n = 6;
    printf("ENCODE\n\n");
    for (j = 1; j < n; j++) {
        for (i = (j + 1); i <= n; i++) {
            (void)dencode(&i, &j, &n, &k);
            printf(" %4d ", k);
        }
    }
    printf("\n\n");
    printf("DECODE\n\n");
    for (int k = 1; k <= n * (n - 1) / 2; k++) {
        (void)ddecode(&k, &n, &i, &j);
        printf("%4d %4d\n", i, j);
    }
    return EXIT_SUCCESS;
}
*/

