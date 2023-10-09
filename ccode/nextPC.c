
void swap(int* x, int i, int j) {
    int temp;
    temp = x[i];
    x[i] = x[j];
    x[j] = temp;
}

void nextPermutation(int* x, int* nn) {
    int i, j, n = *nn;
    i = n - 1;
    while (x[i - 1] >= x[i]) i--;
    if (i == 0) return;
    j = n;
    while (x[j - 1] <= x[i - 1]) j--;
    swap(x, i - 1, j - 1);
    j = n;
    i++;
    while (i < j) {
        swap(x, i - 1, j - 1);
        j--;
        i++;
    }
}

void nextCombination(int* n, int* m, int* next) {
    int i, j, mm = *m - 1, nn = *n;
    for (i = mm; i >= 0; i--) {
        if (next[i] != nn - mm + i) {
            next[i]++;
            if (i < mm) {
                for (j = i + 1; j <= mm; j++) next[j] = next[j - 1] + 1;
            }
            return;
        }
    }
}
