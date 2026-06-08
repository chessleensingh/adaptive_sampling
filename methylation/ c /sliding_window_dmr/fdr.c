#include "fdr.h"
#include <stdlib.h>

typedef struct { double p; size_t idx; } PIdx;

static int cmp_pidx(const void *a, const void *b) {
    double da = ((const PIdx *)a)->p;
    double db = ((const PIdx *)b)->p;
    return (da > db) - (da < db);
}

void bh_fdr(const double *pvals, size_t n, double *qvals) {
    if (n == 0) return;

    PIdx *tmp = malloc(n * sizeof(PIdx));
    double *q_sorted = malloc(n * sizeof(double));
    if (!tmp || !q_sorted) { free(tmp); free(q_sorted); return; }

    for (size_t i = 0; i < n; i++) { tmp[i].p = pvals[i]; tmp[i].idx = i; }
    qsort(tmp, n, sizeof(PIdx), cmp_pidx);

    for (size_t r = 0; r < n; r++)
        q_sorted[r] = tmp[r].p * (double)n / (double)(r + 1);

    for (size_t r = n - 1; r > 0; r--)
        if (q_sorted[r - 1] > q_sorted[r]) q_sorted[r - 1] = q_sorted[r];

    for (size_t r = 0; r < n; r++)
        if (q_sorted[r] > 1.0) q_sorted[r] = 1.0;

    for (size_t r = 0; r < n; r++)
        qvals[tmp[r].idx] = q_sorted[r];

    free(tmp);
    free(q_sorted);
}
