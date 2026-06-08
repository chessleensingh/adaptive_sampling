#pragma once
#include <stddef.h>

/*
 * Benjamini-Hochberg FDR correction.
 * pvals[n] — input p-values (not modified).
 * qvals[n] — output adjusted p-values (caller allocates, same length).
 * qvals[i] corresponds to pvals[i].
 */
void bh_fdr(const double *pvals, size_t n, double *qvals);
