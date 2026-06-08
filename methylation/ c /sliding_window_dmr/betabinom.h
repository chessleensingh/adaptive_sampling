#pragma once
#include "types.h"

/*
 * Beta-binomial log-likelihood for n samples with fixed rho.
 * X[i] = methylated reads, N[i] = total reads for sample i.
 * Samples with N[i]==0 are skipped.
 */
double bb_loglik(double mu, double rho,
                 const int32_t *X, const int32_t *N, int n);

/*
 * MLE of mu (mean methylation) for given rho via golden-section search.
 * Returns NAN if all samples have N==0.
 */
double bb_mle_mu(double rho, const int32_t *X, const int32_t *N, int n);

/*
 * Likelihood ratio test statistic comparing mock (N_MOCK samples) vs
 * treatment (N_TREAT samples) with fixed global dispersion rho.
 * Returns 0 if the test cannot be run (insufficient data).
 */
double bb_lrt(double rho,
              const int32_t xm[N_MOCK], const int32_t nm[N_MOCK],
              const int32_t xt[N_TREAT], const int32_t nt[N_TREAT]);

/* p-value from chi-squared(1): erfc(sqrt(g/2)) */
double bb_pvalue(double g);

/*
 * Estimate global overdispersion rho via method-of-moments pooled
 * across all windows and both groups.
 * aggs: array of n_agg WinAgg structs (from slide_region).
 * Returns rho clamped to [1e-6, 0.5]; falls back to 0.01 if estimation fails.
 */
double estimate_dispersion(const WinAgg *aggs, size_t n_agg);
