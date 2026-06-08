#include "betabinom.h"
#include <math.h>
#include <string.h>

double bb_loglik(double mu, double rho,
                 const int32_t *X, const int32_t *N, int n) {
    if (mu <= 0.0) mu = 1e-9;
    if (mu >= 1.0) mu = 1.0 - 1e-9;
    double alpha = mu * (1.0 - rho) / rho;
    double beta  = (1.0 - mu) * (1.0 - rho) / rho;
    double ab    = alpha + beta;
    double ll    = 0.0;
    for (int i = 0; i < n; i++) {
        if (N[i] == 0) continue;
        ll += lgamma(X[i] + alpha)
            + lgamma(N[i] - X[i] + beta)
            - lgamma(N[i] + ab)
            - lgamma(alpha)
            - lgamma(beta)
            + lgamma(ab);
    }
    return ll;
}

typedef struct { double rho; const int32_t *X; const int32_t *N; int n; } BBCtx;

static double gs_f(double mu, const BBCtx *c) {
    return bb_loglik(mu, c->rho, c->X, c->N, c->n);
}

double bb_mle_mu(double rho, const int32_t *X, const int32_t *N, int n) {
    int total = 0;
    for (int i = 0; i < n; i++) total += N[i];
    if (total == 0) return (double)0.0 / 0.0; /* NAN */

    BBCtx ctx = {rho, X, N, n};
    const double phi = (sqrt(5.0) - 1.0) / 2.0;
    double a = 1e-9, b = 1.0 - 1e-9;
    double c = b - phi * (b - a);
    double d = a + phi * (b - a);
    for (int iter = 0; iter < 100 && (b - a) > 1e-8; iter++) {
        if (gs_f(c, &ctx) > gs_f(d, &ctx)) b = d;
        else                                a = c;
        c = b - phi * (b - a);
        d = a + phi * (b - a);
    }
    return (a + b) / 2.0;
}

double bb_lrt(double rho,
              const int32_t xm[N_MOCK], const int32_t nm[N_MOCK],
              const int32_t xt[N_TREAT], const int32_t nt[N_TREAT]) {
    int sum_m = 0, sum_t = 0;
    for (int i = 0; i < N_MOCK;  i++) sum_m += nm[i];
    for (int i = 0; i < N_TREAT; i++) sum_t += nt[i];
    if (sum_m < MIN_READS || sum_t < MIN_READS) return 0.0;

    int32_t x_all[N_MOCK + N_TREAT], n_all[N_MOCK + N_TREAT];
    memcpy(x_all,          xm, N_MOCK  * sizeof(int32_t));
    memcpy(x_all + N_MOCK, xt, N_TREAT * sizeof(int32_t));
    memcpy(n_all,          nm, N_MOCK  * sizeof(int32_t));
    memcpy(n_all + N_MOCK, nt, N_TREAT * sizeof(int32_t));

    double mu0  = bb_mle_mu(rho, x_all, n_all, N_MOCK + N_TREAT);
    double mu_m = bb_mle_mu(rho, xm, nm, N_MOCK);
    double mu_t = bb_mle_mu(rho, xt, nt, N_TREAT);

    if (mu0 != mu0 || mu_m != mu_m || mu_t != mu_t) return 0.0; /* NAN check */

    double ll_null = bb_loglik(mu0,  rho, x_all, n_all, N_MOCK + N_TREAT);
    double ll_alt  = bb_loglik(mu_m, rho, xm, nm, N_MOCK)
                   + bb_loglik(mu_t, rho, xt, nt, N_TREAT);

    double g = 2.0 * (ll_alt - ll_null);
    return g > 0.0 ? g : 0.0;
}

double bb_pvalue(double g) {
    if (g <= 0.0) return 1.0;
    return erfc(sqrt(g / 2.0));
}

double estimate_dispersion(const WinAgg *aggs, size_t n_agg) {
    double sum_var = 0.0, sum_mu_var = 0.0;
    size_t count = 0;

    for (size_t i = 0; i < n_agg; i++) {
        for (int g = 0; g < 2; g++) {
            int s0 = g * 2, s1 = g * 2 + 1;
            if (aggs[i].n[s0] == 0 || aggs[i].n[s1] == 0) continue;
            double p0 = (double)aggs[i].x[s0] / aggs[i].n[s0];
            double p1 = (double)aggs[i].x[s1] / aggs[i].n[s1];
            double mu = (p0 + p1) / 2.0;
            if (mu < 1e-6 || mu > 1.0 - 1e-6) continue;
            double within_var   = (p0 - p1) * (p0 - p1) / 2.0;
            double binom_noise  = mu * (1.0 - mu) / 2.0
                                * (1.0 / aggs[i].n[s0] + 1.0 / aggs[i].n[s1]);
            sum_var    += within_var - binom_noise;
            sum_mu_var += mu * (1.0 - mu);
            count++;
        }
    }

    if (count == 0 || sum_mu_var < 1e-10) return 0.01;
    double rho = sum_var / sum_mu_var;
    if (rho < 1e-6) rho = 1e-6;
    if (rho > 0.5)  rho = 0.5;
    return rho;
}
