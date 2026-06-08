#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "types.h"
#include "io.h"
#include "window.h"
#include "betabinom.h"
#include "fdr.h"

static const char *SAMPLE_NAMES[N_SAMPLES] = {"M1", "M3", "R1", "R3"};

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s <context> <input_dir> <regions_bed> <output_tsv> [n_threads]\n"
        "  context     : CG | CHG | CHH\n"
        "  input_dir   : directory containing {M1,M3,R1,R3}_{context}_dss.txt\n"
        "  regions_bed : target_regions_buffered.bed\n"
        "  output_tsv  : path for output\n"
        "  n_threads   : number of OpenMP threads (default 12, ignored if no OpenMP)\n",
        prog);
    exit(1);
}

int main(int argc, char *argv[]) {
    if (argc < 5) usage(argv[0]);

    const char *ctx      = argv[1];
    const char *indir    = argv[2];
    const char *bed_path = argv[3];
    const char *out_path = argv[4];
    int n_threads = (argc >= 6) ? atoi(argv[5]) : 12;

#ifdef _OPENMP
    omp_set_num_threads(n_threads);
    fprintf(stderr, "OpenMP enabled: %d threads\n", n_threads);
#else
    (void)n_threads;
    fprintf(stderr, "OpenMP not available: running serially\n");
#endif

    /* 1. Load regions first (needed for filtered site loading) */
    RegionArray regions = read_bed(bed_path);
    fprintf(stderr, "Regions: %zu\n", regions.n);

    /* 2. Load DSS files — filtered to target regions only (memory-efficient) */
    SiteArray samples[N_SAMPLES];
    for (int s = 0; s < N_SAMPLES; s++) {
        char path[512];
        snprintf(path, sizeof(path), "%s/%s_%s_dss.txt", indir, SAMPLE_NAMES[s], ctx);
        fprintf(stderr, "Loading %s ...\n", path);
        samples[s] = read_dss_filtered(path, &regions);
        fprintf(stderr, "  %zu sites in target regions\n", samples[s].n);
    }

    /* 3. Pre-compute per-region window counts and global offset */
    size_t *offsets = malloc((regions.n + 1) * sizeof(size_t));
    if (!offsets) { perror("malloc offsets"); exit(1); }
    offsets[0] = 0;
    for (size_t r = 0; r < regions.n; r++) {
        int32_t len = regions.data[r].end - regions.data[r].start;
        offsets[r + 1] = offsets[r] + n_windows_for_len(len);
    }
    size_t total_win = offsets[regions.n];
    fprintf(stderr, "Total windows: %zu\n", total_win);

    WinAgg    *aggs    = malloc(total_win * sizeof(WinAgg));
    WinResult *results = malloc(total_win * sizeof(WinResult));
    double    *pvals   = malloc(total_win * sizeof(double));
    double    *qvals   = malloc(total_win * sizeof(double));
    if (!aggs || !results || !pvals || !qvals) { perror("malloc arrays"); exit(1); }

    /* 4. Pass 1: aggregate windows (parallel over regions) */
    fprintf(stderr, "Pass 1: aggregating windows...\n");
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int r = 0; r < (int)regions.n; r++) {
        slide_region(samples, &regions.data[r],
                     aggs + offsets[r], offsets[r + 1] - offsets[r]);
    }

    /* 5. Estimate global dispersion */
    double rho = estimate_dispersion(aggs, total_win);
    fprintf(stderr, "Estimated rho (dispersion): %.6f\n", rho);

    /* 6. Pass 2: LRT (parallel over windows) */
    fprintf(stderr, "Pass 2: computing LRT...\n");
    size_t covered = 0;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) reduction(+:covered)
#endif
    for (int w = 0; w < (int)total_win; w++) {
        const WinAgg *a = &aggs[w];
        double g    = bb_lrt(rho, a->x, a->n, a->x + N_MOCK, a->n + N_MOCK);
        double pval = bb_pvalue(g);

        WinResult *res = &results[w];
        strncpy(res->chr, a->chr, MAX_CHR_LEN - 1);
        res->chr[MAX_CHR_LEN - 1] = '\0';
        res->start   = a->start;
        res->end     = a->end;
        res->n_mock  = a->n[0] + a->n[1];
        res->n_treat = a->n[2] + a->n[3];
        res->meth_mock  = res->n_mock  > 0
                        ? (double)(a->x[0] + a->x[1]) / res->n_mock  : 0.0;
        res->meth_treat = res->n_treat > 0
                        ? (double)(a->x[2] + a->x[3]) / res->n_treat : 0.0;
        res->lrt_stat  = g;
        res->p_value   = pval;
        res->direction = (res->meth_treat > res->meth_mock) ? 1 : -1;
        pvals[w] = pval;
        if (res->n_mock > 0 && res->n_treat > 0) covered++;
    }
    fprintf(stderr, "Windows with coverage in both groups: %zu / %zu\n",
            covered, total_win);

    /* 7. BH FDR — only on covered windows to avoid dilution from zero-coverage windows */
    fprintf(stderr, "Applying BH FDR (covered windows only)...\n");
    double *cov_pvals = malloc(covered * sizeof(double));
    size_t *cov_idx   = malloc(covered * sizeof(size_t));
    double *cov_qvals = malloc(covered * sizeof(double));
    if (!cov_pvals || !cov_idx || !cov_qvals) { perror("malloc fdr"); exit(1); }

    size_t ci = 0;
    for (size_t w = 0; w < total_win; w++) {
        if (results[w].n_mock > 0 && results[w].n_treat > 0) {
            cov_pvals[ci] = pvals[w];
            cov_idx[ci]   = w;
            ci++;
        }
    }

    bh_fdr(cov_pvals, covered, cov_qvals);

    for (size_t w = 0; w < total_win; w++) results[w].fdr = 1.0;
    size_t sig05 = 0;
    for (size_t j = 0; j < covered; j++) {
        results[cov_idx[j]].fdr = cov_qvals[j];
        if (cov_qvals[j] < 0.05) sig05++;
    }
    free(cov_pvals); free(cov_idx); free(cov_qvals);
    fprintf(stderr, "Windows with FDR < 0.05: %zu\n", sig05);

    /* 8. Write output */
    fprintf(stderr, "Writing %s ...\n", out_path);
    FILE *out = fopen(out_path, "w");
    if (!out) { perror(out_path); exit(1); }

    fprintf(out, "chr\tstart\tend\tn_mock\tn_treat\t"
                 "meth_mock\tmeth_treat\tlrt_stat\tp_value\tfdr\tdirection\n");

    size_t written = 0;
    for (size_t w = 0; w < total_win; w++) {
        const WinResult *res = &results[w];
        if (res->n_mock == 0 && res->n_treat == 0) continue;
        fprintf(out, "%s\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6g\t%.6g\t%s\n",
                res->chr, res->start, res->end,
                res->n_mock, res->n_treat,
                res->meth_mock, res->meth_treat,
                res->lrt_stat, res->p_value, res->fdr,
                res->direction == 1 ? "hyper" : "hypo");
        written++;
    }
    fclose(out);
    fprintf(stderr, "Wrote %zu windows.\n", written);

    /* Cleanup */
    free(aggs); free(results); free(pvals); free(qvals); free(offsets);
    for (int s = 0; s < N_SAMPLES; s++) free_site_array(&samples[s]);
    free_region_array(&regions);

    fprintf(stderr, "Done.\n");
    return 0;
}
