#pragma once
#include "types.h"
#include "io.h"

/*
 * Aggregate methylated (x_out) and total (n_out) reads from sa->data in
 * [win_start, win_start + WIN_SIZE) on chromosome chr.
 * lo/hi bound the search range within sa (pre-computed per region).
 */
void agg_window(const SiteArray *sa, size_t lo, size_t hi,
                const char *chr, int32_t win_start,
                int32_t *x_out, int32_t *n_out);

/*
 * Slide WIN_SIZE/WIN_STRIDE windows across region reg for all N_SAMPLES
 * site arrays.  Fills aggs[0..n_wins-1]; caller must pre-allocate aggs.
 * Returns number of windows written (== n_wins).
 */
size_t slide_region(const SiteArray samples[N_SAMPLES],
                    const Region   *reg,
                    WinAgg         *aggs,
                    size_t          aggs_cap);

/* Number of windows for a region of given length */
static inline size_t n_windows_for_len(int32_t len) {
    if (len < WIN_SIZE) return 0;
    return (size_t)((len - WIN_SIZE) / WIN_STRIDE) + 1;
}
