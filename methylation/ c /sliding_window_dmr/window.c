#include "window.h"
#include <string.h>

void agg_window(const SiteArray *sa, size_t lo, size_t hi,
                const char *chr, int32_t win_start,
                int32_t *x_out, int32_t *n_out) {
    int32_t win_end = win_start + WIN_SIZE;
    *x_out = 0; *n_out = 0;

    size_t i = bisect_left(sa, chr, lo, hi, win_start);
    for (; i < hi; i++) {
        if (sa->data[i].pos >= win_end) break;
        *x_out += sa->data[i].X;
        *n_out += sa->data[i].N;
    }
}

size_t slide_region(const SiteArray samples[N_SAMPLES],
                    const Region   *reg,
                    WinAgg         *aggs,
                    size_t          aggs_cap) {
    int32_t rlen = reg->end - reg->start;
    size_t  n_win = n_windows_for_len(rlen);
    if (n_win == 0 || n_win > aggs_cap) return 0;

    size_t lo[N_SAMPLES], hi[N_SAMPLES];
    for (int s = 0; s < N_SAMPLES; s++) {
        lo[s] = bisect_left(&samples[s], reg->chr, 0, samples[s].n, reg->start);
        hi[s] = bisect_left(&samples[s], reg->chr, lo[s], samples[s].n, reg->end);
    }

    for (size_t w = 0; w < n_win; w++) {
        int32_t ws = reg->start + (int32_t)(w * WIN_STRIDE);
        WinAgg *a  = &aggs[w];
        strncpy(a->chr, reg->chr, MAX_CHR_LEN - 1);
        a->chr[MAX_CHR_LEN - 1] = '\0';
        a->start = ws;
        a->end   = ws + WIN_SIZE;
        for (int s = 0; s < N_SAMPLES; s++) {
            agg_window(&samples[s], lo[s], hi[s],
                       reg->chr, ws, &a->x[s], &a->n[s]);
        }
    }
    return n_win;
}
