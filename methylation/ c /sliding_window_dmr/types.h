#pragma once
#include <stdint.h>
#include <stddef.h>

#define MAX_CHR_LEN  32
#define N_SAMPLES     4   /* M1, M3, R1, R3 */
#define N_MOCK        2   /* samples 0,1 */
#define N_TREAT       2   /* samples 2,3 */
#define WIN_SIZE    100
#define WIN_STRIDE   10
#define MIN_READS     1   /* skip window if either group has 0 total reads */

typedef struct {
    char    chr[MAX_CHR_LEN];
    int32_t pos;
    int32_t N;   /* total reads */
    int32_t X;   /* methylated reads */
} Site;

typedef struct {
    Site   *data;
    size_t  n;
    size_t  cap;
} SiteArray;

typedef struct {
    char    chr[MAX_CHR_LEN];
    int32_t start;
    int32_t end;
} Region;

typedef struct {
    Region *data;
    size_t  n;
} RegionArray;

/* Per-window raw aggregates across all 4 samples */
typedef struct {
    char    chr[MAX_CHR_LEN];
    int32_t start;
    int32_t end;
    int32_t x[N_SAMPLES];  /* methylated count per sample */
    int32_t n[N_SAMPLES];  /* total reads per sample */
} WinAgg;

/* Final per-window result after LRT + FDR */
typedef struct {
    char    chr[MAX_CHR_LEN];
    int32_t start;
    int32_t end;
    int32_t n_mock;        /* total reads in mock (sum M1+M3) */
    int32_t n_treat;       /* total reads in treatment (sum R1+R3) */
    double  meth_mock;
    double  meth_treat;
    double  lrt_stat;
    double  p_value;
    double  fdr;
    int     direction;     /* 1 = hyper (treat > mock), -1 = hypo */
} WinResult;
