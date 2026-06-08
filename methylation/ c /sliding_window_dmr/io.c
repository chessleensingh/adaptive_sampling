#include "io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

RegionArray read_bed(const char *path) {
    RegionArray ra = {NULL, 0};
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(1); }

    size_t cap = 256;
    ra.data = malloc(cap * sizeof(Region));
    if (!ra.data) { perror("malloc"); exit(1); }

    char line[256];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        Region r;
        if (sscanf(line, "%31s %d %d", r.chr, &r.start, &r.end) != 3) continue;
        if (ra.n == cap) {
            cap *= 2;
            ra.data = realloc(ra.data, cap * sizeof(Region));
            if (!ra.data) { perror("realloc"); exit(1); }
        }
        ra.data[ra.n++] = r;
    }
    fclose(f);
    return ra;
}

SiteArray read_dss(const char *path) {
    SiteArray sa = {NULL, 0, 0};
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(1); }

    size_t cap = 1 << 20;   /* 1M sites initial */
    sa.data = malloc(cap * sizeof(Site));
    if (!sa.data) { perror("malloc"); exit(1); }
    sa.cap = cap;

    char line[256];
    fgets(line, sizeof(line), f);  /* skip header */

    while (fgets(line, sizeof(line), f)) {
        Site s;
        if (sscanf(line, "%31s %d %d %d", s.chr, &s.pos, &s.N, &s.X) != 4)
            continue;
        if (sa.n == sa.cap) {
            sa.cap *= 2;
            sa.data = realloc(sa.data, sa.cap * sizeof(Site));
            if (!sa.data) { perror("realloc"); exit(1); }
        }
        sa.data[sa.n++] = s;
    }
    fclose(f);
    return sa;
}

SiteArray read_dss_filtered(const char *path, const RegionArray *ra) {
    SiteArray sa = {NULL, 0, 0};
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(1); }

    size_t cap = 1 << 15;  /* 32k initial — we keep very few sites */
    sa.data = malloc(cap * sizeof(Site));
    if (!sa.data) { perror("malloc"); exit(1); }
    sa.cap = cap;

    char line[256];
    fgets(line, sizeof(line), f);  /* skip header */

    size_t ri = 0;  /* current region index (both files sorted by chr+pos) */

    while (fgets(line, sizeof(line), f)) {
        if (ri >= ra->n) break;  /* exhausted all regions */

        Site s;
        if (sscanf(line, "%31s %d %d %d", s.chr, &s.pos, &s.N, &s.X) != 4)
            continue;

        /* Advance past regions that end before this site */
        while (ri < ra->n) {
            int cmp = strcmp(ra->data[ri].chr, s.chr);
            if (cmp > 0) break;                              /* region chr > site chr */
            if (cmp == 0 && ra->data[ri].end > s.pos) break; /* region may include site */
            ri++;
        }
        if (ri >= ra->n) break;

        /* Check if site falls within the current region */
        if (strcmp(ra->data[ri].chr, s.chr) == 0 &&
            s.pos >= ra->data[ri].start &&
            s.pos <  ra->data[ri].end) {
            if (sa.n == sa.cap) {
                sa.cap *= 2;
                sa.data = realloc(sa.data, sa.cap * sizeof(Site));
                if (!sa.data) { perror("realloc"); exit(1); }
            }
            sa.data[sa.n++] = s;
        }
    }
    fclose(f);
    return sa;
}

void free_region_array(RegionArray *ra) { free(ra->data); ra->data = NULL; ra->n = 0; }
void free_site_array(SiteArray *sa)     { free(sa->data); sa->data = NULL; sa->n = sa->cap = 0; }

size_t bisect_left(const SiteArray *sa, const char *chr,
                   size_t lo, size_t hi, int32_t target) {
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        int cmp = strcmp(sa->data[mid].chr, chr);
        if (cmp < 0 || (cmp == 0 && sa->data[mid].pos < target))
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo;
}
