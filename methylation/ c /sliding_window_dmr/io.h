#pragma once
#include "types.h"

RegionArray read_bed(const char *path);
SiteArray   read_dss(const char *path);

/*
 * Like read_dss but only keeps sites that fall within any region in ra.
 * Assumes both the DSS file and ra are sorted by chr then position.
 * Drastically reduces memory for large contexts like CHH.
 */
SiteArray   read_dss_filtered(const char *path, const RegionArray *ra);
void        free_region_array(RegionArray *ra);
void        free_site_array(SiteArray *sa);

/*
 * First index i in [lo, hi) where sa->data[i].chr == chr
 * and sa->data[i].pos >= target.  Returns hi if not found.
 */
size_t bisect_left(const SiteArray *sa, const char *chr,
                   size_t lo, size_t hi, int32_t target);
