"""
Extract ONT bedMethyl sites for the 133 buffered target regions.

Input:  results/ont/prepared/barcode{01-04}_{ctx}.bed.gz  (genome-wide, bgzf)
        data/target_regions_buffered.bed

Output: results/ont/prepared_buffered/barcode{01-04}_{ctx}_targets.bed
            - plain text, all 18 bedMethyl cols, no header, only covered sites (col9 > 0)
        results/ont/dss_buffered/barcode{01-04}_{ctx}_dss.txt
            - DSS format: chr  pos(1-based)  N  X
"""

import gzip
from pathlib import Path
from collections import defaultdict

BASE_DIR      = Path("V:/Code/Methylation_analysis/04032026_methylation")
TARGET_BED    = BASE_DIR / "data" / "target_regions_buffered.bed"
PREPARED_DIR  = BASE_DIR / "results" / "ont" / "prepared"
OUT_TARGETS   = BASE_DIR / "results" / "ont" / "prepared_buffered"
OUT_DSS       = BASE_DIR / "results" / "ont" / "dss_buffered"

OUT_TARGETS.mkdir(parents=True, exist_ok=True)
OUT_DSS.mkdir(parents=True, exist_ok=True)

BARCODES  = ["barcode01", "barcode02", "barcode03", "barcode04"]
CONTEXTS  = ["CG", "CHG", "CHH"]


def load_targets(bed_path: Path) -> dict:
    regions = defaultdict(list)
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            chrom, t_start, t_end = parts[0], int(parts[1]), int(parts[2])
            regions[chrom].append((t_start, t_end))
    return dict(regions)


def site_in_targets(chrom: str, start: int, regions: dict) -> bool:
    if chrom not in regions:
        return False
    return any(t_start <= start < t_end for t_start, t_end in regions[chrom])


def extract_barcode(barcode: str, context: str, regions: dict) -> tuple:
    gz_path   = PREPARED_DIR / f"{barcode}_{context}.bed.gz"
    tgt_path  = OUT_TARGETS / f"{barcode}_{context}_targets.bed"
    dss_path  = OUT_DSS / f"{barcode}_{context}_dss.txt"

    if not gz_path.exists():
        print(f"  SKIP: {gz_path.name} not found")
        return 0, 0

    n_sites = n_covered = 0

    with gzip.open(gz_path, "rt") as gz_in, \
         open(tgt_path, "w") as tgt_out, \
         open(dss_path, "w") as dss_out:

        dss_out.write("chr\tpos\tN\tX\n")

        for raw in gz_in:
            if not raw.strip() or raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            chrom = cols[0]
            start = int(cols[1])
            n_cov = int(cols[9])

            if not site_in_targets(chrom, start, regions):
                continue

            n_sites += 1
            if n_cov == 0:
                continue
            n_covered += 1

            tgt_out.write(raw if raw.endswith("\n") else raw + "\n")

            pos = start + 1
            n_mod = int(cols[11])
            dss_out.write(f"{chrom}\t{pos}\t{n_cov}\t{n_mod}\n")

    return n_sites, n_covered


targets = load_targets(TARGET_BED)
print(f"Loaded {sum(len(v) for v in targets.values())} regions across {len(targets)} chromosomes")

for ctx in CONTEXTS:
    print(f"\n--- Context: {ctx} ---")
    for bc in BARCODES:
        n_sites, n_cov = extract_barcode(bc, ctx, targets)
        print(f"  {bc}: {n_sites:,} sites in regions, {n_cov:,} covered (written to DSS)")

print("\nDone.")
