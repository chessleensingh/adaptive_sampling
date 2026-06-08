# Analysis Source Code

Source code for the soybean (*Glycine max* Wm82.a4.v1) epigenomic and transcriptomic analysis accompanying:

> **ONT Adaptive Sampling Reveals Promoter-Associated Differential Methylation Independent of Transcriptional Response During Early *Phytophthora sojae* Infection**
> Sachleen Singh, Sachini Wijeratne — Arkansas State University

---

## Repository Structure

```
src/
├── methylation/
│   ├── c/
│   │   └── sliding_window_dmr/     # C implementation of the sliding-window DMR caller
│   └── ont_to_dss/
│       └── extract_ont_buffered.py         # ONT bedMethyl → DSS format conversion
└── rna/
    ├── deseq2_pairwise.R           # Pairwise DESeq2 contrasts (all timepoints)
    ├── deseq2_timepoint.R          # Timepoint-level DESeq2 analysis
    ├── pca_analysis.R              # PCA of normalised counts
    └── scripts/                   # STAR alignment pipeline (HPC SLURM)
        ├── 00_setup.sh
        ├── 01_fastqc_raw.sh
        ├── 02_trim_galore.sh
        ├── 03_fastqc_trimmed.sh
        ├── 04_star_align.sh
        ├── 05_featurecounts.sh
        ├── 06_multiqc.sh
        ├── deseq2.R
        └── normalize_counts.R
```

---

## Methylation Analysis

### Sliding-Window DMR Caller (`methylation/c/sliding_window_dmr/`)

A multi-threaded C implementation of a beta-binomial likelihood ratio test applied across 100 bp sliding windows (10 bp stride) over the 133 ONT RRMS target regions. Significant windows are called at FDR < 0.05 (Benjamini–Hochberg).

**Dependencies:** GCC with OpenMP, standard C library (`-lm`)

**Build:**
```bash
cd methylation/c/sliding_window_dmr
gcc -O2 -fopenmp main.c betabinom.c fdr.c io.c window.c -lm -o sliding_window_dmr
```

**Usage:**
```bash
./sliding_window_dmr <CG|CHG|CHH> <dss_input_dir> <regions.bed> <output.tsv> [n_threads]
```

**Module overview:**

| File | Description |
|------|-------------|
| `main.c` | Entry point, argument parsing, per-context orchestration |
| `types.h` | Shared struct definitions (`SiteData`, `Window`, `RegionResult`) |
| `io.c / io.h` | DSS-format input parsing (`chr pos N X`) |
| `window.c / window.h` | Sliding window aggregation and replicate averaging |
| `betabinom.c / betabinom.h` | Beta-binomial dispersion estimation and LRT statistic |
| `fdr.c / fdr.h` | Benjamini–Hochberg FDR correction |

---

### ONT to DSS Conversion (`methylation/ont_to_dss/`)

| Script | Purpose |
|--------|---------|
| `extract_ont_buffered.py` | Subsets genome-wide ONT bedMethyl files to the 133 buffered target regions and converts to DSS input format (`chr pos N X`) for downstream differential methylation analysis |

**Dependencies:** Python ≥ 3.8, standard library only (`gzip`, `pathlib`)

---

## RNA-seq Analysis

### Preprocessing Pipeline (`rna/scripts/`)

Shell scripts designed for execution on an HPC cluster with SLURM. Run sequentially or via `run_incremental.sh`.

| Script | Tool | Description |
|--------|------|-------------|
| `00_setup.sh` | — | Creates directory structure and sets paths |
| `01_fastqc_raw.sh` | FastQC | QC on raw FASTQ files |
| `02_trim_galore.sh` | Trim Galore | Adapter trimming and quality filtering |
| `03_fastqc_trimmed.sh` | FastQC | QC on trimmed reads |
| `04_star_align.sh` | STAR | Alignment to *Glycine max* Wm82.a4.v1 |
| `05_featurecounts.sh` | featureCounts | Read quantification against Wm82.a4.v1 GFF3 |
| `06_multiqc.sh` | MultiQC | Aggregated QC report |
| `deseq2.R` | DESeq2 | Differential expression analysis |
| `normalize_counts.R` | DESeq2/edgeR | Count normalisation (VST/TMM) |

### Differential Expression (`rna/`)

| Script | Description |
|--------|-------------|
| `deseq2_pairwise.R` | All pairwise DESeq2 contrasts (Race 1 vs. mock, Race 25 vs. mock, Race 25 vs. Race 1) at 12 hpi and 24 hpi |
| `deseq2_timepoint.R` | Timepoint-level DESeq2 analysis |
| `pca_analysis.R` | Principal component analysis of VST-normalised counts with sample grouping |

**Dependencies:** R ≥ 4.2, `DESeq2`, `ggplot2`, `dplyr`, `data.table`

---

## Data

Raw ONT sequencing data and RNA-seq FASTQ files are deposited at NCBI under BioProject **[accession pending]**.

Reference genome: *Glycine max* Wm82.a4.v1 (Phytozome v13 / JGI)

---

## Contact

Asela Wijeratne — awijeratne@astate.edu
