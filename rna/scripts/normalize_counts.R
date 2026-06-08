library(DESeq2)

base_dir <- "/home/sachleen.fnu/Projects/04282026_RNA_SEQ_Psomagen_Sachini"

# Load counts (skip featureCounts comment line, keep gene rows only)
counts_raw <- read.table(file.path(base_dir, "counts/all_counts.txt"),
                         header=TRUE, sep="\t", comment.char="#", row.names=1)

# Drop annotation columns (Chr, Start, End, Strand, Length), keep count columns
counts <- counts_raw[, 6:ncol(counts_raw)]

# Clean up column names to just sample names
colnames(counts) <- sub(".*\\.aligned\\.(.+)\\.Aligned\\..*", "\\1", colnames(counts))

# Load sample info
sample_info <- read.csv(file.path(base_dir, "sample_info.csv"), row.names=1)
sample_info$time      <- factor(sample_info$time)
sample_info$treatment <- factor(sample_info$treatment, levels=c("Mock","Race1","Race25"))

# Make sure column order matches
counts <- counts[, rownames(sample_info)]

# Remove low-count genes
counts <- counts[rowSums(counts) >= 10, ]

# Build DESeq2 object and estimate size factors
dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=sample_info,
                              design=~group)
dds <- estimateSizeFactors(dds)

# Normalized counts
norm_counts <- counts(dds, normalized=TRUE)
write.csv(norm_counts, file.path(base_dir, "counts/normalized_counts.csv"))
cat("Normalized counts written.\n")

# Glyma.04G121700 expression
gene_id <- grep("Glyma.04G121700", rownames(norm_counts), value=TRUE)
gene_expr <- data.frame(
  sample    = colnames(norm_counts),
  time      = sample_info$time,
  treatment = sample_info$treatment,
  normalized_counts = round(norm_counts[gene_id, ], 2)
)
write.csv(gene_expr, file.path(base_dir, "Glyma.04G121700_normalized.csv"), row.names=FALSE)
cat("Gene expression table written.\n")
print(gene_expr)
