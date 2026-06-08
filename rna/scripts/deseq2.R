library(DESeq2)

base_dir <- "/home/sachleen.fnu/Projects/04282026_RNA_SEQ_Psomagen_Sachini"
out_dir  <- file.path(base_dir, "deseq2")
dir.create(out_dir, showWarnings=FALSE)

# Load counts
counts_raw <- read.table(file.path(base_dir, "counts/all_counts.txt"),
                         header=TRUE, sep="\t", comment.char="#", row.names=1)
counts <- counts_raw[, 6:ncol(counts_raw)]
colnames(counts) <- sub(".*\\.aligned\\.(.+)\\.Aligned\\..*", "\\1", colnames(counts))

# Load sample info
sample_info <- read.csv(file.path(base_dir, "sample_info.csv"), row.names=1)
sample_info$time      <- factor(sample_info$time)
sample_info$treatment <- factor(sample_info$treatment, levels=c("Mock","Race1","Race25"))
sample_info$group     <- factor(sample_info$group)

counts <- counts[, rownames(sample_info)]
counts <- counts[rowSums(counts) >= 10, ]

# Build and run DESeq2
dds <- DESeqDataSetFromMatrix(countData=counts, colData=sample_info, design=~group)
dds <- DESeq(dds)

# Helper: extract results, sort by padj, write CSV
save_results <- function(dds, contrast, label) {
  res <- results(dds, contrast=c("group", contrast[1], contrast[2]),
                 independentFiltering=TRUE)
  res <- res[order(res$padj, na.last=TRUE), ]
  out <- as.data.frame(res)
  out <- cbind(gene=rownames(out), out)
  fname <- file.path(out_dir, paste0(label, ".csv"))
  write.csv(out, fname, row.names=FALSE)
  sig <- sum(res$padj < 0.05 & abs(res$log2FoldChange) >= 1, na.rm=TRUE)
  cat(sprintf("  %-40s  %d DEGs (padj<0.05, |LFC|>=1)\n", label, sig))
}

cat("\n=== DESeq2 results ===\n")
save_results(dds, c("12h_Race1",  "12h_Mock"),  "12h_Race1_vs_Mock")
save_results(dds, c("12h_Race25", "12h_Mock"),  "12h_Race25_vs_Mock")
save_results(dds, c("24h_Race1",  "24h_Mock"),  "24h_Race1_vs_Mock")
save_results(dds, c("24h_Race25", "24h_Mock"),  "24h_Race25_vs_Mock")
save_results(dds, c("12h_Race25", "12h_Race1"), "12h_Race25_vs_Race1")
save_results(dds, c("24h_Race25", "24h_Race1"), "24h_Race25_vs_Race1")

cat("\nResults written to deseq2/\n")
