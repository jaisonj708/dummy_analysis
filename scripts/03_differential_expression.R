#!/usr/bin/env Rscript
# ==============================================================================
# Differential Expression Analysis with DESeq2
# ==============================================================================
# This script performs differential expression analysis:
#   1. Create DESeq2 dataset
#   2. Normalization
#   3. Dispersion estimation
#   4. Statistical testing
#   5. Multiple testing correction
#   6. Export results

cat("\n=================================================================\n")
cat("  Step 3: Differential Expression Analysis\n")
cat("=================================================================\n\n")

# Load configuration
source("config.R")
create_directories()

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
})

# Set random seed
set.seed(RANDOM_SEED)

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading filtered count data...\n")

# Load filtered count matrix
counts <- read_csv(file.path(DATA_DIR, "filtered_counts.csv"), show_col_types = FALSE)
count_matrix <- counts %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

# Load sample metadata
metadata <- read_csv(SAMPLE_METADATA_FILE, show_col_types = FALSE) %>%
  column_to_rownames("sample_id")

# Ensure metadata matches count matrix
metadata <- metadata[colnames(count_matrix), ]

# Convert condition to factor (required for DESeq2)
metadata$condition <- factor(metadata$condition, levels = c("Control", "Treated"))
metadata$batch <- factor(metadata$batch)

cat(sprintf("  - Genes: %d\n", nrow(count_matrix)))
cat(sprintf("  - Samples: %d\n", ncol(count_matrix)))
cat(sprintf("  - Conditions: %s\n", paste(unique(metadata$condition), collapse = ", ")))

# ==============================================================================
# CREATE DESEQ2 DATASET
# ==============================================================================

cat("\nCreating DESeq2 dataset...\n")

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ condition
)

cat(sprintf("  - Design formula: %s\n", deparse(design(dds))))

# ==============================================================================
# RUN DIFFERENTIAL EXPRESSION ANALYSIS
# ==============================================================================

cat("\nRunning DESeq2 analysis...\n")
cat("  - Estimating size factors...\n")
cat("  - Estimating dispersions...\n")
cat("  - Performing Wald test...\n")

# Run the DESeq2 pipeline
dds <- DESeq(dds)

# Get results
cat("\nExtracting results...\n")

res <- results(dds, 
               contrast = c("condition", "Treated", "Control"),
               alpha = ALPHA,
               pAdjustMethod = "BH")

cat(sprintf("  - Comparison: Treated vs Control\n"))
cat(sprintf("  - Alpha threshold: %.2f\n", ALPHA))
cat(sprintf("  - Adjustment method: Benjamini-Hochberg (FDR)\n"))

# ==============================================================================
# SUMMARIZE RESULTS
# ==============================================================================

cat("\n=================================================================\n")
cat("Differential Expression Summary:\n")
cat("=================================================================\n")

summary(res)

# Count significant genes
sig_genes <- res %>%
  as.data.frame() %>%
  filter(!is.na(padj)) %>%
  filter(padj < PADJ_THRESHOLD & abs(log2FoldChange) > LOG2FC_THRESHOLD)

n_sig_up <- sum(sig_genes$log2FoldChange > LOG2FC_THRESHOLD)
n_sig_down <- sum(sig_genes$log2FoldChange < -LOG2FC_THRESHOLD)

cat("\n")
cat(sprintf("Significance thresholds:\n"))
cat(sprintf("  - Adjusted p-value: < %.3f\n", PADJ_THRESHOLD))
cat(sprintf("  - |log2 FC|: > %.2f\n", LOG2FC_THRESHOLD))
cat(sprintf("\nSignificant genes: %d\n", nrow(sig_genes)))
cat(sprintf("  - Upregulated: %d\n", n_sig_up))
cat(sprintf("  - Downregulated: %d\n", n_sig_down))
cat("=================================================================\n\n")

# ==============================================================================
# EXTRACT NORMALIZED COUNTS
# ==============================================================================

cat("Extracting normalized counts...\n")

normalized_counts <- counts(dds, normalized = TRUE)

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

cat("\nSaving results...\n")

# Convert results to data frame
res_df <- res %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  arrange(padj, desc(abs(log2FoldChange)))

# Save full results
write_csv(res_df, DE_RESULTS_FILE)
cat(sprintf("  - Full results: %s\n", DE_RESULTS_FILE))

# Save significant genes only
sig_genes_df <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < PADJ_THRESHOLD & abs(log2FoldChange) > LOG2FC_THRESHOLD)

write_csv(sig_genes_df, DE_SIGNIFICANT_FILE)
cat(sprintf("  - Significant genes: %s\n", DE_SIGNIFICANT_FILE))

# Save normalized counts
normalized_counts_df <- normalized_counts %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")

write_csv(normalized_counts_df, NORMALIZED_COUNTS_FILE)
cat(sprintf("  - Normalized counts: %s\n", NORMALIZED_COUNTS_FILE))

# ==============================================================================
# SAVE DESEQ2 OBJECT
# ==============================================================================

cat("\nSaving DESeq2 object for downstream analysis...\n")

dds_file <- file.path(DE_DIR, "dds_object.rds")
saveRDS(dds, dds_file)
cat(sprintf("  - DESeq2 object: %s\n", dds_file))

# ==============================================================================
# DISPERSION PLOT
# ==============================================================================

cat("\nCreating dispersion plot...\n")

pdf(file.path(QC_DIR, "dispersion_estimates.pdf"), width = 8, height = 6)
plotDispEsts(dds, main = "Dispersion Estimates")
dev.off()
cat("  - Saved dispersion plot\n")

# ==============================================================================
# P-VALUE DISTRIBUTION
# ==============================================================================

cat("\nCreating p-value distribution plot...\n")

pval_data <- data.frame(
  pvalue = res$pvalue
) %>%
  filter(!is.na(pvalue))

p_hist <- ggplot(pval_data, aes(x = pvalue)) +
  geom_histogram(bins = 50, fill = "#3498db", color = "black", alpha = 0.7) +
  theme_classic() +
  labs(title = "P-value Distribution",
       subtitle = "Unadjusted p-values from Wald test",
       x = "P-value",
       y = "Frequency") +
  geom_hline(yintercept = nrow(pval_data) / 50, 
             linetype = "dashed", color = "red", alpha = 0.5)

ggsave(file.path(QC_DIR, "pvalue_distribution.pdf"), p_hist,
       width = FIG_WIDTH, height = FIG_HEIGHT)
cat("  - Saved p-value distribution plot\n")

# ==============================================================================
# LOG2 FOLD CHANGE DISTRIBUTION
# ==============================================================================

cat("\nCreating fold change distribution plot...\n")

fc_data <- data.frame(
  log2FC = res$log2FoldChange,
  significant = !is.na(res$padj) & res$padj < PADJ_THRESHOLD
) %>%
  filter(!is.na(log2FC))

p_fc <- ggplot(fc_data, aes(x = log2FC, fill = significant)) +
  geom_histogram(bins = 100, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("FALSE" = "gray60", "TRUE" = "#e74c3c"),
                    labels = c("FALSE" = "Not significant", "TRUE" = "Significant")) +
  theme_classic() +
  labs(title = "Log2 Fold Change Distribution",
       subtitle = sprintf("Treated vs Control (padj < %.2f)", PADJ_THRESHOLD),
       x = "Log2 Fold Change",
       y = "Frequency",
       fill = "") +
  geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD),
             linetype = "dashed", color = "black", alpha = 0.5) +
  theme(legend.position = "top")

ggsave(file.path(QC_DIR, "fc_distribution.pdf"), p_fc,
       width = FIG_WIDTH, height = FIG_HEIGHT)
cat("  - Saved fold change distribution plot\n")

# ==============================================================================
# ADDITIONAL STATISTICS
# ==============================================================================

cat("\nAdditional statistics:\n")
cat(sprintf("  - Genes with NA padj: %d\n", sum(is.na(res$padj))))
cat(sprintf("  - Mean normalized count (all genes): %.1f\n", mean(rowMeans(normalized_counts))))
cat(sprintf("  - Median |log2FC| (significant genes): %.2f\n", 
            median(abs(sig_genes_df$log2FoldChange))))
cat(sprintf("  - Size factors range: %.3f - %.3f\n", 
            min(sizeFactors(dds)), max(sizeFactors(dds))))

# ==============================================================================
# COMPARE WITH GROUND TRUTH (if available)
# ==============================================================================

truth_file <- file.path(DATA_DIR, "ground_truth.csv")
if (file.exists(truth_file)) {
  cat("\nComparing results with ground truth...\n")
  
  truth <- read_csv(truth_file, show_col_types = FALSE)
  
  # Merge with results
  comparison <- res_df %>%
    left_join(truth, by = "gene_id") %>%
    filter(!is.na(padj))
  
  # Calculate performance metrics
  tp <- sum(comparison$true_DE & comparison$padj < PADJ_THRESHOLD & 
            abs(comparison$log2FoldChange) > LOG2FC_THRESHOLD)
  fp <- sum(!comparison$true_DE & comparison$padj < PADJ_THRESHOLD & 
            abs(comparison$log2FoldChange) > LOG2FC_THRESHOLD)
  fn <- sum(comparison$true_DE & (comparison$padj >= PADJ_THRESHOLD | 
            abs(comparison$log2FoldChange) <= LOG2FC_THRESHOLD))
  tn <- sum(!comparison$true_DE & (comparison$padj >= PADJ_THRESHOLD | 
            abs(comparison$log2FoldChange) <= LOG2FC_THRESHOLD))
  
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  precision <- tp / (tp + fp)
  fdr <- fp / (tp + fp)
  
  cat("\n=================================================================\n")
  cat("Performance Metrics (vs Ground Truth):\n")
  cat("=================================================================\n")
  cat(sprintf("True Positives: %d\n", tp))
  cat(sprintf("False Positives: %d\n", fp))
  cat(sprintf("True Negatives: %d\n", tn))
  cat(sprintf("False Negatives: %d\n", fn))
  cat(sprintf("\nSensitivity (Recall): %.3f\n", sensitivity))
  cat(sprintf("Specificity: %.3f\n", specificity))
  cat(sprintf("Precision: %.3f\n", precision))
  cat(sprintf("FDR: %.3f\n", fdr))
  cat("=================================================================\n\n")
  
  # Save comparison
  comparison_file <- file.path(DE_DIR, "comparison_with_truth.csv")
  write_csv(comparison, comparison_file)
  cat(sprintf("  - Saved comparison: %s\n", comparison_file))
}

cat("\n✓ Differential expression analysis complete!\n\n")

