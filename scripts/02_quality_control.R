#!/usr/bin/env Rscript
# ==============================================================================
# Quality Control and Data Filtering
# ==============================================================================
# This script performs quality control checks and filters low-quality genes:
#   1. Library size distribution
#   2. Gene detection rates
#   3. Sample correlation analysis
#   4. PCA of raw counts
#   5. Filter low-count genes

cat("\n=================================================================\n")
cat("  Step 2: Quality Control and Filtering\n")
cat("=================================================================\n\n")

# Load configuration
source("config.R")
create_directories()

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(DESeq2)
})

# Set random seed
set.seed(RANDOM_SEED)

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading data...\n")

# Load count matrix
counts <- read_csv(RAW_COUNTS_FILE, show_col_types = FALSE)
count_matrix <- counts %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

# Load sample metadata
metadata <- read_csv(SAMPLE_METADATA_FILE, show_col_types = FALSE)

cat(sprintf("  - Loaded %d genes and %d samples\n", 
            nrow(count_matrix), ncol(count_matrix)))

# ==============================================================================
# 1. LIBRARY SIZE DISTRIBUTION
# ==============================================================================

cat("\nAnalyzing library sizes...\n")

library_sizes <- data.frame(
  sample = colnames(count_matrix),
  total_counts = colSums(count_matrix),
  condition = metadata$condition
)

# Plot library sizes
p1 <- ggplot(library_sizes, aes(x = sample, y = total_counts, fill = condition)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = CONDITION_COLORS) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "top") +
  labs(title = "Library Size Distribution",
       x = "Sample",
       y = "Total Counts",
       fill = "Condition") +
  scale_y_continuous(labels = scales::comma)

ggsave(file.path(QC_DIR, "library_sizes.pdf"), p1, 
       width = FIG_WIDTH, height = FIG_HEIGHT)
cat("  - Saved library size plot\n")

# ==============================================================================
# 2. GENE DETECTION RATES
# ==============================================================================

cat("\nAnalyzing gene detection rates...\n")

gene_detection <- data.frame(
  sample = colnames(count_matrix),
  n_genes_detected = colSums(count_matrix > 0),
  condition = metadata$condition
)

p2 <- ggplot(gene_detection, aes(x = sample, y = n_genes_detected, fill = condition)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = CONDITION_COLORS) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "top") +
  labs(title = "Number of Genes Detected per Sample",
       subtitle = "Genes with count > 0",
       x = "Sample",
       y = "Number of Genes Detected",
       fill = "Condition") +
  scale_y_continuous(labels = scales::comma)

ggsave(file.path(QC_DIR, "gene_detection.pdf"), p2,
       width = FIG_WIDTH, height = FIG_HEIGHT)
cat("  - Saved gene detection plot\n")

# ==============================================================================
# 3. SAMPLE CORRELATION ANALYSIS
# ==============================================================================

cat("\nCalculating sample correlations...\n")

# Log-transform counts (add pseudocount)
log_counts <- log2(count_matrix + 1)

# Calculate correlation matrix
cor_matrix <- cor(log_counts, method = "pearson")

# Create annotation for heatmap
annotation_col <- data.frame(
  Condition = metadata$condition,
  Batch = metadata$batch,
  row.names = metadata$sample_id
)

annotation_colors <- list(
  Condition = CONDITION_COLORS,
  Batch = c("Batch1" = "#95a5a6", "Batch2" = "#7f8c8d")
)

# Plot correlation heatmap
pdf(file.path(QC_DIR, "sample_correlation.pdf"), width = 8, height = 7)
pheatmap(cor_matrix,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Sample-to-Sample Correlation",
         fontsize = 10,
         border_color = NA)
dev.off()
cat("  - Saved correlation heatmap\n")

# ==============================================================================
# 4. PCA OF RAW COUNTS
# ==============================================================================

cat("\nPerforming PCA on raw counts...\n")

# Variance stabilizing transformation for PCA
# Filter genes with very low counts first
genes_to_keep <- rowSums(count_matrix) >= 10
filtered_counts <- count_matrix[genes_to_keep, ]

# Use rlog transformation (more conservative than VST)
# Prepare metadata with factors
metadata_for_deseq <- metadata %>%
  column_to_rownames("sample_id")
metadata_for_deseq$condition <- factor(metadata_for_deseq$condition, levels = c("Control", "Treated"))
metadata_for_deseq$batch <- factor(metadata_for_deseq$batch)

dds_temp <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = metadata_for_deseq,
  design = ~ condition
)

rld <- rlog(dds_temp, blind = TRUE)
rld_mat <- assay(rld)

# Perform PCA
pca_result <- prcomp(t(rld_mat), scale. = FALSE)
pca_data <- data.frame(
  sample = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  condition = metadata$condition,
  batch = metadata$batch
)

# Calculate variance explained
var_explained <- summary(pca_result)$importance[2, ] * 100

# Plot PCA
p3 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = batch)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = CONDITION_COLORS) +
  theme_classic() +
  theme(legend.position = "right") +
  labs(title = "PCA of Raw Counts",
       subtitle = "Based on rlog-transformed counts",
       x = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
       y = sprintf("PC2 (%.1f%% variance)", var_explained[2]),
       color = "Condition",
       shape = "Batch")

ggsave(file.path(QC_DIR, "pca_raw.pdf"), p3,
       width = FIG_WIDTH, height = FIG_HEIGHT)
cat("  - Saved PCA plot\n")

# ==============================================================================
# 5. FILTER LOW-COUNT GENES
# ==============================================================================

cat("\nFiltering low-count genes...\n")

# Filter criteria:
# 1. Total counts across all samples >= MIN_TOTAL_COUNTS
# 2. Detected in at least MIN_SAMPLES samples

total_counts <- rowSums(count_matrix)
n_samples_detected <- rowSums(count_matrix >= MIN_COUNT_THRESHOLD)

genes_pass_filter <- (total_counts >= MIN_TOTAL_COUNTS) & 
                     (n_samples_detected >= MIN_SAMPLES)

filtered_count_matrix <- count_matrix[genes_pass_filter, ]

cat(sprintf("  - Genes before filtering: %d\n", nrow(count_matrix)))
cat(sprintf("  - Genes after filtering: %d\n", nrow(filtered_count_matrix)))
cat(sprintf("  - Genes removed: %d (%.1f%%)\n", 
            sum(!genes_pass_filter),
            100 * sum(!genes_pass_filter) / nrow(count_matrix)))

# ==============================================================================
# SAVE FILTERED DATA
# ==============================================================================

cat("\nSaving filtered count matrix...\n")

filtered_counts_df <- as.data.frame(filtered_count_matrix) %>%
  rownames_to_column("gene_id")

filtered_counts_file <- file.path(DATA_DIR, "filtered_counts.csv")
write_csv(filtered_counts_df, filtered_counts_file)
cat(sprintf("  - Filtered counts: %s\n", filtered_counts_file))

# Save filtering statistics
filter_stats <- data.frame(
  gene_id = rownames(count_matrix),
  total_counts = total_counts,
  n_samples_detected = n_samples_detected,
  passed_filter = genes_pass_filter
)

filter_stats_file <- file.path(QC_DIR, "filtering_stats.csv")
write_csv(filter_stats, filter_stats_file)
cat(sprintf("  - Filter statistics: %s\n", filter_stats_file))

# ==============================================================================
# PLOT FILTERING RESULTS
# ==============================================================================

cat("\nCreating filtering summary plots...\n")

# Mean-variance plot
mean_counts <- rowMeans(count_matrix)
var_counts <- apply(count_matrix, 1, var)

filter_plot_data <- data.frame(
  mean = mean_counts,
  variance = var_counts,
  passed = genes_pass_filter
)

p4 <- ggplot(filter_plot_data, aes(x = mean, y = variance, color = passed)) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  scale_color_manual(values = c("TRUE" = "#27ae60", "FALSE" = "#e74c3c"),
                     labels = c("TRUE" = "Kept", "FALSE" = "Removed")) +
  theme_classic() +
  labs(title = "Mean-Variance Relationship",
       subtitle = "Genes colored by filtering status",
       x = "Mean Count (log10)",
       y = "Variance (log10)",
       color = "Filter Status")

ggsave(file.path(QC_DIR, "mean_variance.pdf"), p4,
       width = FIG_WIDTH, height = FIG_HEIGHT)
cat("  - Saved mean-variance plot\n")

# ==============================================================================
# SUMMARY REPORT
# ==============================================================================

cat("\n=================================================================\n")
cat("Quality Control Summary:\n")
cat("=================================================================\n")
cat(sprintf("Initial genes: %d\n", nrow(count_matrix)))
cat(sprintf("Filtered genes: %d\n", nrow(filtered_count_matrix)))
cat(sprintf("Genes removed: %d\n", sum(!genes_pass_filter)))
cat("\nLibrary size range:\n")
cat(sprintf("  - Min: %s\n", format(min(library_sizes$total_counts), big.mark = ",")))
cat(sprintf("  - Max: %s\n", format(max(library_sizes$total_counts), big.mark = ",")))
cat(sprintf("  - Mean: %s\n", format(round(mean(library_sizes$total_counts)), big.mark = ",")))
cat("\nGenes detected per sample:\n")
cat(sprintf("  - Min: %d\n", min(gene_detection$n_genes_detected)))
cat(sprintf("  - Max: %d\n", max(gene_detection$n_genes_detected)))
cat(sprintf("  - Mean: %d\n", round(mean(gene_detection$n_genes_detected))))
cat("\nSample correlation:\n")
cat(sprintf("  - Min: %.3f\n", min(cor_matrix[lower.tri(cor_matrix)])))
cat(sprintf("  - Max: %.3f\n", max(cor_matrix[lower.tri(cor_matrix)])))
cat(sprintf("  - Mean: %.3f\n", mean(cor_matrix[lower.tri(cor_matrix)])))
cat("=================================================================\n\n")

cat("✓ Quality control complete!\n\n")

