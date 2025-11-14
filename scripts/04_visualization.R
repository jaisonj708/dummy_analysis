#!/usr/bin/env Rscript
# ==============================================================================
# Visualization of Differential Expression Results
# ==============================================================================
# This script creates publication-quality visualizations:
#   1. PCA plot of normalized counts
#   2. Heatmap of top differentially expressed genes
#   3. Volcano plot
#   4. MA plot
#   5. Expression patterns of top genes

cat("\n=================================================================\n")
cat("  Step 4: Creating Visualizations\n")
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
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(ggrepel)
  library(tidyr)
})

# Set random seed
set.seed(RANDOM_SEED)

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading data...\n")

# Load DESeq2 object
dds_file <- file.path(DE_DIR, "dds_object.rds")
dds <- readRDS(dds_file)

# Load results
res_df <- read_csv(DE_RESULTS_FILE, show_col_types = FALSE)

# Load metadata
metadata <- read_csv(SAMPLE_METADATA_FILE, show_col_types = FALSE)

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

cat("  - Data loaded successfully\n")

# ==============================================================================
# 1. PCA PLOT
# ==============================================================================

cat("\nCreating PCA plot...\n")

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA
pca_data <- plotPCA(vsd, intgroup = c("condition", "batch"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"), 1)

# Enhanced PCA plot
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = batch)) +
  geom_point(size = 5, alpha = 0.8) +
  scale_color_manual(values = CONDITION_COLORS) +
  scale_shape_manual(values = c(16, 17)) +
  theme_classic() +
  theme(legend.position = "right",
        text = element_text(size = 12),
        plot.title = element_text(face = "bold")) +
  labs(title = "Principal Component Analysis",
       subtitle = "Variance-stabilized counts",
       x = sprintf("PC1: %.1f%% variance", percent_var[1]),
       y = sprintf("PC2: %.1f%% variance", percent_var[2]),
       color = "Condition",
       shape = "Batch") +
  stat_ellipse(aes(group = condition), type = "norm", linetype = 2, 
               show.legend = FALSE, level = 0.95)

ggsave(file.path(PLOTS_DIR, "pca.pdf"), p_pca,
       width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_DPI)
ggsave(file.path(PLOTS_DIR, "pca.png"), p_pca,
       width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_DPI)
cat("  - Saved PCA plot\n")

# ==============================================================================
# 2. VOLCANO PLOT
# ==============================================================================

cat("\nCreating volcano plot...\n")

# Prepare data for volcano plot
volcano_data <- res_df %>%
  mutate(
    significance = case_when(
      is.na(padj) ~ "Not tested",
      padj < PADJ_THRESHOLD & log2FoldChange > LOG2FC_THRESHOLD ~ "Upregulated",
      padj < PADJ_THRESHOLD & log2FoldChange < -LOG2FC_THRESHOLD ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    neg_log10_padj = -log10(padj),
    label = ifelse(padj < PADJ_THRESHOLD & abs(log2FoldChange) > LOG2FC_THRESHOLD * 1.5,
                   gene_id, NA)
  ) %>%
  filter(!is.na(padj), !is.na(log2FoldChange))

# Get top genes to label
top_genes <- volcano_data %>%
  filter(significance %in% c("Upregulated", "Downregulated")) %>%
  arrange(padj) %>%
  head(N_LABEL_GENES)

volcano_data <- volcano_data %>%
  mutate(label = ifelse(gene_id %in% top_genes$gene_id, gene_id, NA))

# Colors
sig_colors <- c(
  "Upregulated" = "#e74c3c",
  "Downregulated" = "#3498db",
  "Not significant" = "gray70",
  "Not tested" = "gray90"
)

# Create volcano plot
p_volcano <- ggplot(volcano_data, aes(x = log2FoldChange, y = neg_log10_padj, 
                                       color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = sig_colors) +
  geom_hline(yintercept = -log10(PADJ_THRESHOLD), linetype = "dashed", 
             color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), 
             linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                  show.legend = FALSE, color = "black") +
  theme_classic() +
  theme(legend.position = "top",
        text = element_text(size = 12),
        plot.title = element_text(face = "bold")) +
  labs(title = "Volcano Plot",
       subtitle = sprintf("Treated vs Control (padj < %.2f, |log2FC| > %.1f)",
                         PADJ_THRESHOLD, LOG2FC_THRESHOLD),
       x = "Log2 Fold Change",
       y = "-Log10 (Adjusted P-value)",
       color = "") +
  xlim(c(-max(abs(volcano_data$log2FoldChange), na.rm = TRUE) - 1,
         max(abs(volcano_data$log2FoldChange), na.rm = TRUE) + 1))

ggsave(file.path(PLOTS_DIR, "volcano_plot.pdf"), p_volcano,
       width = FIG_WIDTH + 1, height = FIG_HEIGHT, dpi = FIG_DPI)
ggsave(file.path(PLOTS_DIR, "volcano_plot.png"), p_volcano,
       width = FIG_WIDTH + 1, height = FIG_HEIGHT, dpi = FIG_DPI)
cat("  - Saved volcano plot\n")

# ==============================================================================
# 3. MA PLOT
# ==============================================================================

cat("\nCreating MA plot...\n")

# Prepare data for MA plot
ma_data <- res_df %>%
  mutate(
    significance = case_when(
      is.na(padj) ~ "Not tested",
      padj < PADJ_THRESHOLD & log2FoldChange > LOG2FC_THRESHOLD ~ "Upregulated",
      padj < PADJ_THRESHOLD & log2FoldChange < -LOG2FC_THRESHOLD ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  ) %>%
  filter(!is.na(baseMean), !is.na(log2FoldChange))

# Create MA plot
p_ma <- ggplot(ma_data, aes(x = log10(baseMean + 1), y = log2FoldChange, 
                             color = significance)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = sig_colors) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", alpha = 0.8) +
  geom_hline(yintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), 
             linetype = "dashed", color = "black", alpha = 0.5) +
  theme_classic() +
  theme(legend.position = "top",
        text = element_text(size = 12),
        plot.title = element_text(face = "bold")) +
  labs(title = "MA Plot",
       subtitle = "Mean expression vs Log2 Fold Change",
       x = "Log10 (Mean Normalized Count + 1)",
       y = "Log2 Fold Change",
       color = "") +
  ylim(c(-max(abs(ma_data$log2FoldChange), na.rm = TRUE) - 1,
         max(abs(ma_data$log2FoldChange), na.rm = TRUE) + 1))

ggsave(file.path(PLOTS_DIR, "ma_plot.pdf"), p_ma,
       width = FIG_WIDTH + 1, height = FIG_HEIGHT, dpi = FIG_DPI)
ggsave(file.path(PLOTS_DIR, "ma_plot.png"), p_ma,
       width = FIG_WIDTH + 1, height = FIG_HEIGHT, dpi = FIG_DPI)
cat("  - Saved MA plot\n")

# ==============================================================================
# 4. HEATMAP OF TOP GENES
# ==============================================================================

cat("\nCreating heatmap of top differentially expressed genes...\n")

# Get top DE genes
top_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < PADJ_THRESHOLD) %>%
  arrange(padj) %>%
  head(N_TOP_GENES) %>%
  pull(gene_id)

if (length(top_genes) > 0) {
  # Get normalized counts for top genes
  top_counts <- normalized_counts[top_genes, ]
  
  # Z-score transformation (scale by row)
  top_counts_scaled <- t(scale(t(top_counts)))
  
  # Prepare annotation
  annotation_col <- data.frame(
    Condition = metadata$condition,
    Batch = metadata$batch,
    row.names = metadata$sample_id
  )
  
  annotation_colors <- list(
    Condition = CONDITION_COLORS,
    Batch = c("Batch1" = "#95a5a6", "Batch2" = "#7f8c8d")
  )
  
  # Create heatmap
  pdf(file.path(PLOTS_DIR, "heatmap_top_genes.pdf"), 
      width = 10, height = max(8, N_TOP_GENES * 0.15))
  
  pheatmap(top_counts_scaled,
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
           scale = "none",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = ifelse(N_TOP_GENES <= 50, TRUE, FALSE),
           show_colnames = TRUE,
           fontsize = 9,
           fontsize_row = 7,
           fontsize_col = 9,
           main = sprintf("Top %d Differentially Expressed Genes", 
                         length(top_genes)),
           border_color = NA)
  
  dev.off()
  
  cat(sprintf("  - Saved heatmap with %d genes\n", length(top_genes)))
} else {
  cat("  - No significant genes found for heatmap\n")
}

# ==============================================================================
# 5. EXPRESSION PATTERNS OF TOP GENES
# ==============================================================================

cat("\nCreating expression pattern plots...\n")

# Get top upregulated and downregulated genes
top_up <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < PADJ_THRESHOLD & log2FoldChange > LOG2FC_THRESHOLD) %>%
  arrange(padj) %>%
  head(6) %>%
  pull(gene_id)

top_down <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < PADJ_THRESHOLD & log2FoldChange < -LOG2FC_THRESHOLD) %>%
  arrange(padj) %>%
  head(6) %>%
  pull(gene_id)

top_genes_pattern <- c(top_up, top_down)

if (length(top_genes_pattern) > 0) {
  # Prepare data for plotting
  pattern_data <- normalized_counts[top_genes_pattern, ] %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "normalized_count") %>%
    left_join(metadata %>% select(sample_id, condition), 
              by = c("sample" = "sample_id"))
  
  # Create boxplot
  p_pattern <- ggplot(pattern_data, aes(x = condition, y = normalized_count, 
                                         fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    scale_fill_manual(values = CONDITION_COLORS) +
    scale_y_continuous(trans = "log10") +
    facet_wrap(~ gene_id, scales = "free_y", ncol = 3) +
    theme_classic() +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "gray90", color = "black"),
          strip.text = element_text(face = "bold")) +
    labs(title = "Expression Patterns of Top DE Genes",
         subtitle = "Normalized counts (log10 scale)",
         x = "Condition",
         y = "Normalized Count")
  
  ggsave(file.path(PLOTS_DIR, "expression_patterns.pdf"), p_pattern,
         width = 10, height = ceiling(length(top_genes_pattern) / 3) * 2.5,
         dpi = FIG_DPI)
  ggsave(file.path(PLOTS_DIR, "expression_patterns.png"), p_pattern,
         width = 10, height = ceiling(length(top_genes_pattern) / 3) * 2.5,
         dpi = FIG_DPI)
  
  cat(sprintf("  - Saved expression patterns for %d genes\n", 
              length(top_genes_pattern)))
} else {
  cat("  - No significant genes found for expression patterns\n")
}

# ==============================================================================
# 6. SUMMARY BARPLOT
# ==============================================================================

cat("\nCreating summary barplot...\n")

# Count genes in each category
summary_data <- data.frame(
  Category = c("Upregulated", "Downregulated", "Not Significant"),
  Count = c(
    sum(res_df$padj < PADJ_THRESHOLD & res_df$log2FoldChange > LOG2FC_THRESHOLD, na.rm = TRUE),
    sum(res_df$padj < PADJ_THRESHOLD & res_df$log2FoldChange < -LOG2FC_THRESHOLD, na.rm = TRUE),
    sum(res_df$padj >= PADJ_THRESHOLD | abs(res_df$log2FoldChange) <= LOG2FC_THRESHOLD, na.rm = TRUE)
  )
)

p_summary <- ggplot(summary_data, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  scale_fill_manual(values = c("Upregulated" = "#e74c3c", 
                                "Downregulated" = "#3498db",
                                "Not Significant" = "gray70")) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 12),
        plot.title = element_text(face = "bold")) +
  labs(title = "Summary of Differential Expression Analysis",
       subtitle = sprintf("Thresholds: padj < %.2f, |log2FC| > %.1f",
                         PADJ_THRESHOLD, LOG2FC_THRESHOLD),
       x = "",
       y = "Number of Genes") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(file.path(PLOTS_DIR, "summary_barplot.pdf"), p_summary,
       width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_DPI)
ggsave(file.path(PLOTS_DIR, "summary_barplot.png"), p_summary,
       width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_DPI)
cat("  - Saved summary barplot\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n=================================================================\n")
cat("Visualization Summary:\n")
cat("=================================================================\n")
cat("Generated plots:\n")
cat("  ✓ PCA plot\n")
cat("  ✓ Volcano plot\n")
cat("  ✓ MA plot\n")
if (length(top_genes) > 0) {
  cat(sprintf("  ✓ Heatmap (%d genes)\n", length(top_genes)))
}
if (length(top_genes_pattern) > 0) {
  cat(sprintf("  ✓ Expression patterns (%d genes)\n", length(top_genes_pattern)))
}
cat("  ✓ Summary barplot\n")
cat(sprintf("\nAll plots saved to: %s\n", PLOTS_DIR))
cat("=================================================================\n\n")

cat("✓ Visualization complete!\n\n")

