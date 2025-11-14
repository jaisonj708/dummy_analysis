#!/usr/bin/env Rscript
# ==============================================================================
# Example Analysis Script
# ==============================================================================
# This script demonstrates how to use the pipeline and explore results
# It can be run after the main pipeline completes

cat("Example Analysis: Exploring RNA-seq Results\n")
cat("============================================\n\n")

# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)

# Load configuration
source("config.R")

# ==============================================================================
# 1. LOAD AND EXPLORE DE RESULTS
# ==============================================================================

cat("1. Loading differential expression results...\n")

# Load full results
de_results <- read_csv(DE_RESULTS_FILE, show_col_types = FALSE)

# Summary statistics
cat("\nBasic statistics:\n")
cat(sprintf("  Total genes tested: %d\n", sum(!is.na(de_results$padj))))
cat(sprintf("  Significant genes (padj < %.2f): %d\n", 
            PADJ_THRESHOLD,
            sum(de_results$padj < PADJ_THRESHOLD, na.rm = TRUE)))

# Top 10 most significant genes
cat("\nTop 10 most significant genes:\n")
top_genes <- de_results %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(10) %>%
  select(gene_id, log2FoldChange, padj)

print(top_genes)

# ==============================================================================
# 2. GENES WITH LARGEST FOLD CHANGES
# ==============================================================================

cat("\n2. Genes with largest fold changes...\n")

# Top upregulated
cat("\nTop 5 upregulated genes:\n")
top_up <- de_results %>%
  filter(!is.na(padj), padj < PADJ_THRESHOLD) %>%
  arrange(desc(log2FoldChange)) %>%
  head(5) %>%
  select(gene_id, log2FoldChange, padj)

print(top_up)

# Top downregulated
cat("\nTop 5 downregulated genes:\n")
top_down <- de_results %>%
  filter(!is.na(padj), padj < PADJ_THRESHOLD) %>%
  arrange(log2FoldChange) %>%
  head(5) %>%
  select(gene_id, log2FoldChange, padj)

print(top_down)

# ==============================================================================
# 3. FILTER GENES BY CUSTOM CRITERIA
# ==============================================================================

cat("\n3. Custom filtering examples...\n")

# Example 1: Highly upregulated genes
highly_upregulated <- de_results %>%
  filter(!is.na(padj),
         padj < 0.01,
         log2FoldChange > 2)

cat(sprintf("\nGenes with padj < 0.01 and log2FC > 2: %d\n", 
            nrow(highly_upregulated)))

# Example 2: Moderately expressed and significant
moderate_de <- de_results %>%
  filter(!is.na(padj),
         padj < PADJ_THRESHOLD,
         baseMean > 100,
         baseMean < 1000,
         abs(log2FoldChange) > 1)

cat(sprintf("Moderately expressed DE genes (100 < baseMean < 1000): %d\n",
            nrow(moderate_de)))

# ==============================================================================
# 4. COMPARE WITH GROUND TRUTH (if available)
# ==============================================================================

truth_file <- file.path(DATA_DIR, "ground_truth.csv")
if (file.exists(truth_file)) {
  cat("\n4. Comparing with ground truth...\n")
  
  truth <- read_csv(truth_file, show_col_types = FALSE)
  
  comparison <- de_results %>%
    left_join(truth, by = "gene_id") %>%
    filter(!is.na(padj))
  
  # Calculate accuracy metrics
  tp <- sum(comparison$true_DE & 
            comparison$padj < PADJ_THRESHOLD & 
            abs(comparison$log2FoldChange) > LOG2FC_THRESHOLD)
  fp <- sum(!comparison$true_DE & 
            comparison$padj < PADJ_THRESHOLD & 
            abs(comparison$log2FoldChange) > LOG2FC_THRESHOLD)
  fn <- sum(comparison$true_DE & 
            (comparison$padj >= PADJ_THRESHOLD | 
             abs(comparison$log2FoldChange) <= LOG2FC_THRESHOLD))
  tn <- sum(!comparison$true_DE & 
            (comparison$padj >= PADJ_THRESHOLD | 
             abs(comparison$log2FoldChange) <= LOG2FC_THRESHOLD))
  
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  precision <- tp / (tp + fp)
  f1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)
  
  cat(sprintf("\nPerformance metrics:\n"))
  cat(sprintf("  Accuracy:    %.3f\n", accuracy))
  cat(sprintf("  Sensitivity: %.3f (recall, TPR)\n", sensitivity))
  cat(sprintf("  Specificity: %.3f (TNR)\n", specificity))
  cat(sprintf("  Precision:   %.3f (PPV)\n", precision))
  cat(sprintf("  F1-score:    %.3f\n", f1_score))
  
  # Plot comparison
  comparison_plot_data <- comparison %>%
    mutate(
      result = case_when(
        true_DE & padj < PADJ_THRESHOLD & abs(log2FoldChange) > LOG2FC_THRESHOLD ~ "True Positive",
        !true_DE & padj < PADJ_THRESHOLD & abs(log2FoldChange) > LOG2FC_THRESHOLD ~ "False Positive",
        true_DE & (padj >= PADJ_THRESHOLD | abs(log2FoldChange) <= LOG2FC_THRESHOLD) ~ "False Negative",
        TRUE ~ "True Negative"
      )
    )
  
  # Correlation of fold changes
  fc_comparison <- comparison_plot_data %>%
    filter(true_DE)
  
  if (nrow(fc_comparison) > 0) {
    fc_cor <- cor(fc_comparison$true_log2FC, fc_comparison$log2FoldChange, 
                  use = "complete.obs")
    
    cat(sprintf("\nCorrelation of true vs estimated log2FC: %.3f\n", fc_cor))
    
    # Create comparison plot
    p_comparison <- ggplot(fc_comparison, 
                          aes(x = true_log2FC, y = log2FoldChange)) +
      geom_point(alpha = 0.5, color = "#3498db") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      theme_classic() +
      labs(title = "True vs Estimated Log2 Fold Changes",
           subtitle = sprintf("Correlation: %.3f", fc_cor),
           x = "True Log2 Fold Change",
           y = "Estimated Log2 Fold Change")
    
    ggsave(file.path(PLOTS_DIR, "fc_comparison.pdf"), p_comparison,
           width = 6, height = 6)
    cat(sprintf("\nSaved fold change comparison plot: %s\n", 
                file.path(PLOTS_DIR, "fc_comparison.pdf")))
  }
}

# ==============================================================================
# 5. EXPORT CUSTOM GENE LISTS
# ==============================================================================

cat("\n5. Creating custom gene lists...\n")

# Create output directory for gene lists
genelist_dir <- file.path(RESULTS_DIR, "gene_lists")
if (!dir.exists(genelist_dir)) {
  dir.create(genelist_dir, recursive = TRUE)
}

# Export different gene categories
categories <- list(
  "all_significant" = de_results %>%
    filter(!is.na(padj), padj < PADJ_THRESHOLD, 
           abs(log2FoldChange) > LOG2FC_THRESHOLD),
  
  "upregulated" = de_results %>%
    filter(!is.na(padj), padj < PADJ_THRESHOLD, 
           log2FoldChange > LOG2FC_THRESHOLD),
  
  "downregulated" = de_results %>%
    filter(!is.na(padj), padj < PADJ_THRESHOLD, 
           log2FoldChange < -LOG2FC_THRESHOLD),
  
  "highly_significant" = de_results %>%
    filter(!is.na(padj), padj < 0.001, 
           abs(log2FoldChange) > LOG2FC_THRESHOLD)
)

for (name in names(categories)) {
  gene_list <- categories[[name]]
  if (nrow(gene_list) > 0) {
    filename <- file.path(genelist_dir, paste0(name, "_genes.csv"))
    write_csv(gene_list, filename)
    cat(sprintf("  Exported %s: %d genes -> %s\n", 
                name, nrow(gene_list), filename))
    
    # Also save just gene IDs (useful for other tools)
    id_filename <- file.path(genelist_dir, paste0(name, "_ids.txt"))
    writeLines(gene_list$gene_id, id_filename)
  }
}

# ==============================================================================
# 6. SUMMARY STATISTICS BY EXPRESSION LEVEL
# ==============================================================================

cat("\n6. DE statistics by expression level...\n")

de_by_expression <- de_results %>%
  filter(!is.na(padj)) %>%
  mutate(
    expression_bin = cut(baseMean, 
                        breaks = c(0, 10, 100, 1000, Inf),
                        labels = c("Low (0-10)", "Medium (10-100)", 
                                  "High (100-1000)", "Very High (>1000)")),
    significant = padj < PADJ_THRESHOLD & abs(log2FoldChange) > LOG2FC_THRESHOLD
  ) %>%
  group_by(expression_bin) %>%
  summarise(
    total = n(),
    significant = sum(significant),
    prop_sig = significant / total,
    .groups = "drop"
  )

cat("\nDE genes by expression level:\n")
print(de_by_expression)

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n============================================\n")
cat("Example analysis complete!\n")
cat("\nAdditional analyses you can perform:\n")
cat("  - Cluster genes by expression patterns\n")
cat("  - Compare multiple contrasts\n")
cat("  - Time-series analysis (if applicable)\n")
cat("  - Integration with other omics data\n")
cat("  - Custom pathway analysis\n")
cat("\nFor more information, consult:\n")
cat("  - DESeq2 vignette: browseVignettes('DESeq2')\n")
cat("  - clusterProfiler book: https://yulab-smu.top/biomedical-knowledge-mining-book/\n")
cat("============================================\n\n")

