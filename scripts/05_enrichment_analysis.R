#!/usr/bin/env Rscript
# ==============================================================================
# Functional Enrichment Analysis
# ==============================================================================
# This script performs functional enrichment analysis:
#   1. Gene Ontology (GO) enrichment
#   2. Pathway analysis
#   3. Visualization of enriched terms

cat("\n=================================================================\n")
cat("  Step 5: Functional Enrichment Analysis\n")
cat("=================================================================\n\n")

# Load configuration
source("config.R")
create_directories()

# Load required libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(enrichplot)
})

# Set random seed
set.seed(RANDOM_SEED)

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading differential expression results...\n")

# Load DE results
res_df <- read_csv(DE_RESULTS_FILE, show_col_types = FALSE)

# Get significant genes
sig_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < PADJ_THRESHOLD & abs(log2FoldChange) > LOG2FC_THRESHOLD)

upregulated_genes <- sig_genes %>%
  filter(log2FoldChange > 0) %>%
  pull(gene_id)

downregulated_genes <- sig_genes %>%
  filter(log2FoldChange < 0) %>%
  pull(gene_id)

all_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  pull(gene_id)

cat(sprintf("  - Total significant genes: %d\n", nrow(sig_genes)))
cat(sprintf("  - Upregulated: %d\n", length(upregulated_genes)))
cat(sprintf("  - Downregulated: %d\n", length(downregulated_genes)))
cat(sprintf("  - Background (all tested genes): %d\n", length(all_genes)))

# ==============================================================================
# NOTE: Using gene symbols as gene IDs
# ==============================================================================
# In a real analysis, you would need to convert gene IDs to appropriate format
# For example, converting gene symbols to Entrez IDs:
# 
# library(AnnotationDbi)
# entrez_ids <- mapIds(org.Hs.eg.db, 
#                      keys = gene_symbols,
#                      column = "ENTREZID",
#                      keytype = "SYMBOL",
#                      multiVals = "first")

cat("\nNote: This analysis uses simulated gene IDs.\n")
cat("In a real analysis, convert your IDs to Entrez Gene IDs for GO analysis.\n\n")

# ==============================================================================
# SIMULATE ENRICHMENT ANALYSIS RESULTS
# ==============================================================================
# Since we're using dummy gene IDs (GENE00001, etc.), we'll simulate
# enrichment results to demonstrate the output format

cat("Simulating GO enrichment analysis results...\n")
cat("(In real analysis, use enrichGO() from clusterProfiler)\n\n")

# Simulate GO enrichment results for upregulated genes
simulate_go_results <- function(gene_list, n_terms = 20) {
  # Sample GO terms
  go_terms <- c(
    "GO:0006355~regulation of transcription, DNA-templated",
    "GO:0045944~positive regulation of transcription from RNA polymerase II promoter",
    "GO:0006351~transcription, DNA-templated",
    "GO:0006367~transcription initiation from RNA polymerase II promoter",
    "GO:0051252~regulation of RNA metabolic process",
    "GO:0016070~RNA metabolic process",
    "GO:0010628~positive regulation of gene expression",
    "GO:0006357~regulation of transcription from RNA polymerase II promoter",
    "GO:0006366~transcription from RNA polymerase II promoter",
    "GO:0045893~positive regulation of transcription, DNA-templated",
    "GO:0008283~cell proliferation",
    "GO:0000122~negative regulation of transcription from RNA polymerase II promoter",
    "GO:0045892~negative regulation of transcription, DNA-templated",
    "GO:0006350~transcription",
    "GO:0006915~apoptotic process",
    "GO:0008219~cell death",
    "GO:0016071~mRNA metabolic process",
    "GO:0006397~mRNA processing",
    "GO:0008380~RNA splicing",
    "GO:0000375~RNA splicing, via transesterification reactions"
  )
  
  n_terms <- min(n_terms, length(go_terms))
  selected_terms <- sample(go_terms, n_terms)
  
  results <- data.frame(
    ID = sapply(strsplit(selected_terms, "~"), `[`, 1),
    Description = sapply(strsplit(selected_terms, "~"), `[`, 2),
    GeneRatio = paste0(sample(5:50, n_terms), "/", length(gene_list)),
    BgRatio = paste0(sample(100:1000, n_terms), "/20000"),
    pvalue = sort(10^(-runif(n_terms, 2, 10))),
    p.adjust = NA,
    qvalue = NA,
    geneID = sapply(1:n_terms, function(i) {
      paste(sample(gene_list, min(10, length(gene_list))), collapse = "/")
    }),
    Count = sample(5:50, n_terms)
  )
  
  results$p.adjust <- p.adjust(results$pvalue, method = "BH")
  results$qvalue <- results$p.adjust
  
  results <- results %>% arrange(pvalue)
  
  return(results)
}

# Create example GO enrichment results
if (length(upregulated_genes) > 0) {
  go_up <- simulate_go_results(upregulated_genes, n_terms = 15)
  cat(sprintf("  - Simulated %d GO terms for upregulated genes\n", nrow(go_up)))
} else {
  go_up <- data.frame()
  cat("  - No upregulated genes for enrichment\n")
}

if (length(downregulated_genes) > 0) {
  go_down <- simulate_go_results(downregulated_genes, n_terms = 15)
  cat(sprintf("  - Simulated %d GO terms for downregulated genes\n", nrow(go_down)))
} else {
  go_down <- data.frame()
  cat("  - No downregulated genes for enrichment\n")
}

# ==============================================================================
# EXAMPLE: Real GO enrichment analysis code
# ==============================================================================
# Uncomment and modify this code for real analysis with proper gene IDs:
#
# if (length(upregulated_genes) > 10) {
#   go_up <- enrichGO(
#     gene = upregulated_genes,
#     universe = all_genes,
#     OrgDb = org.Hs.eg.db,
#     ont = GO_ONTOLOGY,
#     pAdjustMethod = "BH",
#     pvalueCutoff = GO_PVALUE_CUTOFF,
#     qvalueCutoff = GO_QVALUE_CUTOFF,
#     readable = TRUE
#   )
#   
#   go_up <- as.data.frame(go_up)
# }
#
# if (length(downregulated_genes) > 10) {
#   go_down <- enrichGO(
#     gene = downregulated_genes,
#     universe = all_genes,
#     OrgDb = org.Hs.eg.db,
#     ont = GO_ONTOLOGY,
#     pAdjustMethod = "BH",
#     pvalueCutoff = GO_PVALUE_CUTOFF,
#     qvalueCutoff = GO_QVALUE_CUTOFF,
#     readable = TRUE
#   )
#   
#   go_down <- as.data.frame(go_down)
# }

# ==============================================================================
# SAVE ENRICHMENT RESULTS
# ==============================================================================

cat("\nSaving enrichment results...\n")

if (nrow(go_up) > 0) {
  go_up_file <- file.path(ENRICHMENT_DIR, "go_enrichment_upregulated.csv")
  write_csv(go_up, go_up_file)
  cat(sprintf("  - Upregulated GO terms: %s\n", go_up_file))
}

if (nrow(go_down) > 0) {
  go_down_file <- file.path(ENRICHMENT_DIR, "go_enrichment_downregulated.csv")
  write_csv(go_down, go_down_file)
  cat(sprintf("  - Downregulated GO terms: %s\n", go_down_file))
}

# Combined results
go_combined <- bind_rows(
  if (nrow(go_up) > 0) mutate(go_up, Direction = "Upregulated") else NULL,
  if (nrow(go_down) > 0) mutate(go_down, Direction = "Downregulated") else NULL
)

if (nrow(go_combined) > 0) {
  go_combined_file <- file.path(ENRICHMENT_DIR, "go_enrichment_all.csv")
  write_csv(go_combined, go_combined_file)
  cat(sprintf("  - Combined GO terms: %s\n", go_combined_file))
}

# ==============================================================================
# VISUALIZE GO ENRICHMENT RESULTS
# ==============================================================================

cat("\nCreating enrichment visualizations...\n")

# 1. Barplot of top GO terms
if (nrow(go_up) > 0) {
  plot_data_up <- go_up %>%
    head(N_GO_TERMS) %>%
    mutate(Description = factor(Description, levels = rev(Description)))
  
  p_go_up <- ggplot(plot_data_up, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "#e74c3c", high = "#f39c12", 
                       name = "Adjusted\nP-value") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 9),
          plot.title = element_text(face = "bold")) +
    labs(title = "GO Enrichment - Upregulated Genes",
         subtitle = sprintf("Top %d terms (Biological Process)", 
                           min(N_GO_TERMS, nrow(plot_data_up))),
         x = "Gene Count",
         y = "")
  
  ggsave(file.path(ENRICHMENT_DIR, "go_barplot_upregulated.pdf"), p_go_up,
         width = 10, height = max(6, nrow(plot_data_up) * 0.3))
  ggsave(file.path(ENRICHMENT_DIR, "go_barplot_upregulated.png"), p_go_up,
         width = 10, height = max(6, nrow(plot_data_up) * 0.3), dpi = FIG_DPI)
  cat("  - Saved GO barplot (upregulated)\n")
}

if (nrow(go_down) > 0) {
  plot_data_down <- go_down %>%
    head(N_GO_TERMS) %>%
    mutate(Description = factor(Description, levels = rev(Description)))
  
  p_go_down <- ggplot(plot_data_down, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "#3498db", high = "#5dade2",
                       name = "Adjusted\nP-value") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 9),
          plot.title = element_text(face = "bold")) +
    labs(title = "GO Enrichment - Downregulated Genes",
         subtitle = sprintf("Top %d terms (Biological Process)", 
                           min(N_GO_TERMS, nrow(plot_data_down))),
         x = "Gene Count",
         y = "")
  
  ggsave(file.path(ENRICHMENT_DIR, "go_barplot_downregulated.pdf"), p_go_down,
         width = 10, height = max(6, nrow(plot_data_down) * 0.3))
  ggsave(file.path(ENRICHMENT_DIR, "go_barplot_downregulated.png"), p_go_down,
         width = 10, height = max(6, nrow(plot_data_down) * 0.3), dpi = FIG_DPI)
  cat("  - Saved GO barplot (downregulated)\n")
}

# 2. Dot plot
if (nrow(go_combined) > 0) {
  plot_data_combined <- go_combined %>%
    group_by(Direction) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    mutate(
      Description = factor(Description, levels = rev(unique(Description))),
      GeneRatio_numeric = sapply(strsplit(GeneRatio, "/"), 
                                 function(x) as.numeric(x[1]) / as.numeric(x[2]))
    )
  
  p_dot <- ggplot(plot_data_combined, 
                  aes(x = GeneRatio_numeric, y = Description, 
                      color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_gradient(low = "#e74c3c", high = "#3498db",
                        name = "Adjusted\nP-value") +
    scale_size_continuous(name = "Gene\nCount", range = c(3, 10)) +
    facet_wrap(~ Direction, scales = "free_y", ncol = 1) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 9),
          strip.background = element_rect(fill = "gray90", color = "black"),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold")) +
    labs(title = "GO Enrichment Overview",
         subtitle = "Top enriched terms by direction",
         x = "Gene Ratio",
         y = "")
  
  ggsave(file.path(ENRICHMENT_DIR, "go_dotplot.pdf"), p_dot,
         width = 10, height = max(8, nrow(plot_data_combined) * 0.25))
  ggsave(file.path(ENRICHMENT_DIR, "go_dotplot.png"), p_dot,
         width = 10, height = max(8, nrow(plot_data_combined) * 0.25), dpi = FIG_DPI)
  cat("  - Saved GO dotplot\n")
}

# 3. Summary plot
if (nrow(go_combined) > 0) {
  summary_enrich <- go_combined %>%
    group_by(Direction) %>%
    summarise(
      n_terms = n(),
      mean_count = mean(Count),
      .groups = "drop"
    )
  
  p_summary <- ggplot(summary_enrich, aes(x = Direction, y = n_terms, fill = Direction)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    scale_fill_manual(values = c("Upregulated" = "#e74c3c", 
                                  "Downregulated" = "#3498db")) +
    geom_text(aes(label = n_terms), vjust = -0.5, size = 5) +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold")) +
    labs(title = "Enriched GO Terms by Direction",
         subtitle = sprintf("p.adjust < %.2f", GO_PVALUE_CUTOFF),
         x = "",
         y = "Number of Enriched Terms") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  ggsave(file.path(ENRICHMENT_DIR, "enrichment_summary.pdf"), p_summary,
         width = FIG_WIDTH, height = FIG_HEIGHT)
  ggsave(file.path(ENRICHMENT_DIR, "enrichment_summary.png"), p_summary,
         width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_DPI)
  cat("  - Saved enrichment summary plot\n")
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n=================================================================\n")
cat("Enrichment Analysis Summary:\n")
cat("=================================================================\n")
if (nrow(go_up) > 0) {
  cat(sprintf("Upregulated genes:\n"))
  cat(sprintf("  - Enriched GO terms: %d\n", nrow(go_up)))
  cat(sprintf("  - Most significant: %s\n", go_up$Description[1]))
}
if (nrow(go_down) > 0) {
  cat(sprintf("Downregulated genes:\n"))
  cat(sprintf("  - Enriched GO terms: %d\n", nrow(go_down)))
  cat(sprintf("  - Most significant: %s\n", go_down$Description[1]))
}
cat(sprintf("\nResults saved to: %s\n", ENRICHMENT_DIR))
cat("=================================================================\n\n")

cat("NOTE: This analysis uses simulated enrichment results.\n")
cat("For real data:\n")
cat("  1. Convert gene IDs to Entrez Gene IDs\n")
cat("  2. Use enrichGO() or enrichKEGG() from clusterProfiler\n")
cat("  3. Adjust parameters in config.R as needed\n\n")

cat("✓ Enrichment analysis complete!\n\n")

