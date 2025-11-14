# Configuration file for Bulk RNA-seq Analysis Pipeline
# Modify these parameters to customize your analysis

# ==============================================================================
# GENERAL SETTINGS
# ==============================================================================

# Random seed for reproducibility
RANDOM_SEED <- 42

# Number of CPU cores to use (set to NULL for automatic detection)
N_CORES <- 4

# ==============================================================================
# DATA GENERATION PARAMETERS (for dummy data)
# ==============================================================================

# Number of genes to simulate
N_GENES <- 10000

# Number of samples per condition
N_SAMPLES_PER_CONDITION <- 6

# Proportion of differentially expressed genes
PROP_DE_GENES <- 0.15

# Mean library size (total counts per sample)
MEAN_LIBRARY_SIZE <- 1e6

# ==============================================================================
# FILTERING PARAMETERS
# ==============================================================================

# Minimum total counts across all samples
MIN_TOTAL_COUNTS <- 10

# Minimum number of samples a gene must be detected in
MIN_SAMPLES <- 3

# Minimum counts to consider a gene "detected" in a sample
MIN_COUNT_THRESHOLD <- 1

# ==============================================================================
# DIFFERENTIAL EXPRESSION PARAMETERS
# ==============================================================================

# Adjusted p-value threshold
PADJ_THRESHOLD <- 0.05

# Log2 fold change threshold
LOG2FC_THRESHOLD <- 1

# Independent filtering alpha (DESeq2)
ALPHA <- 0.05

# ==============================================================================
# VISUALIZATION PARAMETERS
# ==============================================================================

# Number of top genes to show in heatmap
N_TOP_GENES <- 50

# Number of genes to label in volcano plot
N_LABEL_GENES <- 10

# PCA components to plot
PCA_COMPONENTS <- c(1, 2)

# Figure width and height (inches)
FIG_WIDTH <- 8
FIG_HEIGHT <- 6

# DPI for saved figures
FIG_DPI <- 300

# Color palette for conditions
CONDITION_COLORS <- c("Control" = "#3498db", "Treated" = "#e74c3c")

# ==============================================================================
# ENRICHMENT ANALYSIS PARAMETERS
# ==============================================================================

# Organism database (for human)
ORGDB <- "org.Hs.eg.db"

# GO term analysis parameters
GO_PVALUE_CUTOFF <- 0.05
GO_QVALUE_CUTOFF <- 0.2
GO_ONTOLOGY <- "BP"  # Biological Process (can be "BP", "MF", "CC")

# Number of GO terms to show in plots
N_GO_TERMS <- 20

# ==============================================================================
# DIRECTORY PATHS
# ==============================================================================

# Data directories
DATA_DIR <- "data"
RESULTS_DIR <- "results"

# Output subdirectories
QC_DIR <- file.path(RESULTS_DIR, "qc")
DE_DIR <- file.path(RESULTS_DIR, "de")
PLOTS_DIR <- file.path(RESULTS_DIR, "plots")
ENRICHMENT_DIR <- file.path(RESULTS_DIR, "enrichment")

# Create directories if they don't exist
create_directories <- function() {
  dirs <- c(DATA_DIR, RESULTS_DIR, QC_DIR, DE_DIR, PLOTS_DIR, ENRICHMENT_DIR)
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat(sprintf("Created directory: %s\n", dir))
    }
  }
}

# ==============================================================================
# FILE PATHS
# ==============================================================================

# Input files
RAW_COUNTS_FILE <- file.path(DATA_DIR, "raw_counts.csv")
SAMPLE_METADATA_FILE <- file.path(DATA_DIR, "sample_metadata.csv")

# Output files
NORMALIZED_COUNTS_FILE <- file.path(DE_DIR, "normalized_counts.csv")
DE_RESULTS_FILE <- file.path(DE_DIR, "de_results.csv")
DE_SIGNIFICANT_FILE <- file.path(DE_DIR, "de_significant.csv")

# ==============================================================================
# PRINT CONFIGURATION
# ==============================================================================

print_config <- function() {
  cat("\n")
  cat("=================================================================\n")
  cat("  Bulk RNA-seq Analysis Pipeline Configuration\n")
  cat("=================================================================\n")
  cat(sprintf("Random seed: %d\n", RANDOM_SEED))
  cat(sprintf("Number of cores: %s\n", ifelse(is.null(N_CORES), "auto", N_CORES)))
  cat("\n")
  cat("Filtering thresholds:\n")
  cat(sprintf("  - Min total counts: %d\n", MIN_TOTAL_COUNTS))
  cat(sprintf("  - Min samples: %d\n", MIN_SAMPLES))
  cat("\n")
  cat("DE significance thresholds:\n")
  cat(sprintf("  - Adjusted p-value: %.3f\n", PADJ_THRESHOLD))
  cat(sprintf("  - |log2 FC|: %.2f\n", LOG2FC_THRESHOLD))
  cat("\n")
  cat("Output directories:\n")
  cat(sprintf("  - QC: %s\n", QC_DIR))
  cat(sprintf("  - DE: %s\n", DE_DIR))
  cat(sprintf("  - Plots: %s\n", PLOTS_DIR))
  cat(sprintf("  - Enrichment: %s\n", ENRICHMENT_DIR))
  cat("=================================================================\n")
  cat("\n")
}

