#!/usr/bin/env Rscript
# ==============================================================================
# Generate Dummy RNA-seq Count Data
# ==============================================================================
# This script generates realistic dummy RNA-seq data for testing the pipeline
# It creates:
#   1. A count matrix with realistic properties
#   2. Sample metadata with condition information
#   3. Simulated differential expression between conditions

cat("\n=================================================================\n")
cat("  Step 1: Generating Dummy RNA-seq Data\n")
cat("=================================================================\n\n")

# Load configuration
source("config.R")
create_directories()

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})

# Set random seed for reproducibility
set.seed(RANDOM_SEED)

# ==============================================================================
# SIMULATE COUNT DATA
# ==============================================================================

cat("Simulating count data...\n")
cat(sprintf("  - Number of genes: %d\n", N_GENES))
cat(sprintf("  - Samples per condition: %d\n", N_SAMPLES_PER_CONDITION))
cat(sprintf("  - Proportion of DE genes: %.1f%%\n", PROP_DE_GENES * 100))

# Total number of samples
n_samples <- N_SAMPLES_PER_CONDITION * 2

# Generate gene names
gene_names <- paste0("GENE", sprintf("%05d", 1:N_GENES))

# Generate sample names and conditions
sample_names <- c(
  paste0("Control_", 1:N_SAMPLES_PER_CONDITION),
  paste0("Treated_", 1:N_SAMPLES_PER_CONDITION)
)

conditions <- factor(
  c(rep("Control", N_SAMPLES_PER_CONDITION), 
    rep("Treated", N_SAMPLES_PER_CONDITION)),
  levels = c("Control", "Treated")
)

# ==============================================================================
# Generate baseline expression levels
# ==============================================================================

# Mean expression follows a log-normal distribution (typical for RNA-seq)
baseline_mu <- rnorm(N_GENES, mean = 6, sd = 2)
baseline_mu <- exp(baseline_mu)

# Add some low-expressed genes and highly-expressed genes
baseline_mu[1:floor(N_GENES * 0.2)] <- baseline_mu[1:floor(N_GENES * 0.2)] * 0.1
baseline_mu[(N_GENES - floor(N_GENES * 0.05)):N_GENES] <- 
  baseline_mu[(N_GENES - floor(N_GENES * 0.05)):N_GENES] * 10

# ==============================================================================
# Simulate differential expression
# ==============================================================================

# Determine which genes are differentially expressed
n_de_genes <- floor(N_GENES * PROP_DE_GENES)
de_genes <- sample(1:N_GENES, n_de_genes)

# Log2 fold changes for DE genes (mixture of up and down)
log2fc <- rep(0, N_GENES)
log2fc[de_genes] <- c(
  rnorm(floor(n_de_genes/2), mean = 2.5, sd = 1),    # Upregulated
  rnorm(ceiling(n_de_genes/2), mean = -2.5, sd = 1)  # Downregulated
)

# Ensure some genes have strong fold changes
strong_de <- sample(de_genes, floor(n_de_genes * 0.3))
log2fc[strong_de] <- log2fc[strong_de] * 1.5

# Convert to linear fold change
fold_change <- 2^log2fc

# ==============================================================================
# Generate count matrix
# ==============================================================================

count_matrix <- matrix(0, nrow = N_GENES, ncol = n_samples)
rownames(count_matrix) <- gene_names
colnames(count_matrix) <- sample_names

# Biological coefficient of variation (typical values for RNA-seq)
bcv <- 0.4  # Higher BCV = more biological variability

for (i in 1:n_samples) {
  # Sample-specific library size
  lib_size <- rnorm(1, mean = MEAN_LIBRARY_SIZE, sd = MEAN_LIBRARY_SIZE * 0.1)
  
  # Apply fold change for treated samples
  if (conditions[i] == "Treated") {
    mu <- baseline_mu * fold_change
  } else {
    mu <- baseline_mu
  }
  
  # Scale by library size
  mu <- mu * lib_size / sum(mu)
  
  # Add biological variability using negative binomial distribution
  # Dispersion parameter (higher = more variance)
  dispersion <- bcv^2 + 1/sqrt(mu)
  size_param <- 1 / dispersion
  
  # Generate counts
  count_matrix[, i] <- rnbinom(N_GENES, mu = mu, size = size_param)
}

# ==============================================================================
# Add some technical noise and batch effects (subtle)
# ==============================================================================

# Small batch effect between first and second half of samples
batch_effect <- c(rep(1, n_samples/2), rep(1.05, n_samples/2))
for (i in 1:n_samples) {
  count_matrix[, i] <- count_matrix[, i] * batch_effect[i]
}

# Round counts to integers
count_matrix <- round(count_matrix)

# Ensure no negative counts
count_matrix[count_matrix < 0] <- 0

# ==============================================================================
# Create sample metadata
# ==============================================================================

cat("\nCreating sample metadata...\n")

sample_metadata <- data.frame(
  sample_id = sample_names,
  condition = conditions,
  replicate = c(1:N_SAMPLES_PER_CONDITION, 1:N_SAMPLES_PER_CONDITION),
  batch = factor(c(
    rep("Batch1", N_SAMPLES_PER_CONDITION/2),
    rep("Batch2", N_SAMPLES_PER_CONDITION/2),
    rep("Batch1", N_SAMPLES_PER_CONDITION/2),
    rep("Batch2", N_SAMPLES_PER_CONDITION/2)
  )),
  library_size = colSums(count_matrix),
  stringsAsFactors = FALSE
)

# ==============================================================================
# Save data to files
# ==============================================================================

cat("\nSaving data files...\n")

# Save count matrix
count_df <- as.data.frame(count_matrix) %>%
  rownames_to_column("gene_id")

write_csv(count_df, RAW_COUNTS_FILE)
cat(sprintf("  - Count matrix: %s\n", RAW_COUNTS_FILE))

# Save sample metadata
write_csv(sample_metadata, SAMPLE_METADATA_FILE)
cat(sprintf("  - Sample metadata: %s\n", SAMPLE_METADATA_FILE))

# ==============================================================================
# Create a truth table for evaluation (which genes are truly DE)
# ==============================================================================

truth_table <- data.frame(
  gene_id = gene_names,
  true_DE = gene_names %in% gene_names[de_genes],
  true_log2FC = log2fc,
  baseline_expression = baseline_mu
)

truth_file <- file.path(DATA_DIR, "ground_truth.csv")
write_csv(truth_table, truth_file)
cat(sprintf("  - Ground truth: %s\n", truth_file))

# ==============================================================================
# Print summary statistics
# ==============================================================================

cat("\n=================================================================\n")
cat("Data Generation Summary:\n")
cat("=================================================================\n")
cat(sprintf("Total genes: %d\n", N_GENES))
cat(sprintf("Total samples: %d\n", n_samples))
cat(sprintf("  - Control: %d\n", N_SAMPLES_PER_CONDITION))
cat(sprintf("  - Treated: %d\n", N_SAMPLES_PER_CONDITION))
cat(sprintf("\nTrue DE genes: %d (%.1f%%)\n", n_de_genes, PROP_DE_GENES * 100))
cat(sprintf("  - Upregulated: %d\n", sum(log2fc > 0)))
cat(sprintf("  - Downregulated: %d\n", sum(log2fc < 0)))
cat("\nLibrary size statistics:\n")
cat(sprintf("  - Mean: %.0f\n", mean(sample_metadata$library_size)))
cat(sprintf("  - Min: %.0f\n", min(sample_metadata$library_size)))
cat(sprintf("  - Max: %.0f\n", max(sample_metadata$library_size)))
cat("\nGene expression statistics:\n")
cat(sprintf("  - Mean counts per gene: %.1f\n", mean(rowMeans(count_matrix))))
cat(sprintf("  - Genes with zero counts across all samples: %d\n", 
            sum(rowSums(count_matrix) == 0)))
cat(sprintf("  - Genes with >1000 mean counts: %d\n", 
            sum(rowMeans(count_matrix) > 1000)))
cat("=================================================================\n\n")

cat("✓ Dummy data generation complete!\n\n")

