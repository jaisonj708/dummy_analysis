# Install Required R Packages for Bulk RNA-seq Analysis Pipeline
# Run this script once to install all dependencies

cat("Installing required packages for Bulk RNA-seq Analysis Pipeline...\n\n")

# ==============================================================================
# Install BiocManager if not already installed
# ==============================================================================
if (!require("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager")
}

# ==============================================================================
# Bioconductor Packages
# ==============================================================================
bioc_packages <- c(
  "DESeq2",              # Differential expression analysis
  "edgeR",               # Alternative DE analysis method
  "clusterProfiler",     # Functional enrichment analysis
  "org.Hs.eg.db",        # Human genome annotation
  "EnhancedVolcano",     # Enhanced volcano plots
  "ComplexHeatmap",      # Advanced heatmap plotting
  "limma",               # Linear models for microarray/RNA-seq
  "AnnotationDbi"        # Annotation database interface
)

cat("Installing Bioconductor packages...\n")
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(sprintf("  %s already installed\n", pkg))
  }
}

# ==============================================================================
# CRAN Packages
# ==============================================================================
cran_packages <- c(
  "ggplot2",             # Publication-quality graphics
  "pheatmap",            # Pretty heatmaps
  "RColorBrewer",        # Color palettes
  "dplyr",               # Data manipulation
  "tidyr",               # Data tidying
  "readr",               # Fast CSV reading/writing
  "viridis",             # Colorblind-friendly palettes
  "ggrepel",             # Better text labels in ggplot2
  "tibble",              # Modern data frames
  "cowplot",             # Combining plots
  "scales",              # Scale functions for visualization
  "ggpubr"               # Publication-ready plots
)

cat("\nInstalling CRAN packages...\n")
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  } else {
    cat(sprintf("  %s already installed\n", pkg))
  }
}

# ==============================================================================
# Verify Installation
# ==============================================================================
cat("\n=================================================================\n")
cat("Verifying package installation...\n")
cat("=================================================================\n\n")

all_packages <- c(bioc_packages, cran_packages)
failed_packages <- c()

for (pkg in all_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("✓ %s\n", pkg))
  } else {
    cat(sprintf("✗ %s FAILED\n", pkg))
    failed_packages <- c(failed_packages, pkg)
  }
}

cat("\n=================================================================\n")
if (length(failed_packages) == 0) {
  cat("SUCCESS! All packages installed successfully.\n")
  cat("You can now run the RNA-seq analysis pipeline.\n")
} else {
  cat("WARNING! The following packages failed to install:\n")
  for (pkg in failed_packages) {
    cat(sprintf("  - %s\n", pkg))
  }
  cat("\nPlease install these packages manually.\n")
}
cat("=================================================================\n\n")

# ==============================================================================
# Print R and Package Version Information
# ==============================================================================
cat("\nR and Package Version Information:\n")
cat("=================================================================\n")
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("Bioconductor version: %s\n", BiocManager::version()))
cat("\nKey package versions:\n")
key_packages <- c("DESeq2", "ggplot2", "clusterProfiler")
for (pkg in key_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
  }
}
cat("=================================================================\n")

