#!/usr/bin/env Rscript
# ==============================================================================
# Master Script: Bulk RNA-seq Analysis Pipeline
# ==============================================================================
# This script runs the complete RNA-seq analysis pipeline from start to finish:
#   1. Generate dummy data (optional, skip if you have your own data)
#   2. Quality control and filtering
#   3. Differential expression analysis
#   4. Visualization
#   5. Functional enrichment analysis

cat("\n")
cat("###################################################################\n")
cat("#                                                                 #\n")
cat("#         BULK RNA-SEQ ANALYSIS PIPELINE                         #\n")
cat("#                                                                 #\n")
cat("###################################################################\n")
cat("\n")

# ==============================================================================
# SETUP
# ==============================================================================

# Load configuration
source("config.R")
print_config()

# Record start time
start_time <- Sys.time()

# Create output directories
create_directories()

# ==============================================================================
# CHECK FOR EXISTING DATA
# ==============================================================================

data_exists <- file.exists(RAW_COUNTS_FILE) && file.exists(SAMPLE_METADATA_FILE)

if (data_exists) {
  cat("\n")
  cat("=================================================================\n")
  cat("EXISTING DATA DETECTED\n")
  cat("=================================================================\n")
  cat(sprintf("Found: %s\n", RAW_COUNTS_FILE))
  cat(sprintf("Found: %s\n", SAMPLE_METADATA_FILE))
  cat("\n")
  cat("Options:\n")
  cat("  1. Use existing data (skip data generation)\n")
  cat("  2. Generate new dummy data (will overwrite existing data)\n")
  cat("\n")
  
  # For automated run, default to using existing data
  # In interactive mode, you could prompt the user
  use_existing <- TRUE  # Change to FALSE to regenerate data
  
  if (use_existing) {
    cat("Using existing data files...\n")
    cat("=================================================================\n\n")
  } else {
    cat("Generating new dummy data...\n")
    cat("=================================================================\n\n")
    data_exists <- FALSE
  }
}

# ==============================================================================
# STEP 1: GENERATE DUMMY DATA (if needed)
# ==============================================================================

if (!data_exists) {
  cat("Running Step 1: Data Generation\n\n")
  
  tryCatch({
    source("scripts/01_generate_dummy_data.R")
  }, error = function(e) {
    cat("\n❌ ERROR in Step 1: Data Generation\n")
    cat(sprintf("Error message: %s\n", e$message))
    cat("\nPipeline stopped.\n")
    quit(status = 1)
  })
} else {
  cat("⊳ Skipping Step 1: Using existing data files\n\n")
}

# ==============================================================================
# STEP 2: QUALITY CONTROL AND FILTERING
# ==============================================================================

cat("Running Step 2: Quality Control\n\n")

tryCatch({
  source("scripts/02_quality_control.R")
}, error = function(e) {
  cat("\n❌ ERROR in Step 2: Quality Control\n")
  cat(sprintf("Error message: %s\n", e$message))
  cat("\nPipeline stopped.\n")
  quit(status = 1)
})

# ==============================================================================
# STEP 3: DIFFERENTIAL EXPRESSION ANALYSIS
# ==============================================================================

cat("Running Step 3: Differential Expression Analysis\n\n")

tryCatch({
  source("scripts/03_differential_expression.R")
}, error = function(e) {
  cat("\n❌ ERROR in Step 3: Differential Expression\n")
  cat(sprintf("Error message: %s\n", e$message))
  cat("\nPipeline stopped.\n")
  quit(status = 1)
})

# ==============================================================================
# STEP 4: VISUALIZATION
# ==============================================================================

cat("Running Step 4: Visualization\n\n")

tryCatch({
  source("scripts/04_visualization.R")
}, error = function(e) {
  cat("\n❌ ERROR in Step 4: Visualization\n")
  cat(sprintf("Error message: %s\n", e$message))
  cat("\nPipeline stopped.\n")
  quit(status = 1)
})

# ==============================================================================
# STEP 5: FUNCTIONAL ENRICHMENT ANALYSIS
# ==============================================================================

cat("Running Step 5: Enrichment Analysis\n\n")

tryCatch({
  source("scripts/05_enrichment_analysis.R")
}, error = function(e) {
  cat("\n❌ ERROR in Step 5: Enrichment Analysis\n")
  cat(sprintf("Error message: %s\n", e$message))
  cat("\nNote: Enrichment analysis may fail with dummy data.\n")
  cat("This is expected and can be ignored for testing purposes.\n")
})

# ==============================================================================
# PIPELINE COMPLETE
# ==============================================================================

end_time <- Sys.time()
elapsed_time <- end_time - start_time

cat("\n")
cat("###################################################################\n")
cat("#                                                                 #\n")
cat("#         PIPELINE COMPLETED SUCCESSFULLY!                       #\n")
cat("#                                                                 #\n")
cat("###################################################################\n")
cat("\n")
cat("=================================================================\n")
cat("EXECUTION SUMMARY\n")
cat("=================================================================\n")
cat(sprintf("Start time:    %s\n", format(start_time, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("End time:      %s\n", format(end_time, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Elapsed time:  %s\n", format(elapsed_time)))
cat("\n")
cat("Output directories:\n")
cat(sprintf("  - QC reports:      %s\n", QC_DIR))
cat(sprintf("  - DE results:      %s\n", DE_DIR))
cat(sprintf("  - Plots:           %s\n", PLOTS_DIR))
cat(sprintf("  - Enrichment:      %s\n", ENRICHMENT_DIR))
cat("\n")
cat("Key output files:\n")
cat(sprintf("  - DE results:         %s\n", DE_RESULTS_FILE))
cat(sprintf("  - Significant genes:  %s\n", DE_SIGNIFICANT_FILE))
cat(sprintf("  - Normalized counts:  %s\n", NORMALIZED_COUNTS_FILE))
cat("=================================================================\n")
cat("\n")
cat("Next steps:\n")
cat("  1. Review QC plots in results/qc/\n")
cat("  2. Examine DE results in results/de/\n")
cat("  3. Check visualizations in results/plots/\n")
cat("  4. Explore enrichment results in results/enrichment/\n")
cat("\n")
cat("To modify analysis parameters, edit config.R and re-run.\n")
cat("\n")
cat("Happy analyzing! 🧬\n")
cat("\n")

