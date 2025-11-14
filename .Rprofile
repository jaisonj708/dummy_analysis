# .Rprofile for Bulk RNA-seq Analysis Pipeline
# This file is automatically sourced when R starts in this directory

# Set CRAN mirror
local({
  r <- getOption("repos")
  r["CRAN"] <- "https://cloud.r-project.org/"
  options(repos = r)
})

# Set default Bioconductor version
options(BioC_mirror = "https://bioconductor.org")

# Increase timeout for large package downloads
options(timeout = 300)

# Better printing of tibbles
options(
  tibble.print_max = 20,
  tibble.print_min = 10
)

# Disable scientific notation for readability
options(scipen = 999)

# Set number of cores for parallel processing
# Detect number of cores and use all but one
n_cores <- parallel::detectCores() - 1
if (n_cores < 1) n_cores <- 1
options(mc.cores = n_cores)

# Startup message
cat("\n")
cat("========================================\n")
cat("  Bulk RNA-seq Analysis Pipeline\n")
cat("========================================\n")
cat(sprintf("Working directory: %s\n", getwd()))
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("Cores available: %d\n", n_cores))
cat("\n")
cat("Quick start:\n")
cat("  source('scripts/run_pipeline.R')  # Run complete pipeline\n")
cat("  source('requirements.R')          # Install packages\n")
cat("\n")
cat("For help, see README.md or QUICKSTART.md\n")
cat("========================================\n")
cat("\n")

