# Troubleshooting Guide

Common issues and their solutions for the Bulk RNA-seq Analysis Pipeline.

## Installation Issues

### Issue: "Package X is not available for R version Y.Z"

**Solution:**
```R
# Update R to the latest version
# Then update Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")  # or latest version

# Reinstall packages
source("requirements.R")
```

### Issue: "Cannot install DESeq2 / clusterProfiler"

**Solution - macOS:**
```bash
# Install system dependencies
brew install openssl libxml2

# In R:
BiocManager::install("DESeq2")
```

**Solution - Ubuntu/Debian:**
```bash
# Install system dependencies
sudo apt-get update
sudo apt-get install -y libssl-dev libxml2-dev libcurl4-openssl-dev

# In R:
BiocManager::install("DESeq2")
```

**Solution - Windows:**
1. Install Rtools: https://cran.r-project.org/bin/windows/Rtools/
2. Restart R/RStudio
3. Run `source("requirements.R")`

### Issue: "Compilation failed for package X"

**Solution:**
```R
# Try installing pre-compiled binaries
install.packages("X", type = "binary")

# Or for Bioconductor packages:
BiocManager::install("X", type = "binary")
```

## Data Loading Issues

### Issue: "Cannot open file 'data/raw_counts.csv'"

**Possible causes:**
1. Data generation step not run
2. Working directory is incorrect
3. File permissions issue

**Solution:**
```R
# Check working directory
getwd()

# Should be: /path/to/bulkseq
# If not, set it:
setwd("/path/to/bulkseq")

# Generate data
source("scripts/01_generate_dummy_data.R")
```

### Issue: "Error in read_csv: subscript out of bounds"

**Possible causes:**
1. Corrupted CSV file
2. Mismatch between count matrix and metadata

**Solution:**
```R
# Regenerate data
source("scripts/01_generate_dummy_data.R")

# Or check your data format:
counts <- read.csv("data/raw_counts.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("data/sample_metadata.csv")

# Ensure sample names match
colnames(counts) == metadata$sample_id
```

## Analysis Issues

### Issue: "Every gene had p-value set to NA"

**Possible causes:**
1. All genes filtered out
2. No variation in counts
3. Wrong design formula

**Solution:**
```R
# Check filtering
counts <- read.csv("data/filtered_counts.csv")
dim(counts)  # Should have > 0 genes

# Check for zero variance
apply(counts[,-1], 1, var) > 0  # Should be mostly TRUE

# Adjust filtering parameters in config.R
MIN_TOTAL_COUNTS <- 5  # Less stringent
MIN_SAMPLES <- 2
```

### Issue: "Size factor estimation failed"

**Possible causes:**
1. All genes are zeros in some samples
2. Extreme outliers

**Solution:**
```R
# Remove problematic samples
# Or use alternative normalization:
dds <- estimateSizeFactors(dds, type = "poscounts")
```

### Issue: "Dispersion estimation failed"

**Possible causes:**
1. Too few replicates
2. Extreme outliers
3. Wrong design

**Solution:**
```R
# For very small sample sizes:
dds <- estimateDispersions(dds, fitType = "mean")

# Or use a simpler approach:
library(edgeR)
# See edgeR documentation for alternative methods
```

## Visualization Issues

### Issue: "Cannot create PDF/PNG files"

**Possible causes:**
1. Directory doesn't exist
2. Permission issues
3. File is open in another program

**Solution:**
```R
# Create directories
source("config.R")
create_directories()

# Close any open PDF viewers

# Try saving to a different location
ggsave("~/Desktop/test_plot.pdf", plot)
```

### Issue: "Heatmap shows no genes"

**Possible causes:**
1. No significant genes found
2. Too stringent thresholds

**Solution:**
```R
# Adjust thresholds in config.R
PADJ_THRESHOLD <- 0.1  # Less stringent
LOG2FC_THRESHOLD <- 0.5

# Or manually create heatmap with top varying genes
top_var_genes <- head(order(apply(counts, 1, var), decreasing = TRUE), 50)
```

### Issue: "PCA plot shows no separation"

**Possible interpretation:**
1. This may be expected if:
   - True biological differences are small
   - Batch effects dominate
   - Wrong samples are being compared

**Solution:**
```R
# Check metadata
metadata

# Consider batch correction
library(sva)
mat <- assay(vsd)
batch <- metadata$batch
condition <- metadata$condition
mat_corrected <- ComBat(mat, batch = batch, mod = model.matrix(~condition))

# Re-run PCA with corrected data
```

## Memory Issues

### Issue: "Cannot allocate vector of size X GB"

**Solution:**
```R
# Reduce number of genes in config.R
N_GENES <- 5000  # Instead of 10000

# Process data in chunks
# Or use a machine with more RAM

# Close other applications

# Increase R memory limit (Windows only)
memory.limit(size = 16000)  # 16 GB
```

## Enrichment Analysis Issues

### Issue: "No enrichment results found"

**Expected with dummy data** - the simulated gene IDs (GENE00001, etc.) don't correspond to real genes.

**For real data:**
```R
# Ensure proper gene ID format
library(AnnotationDbi)

# Convert gene symbols to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = gene_symbols,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")

# Remove NAs
entrez_ids <- entrez_ids[!is.na(entrez_ids)]

# Run enrichment
ego <- enrichGO(gene = entrez_ids,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)
```

### Issue: "org.Hs.eg.db is for humans, but I have mouse/rat/other data"

**Solution:**
```R
# Install appropriate organism database
# Mouse:
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

# Rat:
BiocManager::install("org.Rn.eg.db")
library(org.Rn.eg.db)

# Others: check https://bioconductor.org/packages/release/BiocViews.html#___OrgDb

# Update config.R
ORGDB <- "org.Mm.eg.db"  # for mouse
```

## Performance Issues

### Issue: "Pipeline is very slow"

**Solution:**
```R
# Enable parallel processing
library(BiocParallel)
register(MulticoreParam(workers = 4))  # Use 4 cores

# Or in config.R
N_CORES <- 4

# Reduce data size for testing
N_GENES <- 1000
N_SAMPLES_PER_CONDITION <- 3
```

## Output Issues

### Issue: "Results look strange / unexpected"

**Diagnostic steps:**
```R
# 1. Check input data
head(counts)
summary(metadata)

# 2. Check filtering
table(genes_pass_filter)

# 3. Check size factors
sizeFactors(dds)
# Should be around 1, not too extreme

# 4. Check dispersion estimates
plotDispEsts(dds)
# Should show decreasing trend

# 5. Check results
summary(res)
# Should have reasonable numbers

# 6. Compare with ground truth (for dummy data)
truth <- read_csv("data/ground_truth.csv")
```

## Platform-Specific Issues

### macOS: "ld: library not found"

```bash
# Install command line tools
xcode-select --install

# Install homebrew dependencies
brew install gcc
```

### Windows: "Rtools is required"

1. Download Rtools from: https://cran.r-project.org/bin/windows/Rtools/
2. Install (default location is fine)
3. Restart R/RStudio
4. Verify: `Sys.which("make")`

### Linux: "Cannot find X11"

```bash
# Ubuntu/Debian
sudo apt-get install libx11-dev

# CentOS/RedHat
sudo yum install libX11-devel
```

## Getting More Help

If your issue isn't covered here:

1. **Check package documentation:**
   ```R
   ?DESeqDataSetFromMatrix
   browseVignettes("DESeq2")
   ```

2. **Search for error messages:**
   - Bioconductor support: https://support.bioconductor.org/
   - Stack Overflow: https://stackoverflow.com/questions/tagged/deseq2

3. **Check package versions:**
   ```R
   sessionInfo()
   packageVersion("DESeq2")
   ```

4. **Enable debugging:**
   ```R
   options(error = recover)  # Interactive debugging
   traceback()  # See where error occurred
   ```

5. **Minimal reproducible example:**
   - Isolate the problem
   - Use dummy data if possible
   - Share code and error message

## Useful Commands

```R
# Check R and package versions
sessionInfo()

# Clear workspace
rm(list = ls())

# Clear plots
dev.off()

# Check available memory
memory.size()  # Windows only
gc()  # Garbage collection

# Reset graphics device
dev.off()
while (!is.null(dev.list())) dev.off()

# Reinstall corrupted package
remove.packages("PackageName")
install.packages("PackageName")

# Update all packages
update.packages(ask = FALSE)
BiocManager::install()
```

---

**Still having issues?** Double-check that:
- Your R version is >= 4.0
- All packages are installed via `source("requirements.R")`
- You're in the correct working directory
- The data files exist and are properly formatted

