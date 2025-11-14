# Complete Usage Guide

This guide provides detailed instructions for using the Bulk RNA-seq Analysis Pipeline.

## Table of Contents

1. [Installation](#installation)
2. [Running the Pipeline](#running-the-pipeline)
3. [Understanding the Output](#understanding-the-output)
4. [Customizing the Analysis](#customizing-the-analysis)
5. [Using Your Own Data](#using-your-own-data)
6. [Advanced Usage](#advanced-usage)

---

## Installation

### Step 1: Install R

Download and install R (version >= 4.0) from: https://www.r-project.org/

**Recommended:** Also install RStudio from: https://posit.co/download/rstudio-desktop/

### Step 2: Install Required Packages

Open R in the project directory and run:

```R
source("requirements.R")
```

This will install all necessary packages. **Expected time: 10-20 minutes**

Package categories:
- **Bioconductor packages:** DESeq2, clusterProfiler, org.Hs.eg.db, etc.
- **CRAN packages:** ggplot2, dplyr, pheatmap, etc.

---

## Running the Pipeline

### Method 1: Complete Pipeline (Recommended)

Run everything in one go:

```R
source("scripts/run_pipeline.R")
```

This executes all 5 steps sequentially:
1. Generate dummy data
2. Quality control
3. Differential expression
4. Visualization
5. Enrichment analysis

**Expected time:** 3-5 minutes

### Method 2: Step-by-Step Execution

Run individual scripts for more control:

```R
# Step 1: Generate dummy data
source("scripts/01_generate_dummy_data.R")

# Step 2: Quality control
source("scripts/02_quality_control.R")

# Step 3: Differential expression analysis
source("scripts/03_differential_expression.R")

# Step 4: Create visualizations
source("scripts/04_visualization.R")

# Step 5: Enrichment analysis
source("scripts/05_enrichment_analysis.R")
```

### Method 3: Command Line Execution

From terminal/command prompt:

```bash
cd /path/to/bulkseq
Rscript scripts/run_pipeline.R
```

Or individual scripts:

```bash
Rscript scripts/01_generate_dummy_data.R
Rscript scripts/02_quality_control.R
# ... etc
```

---

## Understanding the Output

### Directory Structure After Running

```
results/
├── qc/                          # Quality control outputs
│   ├── library_sizes.pdf
│   ├── gene_detection.pdf
│   ├── sample_correlation.pdf
│   ├── pca_raw.pdf
│   ├── dispersion_estimates.pdf
│   ├── pvalue_distribution.pdf
│   ├── fc_distribution.pdf
│   ├── mean_variance.pdf
│   └── filtering_stats.csv
│
├── de/                          # Differential expression results
│   ├── de_results.csv          # Full results table
│   ├── de_significant.csv      # Significant genes only
│   ├── normalized_counts.csv   # Normalized expression matrix
│   ├── dds_object.rds          # DESeq2 object
│   └── comparison_with_truth.csv
│
├── plots/                       # Visualizations
│   ├── pca.pdf
│   ├── volcano_plot.pdf
│   ├── ma_plot.pdf
│   ├── heatmap_top_genes.pdf
│   ├── expression_patterns.pdf
│   └── summary_barplot.pdf
│
└── enrichment/                  # Functional enrichment
    ├── go_enrichment_all.csv
    ├── go_enrichment_upregulated.csv
    ├── go_enrichment_downregulated.csv
    ├── go_barplot_upregulated.pdf
    ├── go_barplot_downregulated.pdf
    ├── go_dotplot.pdf
    └── enrichment_summary.pdf
```

### Key Output Files

#### 1. `de_results.csv` - Complete DE Results

Columns:
- `gene_id`: Gene identifier
- `baseMean`: Mean normalized counts across all samples
- `log2FoldChange`: Log2 fold change (Treated vs Control)
- `lfcSE`: Standard error of log2FC
- `stat`: Wald test statistic
- `pvalue`: Nominal p-value
- `padj`: Benjamini-Hochberg adjusted p-value (FDR)

**How to read:**
- Positive log2FC = upregulated in Treated
- Negative log2FC = downregulated in Treated
- padj < 0.05 = statistically significant

#### 2. `de_significant.csv` - Filtered Significant Genes

Contains only genes meeting significance thresholds:
- padj < 0.05 (configurable)
- |log2FC| > 1 (configurable)

**Use this file for:**
- Downstream pathway analysis
- Experimental validation
- Follow-up studies

#### 3. `normalized_counts.csv` - Normalized Expression Matrix

DESeq2-normalized counts for all genes and samples.

**Uses:**
- Visualization
- Clustering analysis
- Expression comparisons
- Input for other tools

---

## Customizing the Analysis

All parameters are in `config.R`. Edit and re-run the pipeline.

### Common Customizations

#### 1. Change Significance Thresholds

```R
# More stringent
PADJ_THRESHOLD <- 0.01      # Default: 0.05
LOG2FC_THRESHOLD <- 2       # Default: 1

# Less stringent
PADJ_THRESHOLD <- 0.1
LOG2FC_THRESHOLD <- 0.5
```

#### 2. Modify Data Generation

```R
# More genes
N_GENES <- 20000            # Default: 10000

# More samples
N_SAMPLES_PER_CONDITION <- 10  # Default: 6

# More DE genes
PROP_DE_GENES <- 0.30       # Default: 0.15 (15%)

# Larger libraries
MEAN_LIBRARY_SIZE <- 5e6    # Default: 1e6
```

#### 3. Change Filtering Parameters

```R
# More stringent filtering
MIN_TOTAL_COUNTS <- 20      # Default: 10
MIN_SAMPLES <- 5            # Default: 3

# Less stringent filtering
MIN_TOTAL_COUNTS <- 5
MIN_SAMPLES <- 2
```

#### 4. Adjust Visualization

```R
# More genes in heatmap
N_TOP_GENES <- 100          # Default: 50

# More labels in volcano plot
N_LABEL_GENES <- 20         # Default: 10

# Larger figures
FIG_WIDTH <- 10             # Default: 8
FIG_HEIGHT <- 8             # Default: 6
FIG_DPI <- 600              # Default: 300
```

---

## Using Your Own Data

### Data Format Requirements

#### 1. Count Matrix (`data/raw_counts.csv`)

Format:
```
gene_id,Sample1,Sample2,Sample3,...
GENE1,100,150,120,...
GENE2,50,60,55,...
GENE3,0,5,10,...
```

**Requirements:**
- First column: gene IDs (unique)
- Remaining columns: integer counts for each sample
- No missing values (use 0 for undetected genes)
- Raw counts (not normalized, not TPM/FPKM)

#### 2. Sample Metadata (`data/sample_metadata.csv`)

Format:
```
sample_id,condition,replicate,batch
Sample1,Control,1,Batch1
Sample2,Control,2,Batch1
Sample3,Treated,1,Batch1
```

**Requirements:**
- `sample_id`: Must match column names in count matrix
- `condition`: Treatment groups (e.g., Control, Treated)
- `replicate`: Sample number within each condition
- `batch`: Optional batch information
- Additional columns allowed

### Converting from Other Formats

#### From FASTQ files

Use tools like:
- **Alignment:** STAR, HISAT2
- **Quantification:** featureCounts, HTSeq, RSEM

Example workflow:
```bash
# Align with STAR
STAR --genomeDir genome_index --readFilesIn sample1_R1.fastq.gz

# Count with featureCounts
featureCounts -a genes.gtf -o counts.txt aligned.bam
```

Then convert counts to CSV format.

#### From Salmon/Kallisto output

```R
library(tximport)

# Import transcript-level counts
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# Extract gene-level counts
counts <- txi$counts

# Save
write.csv(counts, "data/raw_counts.csv")
```

#### From GEO/SRA datasets

```R
library(GEOquery)

# Download from GEO
gse <- getGEO("GSE12345")
counts <- exprs(gse[[1]])

# Process and save
# (Format depends on dataset)
```

### Modifying the Design

Edit `scripts/03_differential_expression.R`:

```R
# Change design formula (line ~60)
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ batch + condition  # Add batch effect
)

# Change contrast (line ~90)
res <- results(dds, 
               contrast = c("condition", "Treatment1", "Treatment2"))
```

---

## Advanced Usage

### Multiple Comparisons

To compare multiple conditions:

```R
# In 03_differential_expression.R, after running DESeq():

# Comparison 1: Treated vs Control
res1 <- results(dds, contrast = c("condition", "Treated", "Control"))

# Comparison 2: Another comparison
res2 <- results(dds, contrast = c("condition", "TreatmentA", "TreatmentB"))

# Save both
write_csv(as.data.frame(res1), "results/de/treated_vs_control.csv")
write_csv(as.data.frame(res2), "results/de/treatmentA_vs_B.csv")
```

### Batch Effect Correction

```R
# Option 1: Include batch in design
design(dds) <- ~ batch + condition

# Option 2: Use limma removeBatchEffect (for visualization only)
library(limma)
mat <- assay(vsd)
mat_corrected <- removeBatchEffect(mat, batch = metadata$batch)
```

### Time Series Analysis

For time-course experiments:

```R
# Likelihood ratio test
dds <- DESeq(dds, test = "LRT", reduced = ~ 1)

# Or use specific contrasts at each timepoint
```

### Weighted Analysis

For samples with different qualities:

```R
# Calculate weights (e.g., based on RNA quality)
weights <- metadata$RIN / 10  # Normalize to 0-1

# Not directly supported in DESeq2
# Consider using limma-voom with weights instead
```

### Gene Set Testing

Custom gene sets:

```R
library(fgsea)

# Prepare ranked gene list
ranks <- setNames(res_df$stat, res_df$gene_id)

# Load custom gene sets
pathways <- gmtPathways("custom_pathways.gmt")

# Run GSEA
fgsea_results <- fgsea(pathways, ranks, nperm = 10000)
```

### Interactive Exploration

```R
# Create interactive volcano plot
library(plotly)

p <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj),
                              text = gene_id)) +
  geom_point(aes(color = significance))

ggplotly(p)

# Or use specialized tools
library(Glimma)
glMDPlot(dds)
```

### Exporting for Other Tools

```R
# For Cytoscape
write.table(sig_genes_df[, c("gene_id", "log2FoldChange")],
            "results/cytoscape_input.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# For GSEA software
library(DESeq2)
# Export ranked list
rnk <- data.frame(gene = res_df$gene_id, rank = res_df$stat)
write.table(rnk, "results/gsea_ranked.rnk", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# For IPA, Metascape, etc.
# Usually just need significant gene list
writeLines(sig_genes_df$gene_id, "results/gene_list.txt")
```

### Automation and Scripting

Create a custom wrapper:

```R
# run_custom_analysis.R
run_bulkseq_pipeline <- function(
  count_file,
  metadata_file,
  output_dir,
  padj_threshold = 0.05,
  lfc_threshold = 1
) {
  # Load data
  counts <- read_csv(count_file)
  metadata <- read_csv(metadata_file)
  
  # Run analysis
  # ... (adapt from main scripts)
  
  # Return results
  return(list(
    de_results = res_df,
    plots = plots,
    enrichment = enrichment
  ))
}

# Use it
results <- run_bulkseq_pipeline(
  "data/my_counts.csv",
  "data/my_metadata.csv",
  "results/my_analysis"
)
```

---

## Tips and Best Practices

### 1. Always Check QC First

Before trusting DE results:
- ✓ PCA shows expected clustering
- ✓ Sample correlation is high within groups
- ✓ Library sizes are similar
- ✓ No obvious outliers

### 2. Biological Replicates Are Essential

- Minimum: 3 replicates per condition
- Recommended: 5-6 replicates
- More replicates = more statistical power

### 3. Filter Appropriately

- Too lenient = noise and false positives
- Too stringent = miss low-expressed genes
- Default settings usually work well

### 4. Multiple Testing Correction

- Always use adjusted p-values (padj)
- Never use raw p-values for significance
- FDR < 0.05 is standard

### 5. Effect Size Matters

- Consider both statistical significance AND biological relevance
- A gene with padj = 0.001 but log2FC = 0.1 may not be biologically important
- Set meaningful fold change thresholds

### 6. Validate Results

- Check a few genes by qRT-PCR
- Look for known markers
- Compare with literature
- Examine expression patterns visually

### 7. Document Your Analysis

- Save your modified config.R
- Note any manual changes
- Record R session info: `sessionInfo()`
- Keep a lab notebook entry

---

## Common Workflows

### Workflow 1: Standard Two-Group Comparison

```R
# 1. Setup
source("requirements.R")
source("config.R")

# 2. Prepare data (your own or dummy)
source("scripts/01_generate_dummy_data.R")

# 3. Complete analysis
source("scripts/run_pipeline.R")

# 4. Explore results
source("example_analysis.R")
```

### Workflow 2: Custom Analysis with Your Data

```R
# 1. Place your files
# - data/raw_counts.csv
# - data/sample_metadata.csv

# 2. Adjust parameters
# - Edit config.R

# 3. Run analysis (skip data generation)
source("scripts/02_quality_control.R")
source("scripts/03_differential_expression.R")
source("scripts/04_visualization.R")
source("scripts/05_enrichment_analysis.R")

# 4. Custom exploration
source("example_analysis.R")
```

### Workflow 3: Iterative Analysis

```R
# 1. Run initial analysis
source("scripts/run_pipeline.R")

# 2. Check QC, identify issues
# - Remove outlier samples
# - Adjust filtering

# 3. Re-run relevant steps
source("scripts/02_quality_control.R")
source("scripts/03_differential_expression.R")

# 4. Compare results
# - Compare de_results_v1.csv vs de_results_v2.csv
```

---

## Next Steps

After completing the analysis:

1. **Biological Interpretation**
   - What processes are affected?
   - Are results consistent with hypothesis?
   - What are the top hits?

2. **Validation**
   - Select genes for qRT-PCR validation
   - Plan follow-up experiments

3. **Deeper Analysis**
   - Network analysis
   - Integration with other omics
   - Single-cell follow-up

4. **Publication**
   - Prepare figures
   - Write methods section
   - Share code and data

---

## Resources

### Documentation
- **DESeq2:** https://bioconductor.org/packages/DESeq2
- **clusterProfiler:** https://yulab-smu.top/biomedical-knowledge-mining-book/

### Learning
- DESeq2 vignette: `browseVignettes("DESeq2")`
- Bioconductor workflows: https://bioconductor.org/packages/release/BiocViews.html#___Workflow

### Community
- Bioconductor support: https://support.bioconductor.org/
- Stack Overflow: https://stackoverflow.com/questions/tagged/deseq2

### Citation
If you use this pipeline, please cite the key packages (see README.md).

---

**Happy Analyzing!** 🧬

