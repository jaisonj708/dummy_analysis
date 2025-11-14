# Quick Start Guide

This guide will help you run the bulk RNA-seq analysis pipeline in under 5 minutes.

## Prerequisites

- R version 4.0 or higher
- At least 2 GB of free disk space

## Installation (One-time setup)

### 1. Install R Packages

Open R or RStudio and run:

```R
source("requirements.R")
```

This will install all necessary packages. It may take 10-15 minutes depending on your internet connection.

**Note:** If you encounter any issues, you may need to install system dependencies:
- **macOS:** `brew install openssl libxml2`
- **Ubuntu/Debian:** `sudo apt-get install libssl-dev libxml2-dev libcurl4-openssl-dev`

## Running the Pipeline

### Option 1: Run Everything at Once (Recommended for first-time users)

```R
source("scripts/run_pipeline.R")
```

This will:
1. Generate dummy RNA-seq data (10,000 genes, 12 samples)
2. Perform quality control
3. Run differential expression analysis
4. Create visualizations
5. Perform enrichment analysis

**Expected runtime:** 2-5 minutes

### Option 2: Run Steps Individually

If you want more control or need to debug:

```R
# Step 1: Generate dummy data
source("scripts/01_generate_dummy_data.R")

# Step 2: Quality control
source("scripts/02_quality_control.R")

# Step 3: Differential expression
source("scripts/03_differential_expression.R")

# Step 4: Visualizations
source("scripts/04_visualization.R")

# Step 5: Enrichment analysis
source("scripts/05_enrichment_analysis.R")
```

## Viewing Results

After the pipeline completes, explore the results:

### 1. Quality Control Reports
```
results/qc/
├── library_sizes.pdf          # Total reads per sample
├── gene_detection.pdf         # Genes detected per sample
├── sample_correlation.pdf     # Sample similarity
├── pca_raw.pdf               # PCA of raw counts
├── dispersion_estimates.pdf  # DESeq2 dispersion plot
├── pvalue_distribution.pdf   # P-value histogram
└── fc_distribution.pdf       # Fold change distribution
```

### 2. Differential Expression Results
```
results/de/
├── de_results.csv            # Full results table
├── de_significant.csv        # Significant genes only
└── normalized_counts.csv     # Normalized expression matrix
```

### 3. Visualizations
```
results/plots/
├── pca.pdf                   # PCA of samples
├── volcano_plot.pdf          # Volcano plot
├── ma_plot.pdf              # MA plot
├── heatmap_top_genes.pdf    # Heatmap of top DE genes
├── expression_patterns.pdf   # Expression of individual genes
└── summary_barplot.pdf       # Summary statistics
```

### 4. Enrichment Analysis
```
results/enrichment/
├── go_enrichment_all.csv           # All GO terms
├── go_enrichment_upregulated.csv   # Upregulated genes
├── go_enrichment_downregulated.csv # Downregulated genes
├── go_barplot_upregulated.pdf     # Visualization
└── go_dotplot.pdf                 # Overview plot
```

## Key Results to Check

### 1. How many genes are differentially expressed?

Look at `results/plots/summary_barplot.pdf` or check the console output from step 3.

### 2. What are the top differentially expressed genes?

Open `results/de/de_significant.csv` and sort by adjusted p-value.

### 3. Do samples cluster by condition?

Check `results/plots/pca.pdf` - control and treated samples should separate.

### 4. What biological processes are affected?

Review `results/enrichment/go_barplot_upregulated.pdf` and the corresponding CSV files.

## Customization

To modify analysis parameters, edit `config.R`:

```R
# Example: Change significance thresholds
PADJ_THRESHOLD <- 0.01      # More stringent (default: 0.05)
LOG2FC_THRESHOLD <- 2       # Larger fold changes (default: 1)

# Example: Generate more/fewer genes
N_GENES <- 20000            # More genes (default: 10000)
N_SAMPLES_PER_CONDITION <- 3  # Fewer samples (default: 6)
```

After editing `config.R`, re-run the pipeline.

## Using Your Own Data

To analyze your own RNA-seq data:

1. Prepare your count matrix as `data/raw_counts.csv`:
   ```
   gene_id,Sample1,Sample2,Sample3,...
   GENE1,100,150,120,...
   GENE2,50,60,55,...
   ```

2. Prepare sample metadata as `data/sample_metadata.csv`:
   ```
   sample_id,condition,replicate,batch
   Sample1,Control,1,Batch1
   Sample2,Control,2,Batch1
   Sample3,Treated,1,Batch1
   ```

3. Edit `config.R` to skip data generation, or simply run:
   ```R
   # Skip step 1, run steps 2-5
   source("scripts/02_quality_control.R")
   source("scripts/03_differential_expression.R")
   source("scripts/04_visualization.R")
   source("scripts/05_enrichment_analysis.R")
   ```

## Troubleshooting

### Issue: "Package X is not available"
- Run `source("requirements.R")` again
- Some packages require Bioconductor: `BiocManager::install("PackageName")`

### Issue: "Cannot open file 'data/raw_counts.csv'"
- Make sure you ran step 1 first: `source("scripts/01_generate_dummy_data.R")`
- Or provide your own data files

### Issue: Enrichment analysis fails
- This is expected with dummy data (genes are simulated IDs)
- For real data, ensure gene IDs are in the correct format (e.g., Entrez IDs)

### Issue: Out of memory
- Reduce `N_GENES` in `config.R`
- Close other applications
- Use a machine with more RAM

## Getting Help

1. Check the main README.md for detailed documentation
2. Review individual script comments for specific questions
3. Consult package documentation:
   - DESeq2: https://bioconductor.org/packages/DESeq2
   - clusterProfiler: https://yulab-smu.top/biomedical-knowledge-mining-book/

## Citation

If you use this pipeline for your research, please cite the key packages:

- **DESeq2:** Love MI, Huber W, Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome Biology, 15, 550.

- **clusterProfiler:** Yu G, Wang LG, Han Y, He QY (2012). "clusterProfiler: an R package for comparing biological themes among gene clusters." OMICS, 16(5), 284-287.

---

**Happy analyzing!** 🧬📊

