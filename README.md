# Bulk RNA-seq Analysis Pipeline

A complete R-based pipeline for bulk RNA-seq data analysis, from quality control to differential expression and functional enrichment analysis.

## Project Structure

```
bulkseq/
├── README.md
├── requirements.R                    # R package dependencies
├── data/
│   ├── raw_counts.csv               # Raw count matrix
│   └── sample_metadata.csv          # Sample information
├── scripts/
│   ├── 01_generate_dummy_data.R     # Generate dummy data
│   ├── 02_quality_control.R         # QC and filtering
│   ├── 03_differential_expression.R # DE analysis with DESeq2
│   ├── 04_visualization.R           # Plots and figures
│   ├── 05_enrichment_analysis.R     # Functional enrichment
│   └── run_pipeline.R               # Master script
├── results/
│   ├── qc/                          # QC outputs
│   ├── de/                          # DE results
│   ├── plots/                       # Figures
│   └── enrichment/                  # Enrichment results
└── config.R                         # Configuration parameters
```

## Pipeline Overview

### Step 1: Data Generation
Creates dummy RNA-seq count data with:
- 10,000 genes
- 12 samples (6 control, 6 treated)
- Simulated differential expression

### Step 2: Quality Control
- Library size distribution
- Gene detection rates
- Sample correlation analysis
- Filtering low-count genes

### Step 3: Differential Expression Analysis
- Normalization with DESeq2
- Statistical testing for differential expression
- Multiple testing correction
- Export results tables

### Step 4: Visualization
- PCA plot for sample clustering
- Heatmap of top differentially expressed genes
- Volcano plot of DE results
- MA plot

### Step 5: Functional Enrichment
- GO term enrichment analysis
- Pathway analysis
- Visualization of enriched terms

## Installation

### 1. Install R (version >= 4.0)
Download from: https://www.r-project.org/

### 2. Install Required Packages
Run in R console:
```R
source("requirements.R")
```

Or install manually:
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "edgeR", "clusterProfiler", "org.Hs.eg.db", 
                       "EnhancedVolcano", "ComplexHeatmap"))
install.packages(c("ggplot2", "pheatmap", "RColorBrewer", "dplyr", "tidyr", 
                   "readr", "viridis", "ggrepel"))
```

## Usage

### Quick Start - Run Complete Pipeline
```R
source("scripts/run_pipeline.R")
```

### Run Individual Steps

#### 1. Generate Dummy Data
```R
source("scripts/01_generate_dummy_data.R")
```

#### 2. Quality Control
```R
source("scripts/02_quality_control.R")
```

#### 3. Differential Expression Analysis
```R
source("scripts/03_differential_expression.R")
```

#### 4. Generate Visualizations
```R
source("scripts/04_visualization.R")
```

#### 5. Enrichment Analysis
```R
source("scripts/05_enrichment_analysis.R")
```

## Configuration

Edit `config.R` to modify:
- Significance thresholds (p-value, log2FC)
- Number of genes/samples
- Plot parameters
- Output directories

## Output Files

### QC Results (`results/qc/`)
- `library_sizes.pdf` - Distribution of read counts per sample
- `gene_detection.pdf` - Number of detected genes per sample
- `sample_correlation.pdf` - Sample-to-sample correlation heatmap
- `pca_raw.pdf` - PCA of raw counts

### DE Results (`results/de/`)
- `de_results.csv` - Full differential expression results
- `de_significant.csv` - Filtered significant genes
- `normalized_counts.csv` - Normalized expression matrix

### Plots (`results/plots/`)
- `pca.pdf` - PCA of normalized counts
- `heatmap_top_genes.pdf` - Heatmap of top DE genes
- `volcano_plot.pdf` - Volcano plot
- `ma_plot.pdf` - MA plot

### Enrichment Results (`results/enrichment/`)
- `go_enrichment.csv` - GO term enrichment results
- `go_plot.pdf` - Visualization of enriched GO terms

## Methods

### Differential Expression
- **Tool**: DESeq2
- **Design**: ~ condition
- **Normalization**: median-of-ratios
- **Statistical test**: Wald test
- **Multiple testing correction**: Benjamini-Hochberg (FDR)

### Filtering Criteria
- Genes with < 10 total counts across all samples are removed
- Genes detected in < 3 samples are removed

### Significance Thresholds
- Adjusted p-value < 0.05
- |log2 Fold Change| > 1

## References

- Love MI, Huber W, Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." *Genome Biology*, 15, 550.
- Robinson MD, McCarthy DJ, Smyth GK (2010). "edgeR: a Bioconductor package for differential expression analysis of digital gene expression data." *Bioinformatics*, 26(1), 139-140.
- Yu G, Wang LG, Han Y, He QY (2012). "clusterProfiler: an R package for comparing biological themes among gene clusters." *OMICS*, 16(5), 284-287.

## Support

For issues or questions, please check:
1. R package documentation
2. Bioconductor support site
3. GitHub issues for specific packages

## License

MIT License

