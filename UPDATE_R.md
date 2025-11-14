# Update R to Latest Version

## Current Status

- **Current R version:** 4.2.0 (from April 2022)
- **Latest R version:** 4.4.2 (November 2024)
- **Installer downloaded:** `/tmp/R-latest.pkg` (ready to install!)

## Why Update?

Your R 4.2.0 is over 2.5 years old. Updating will:
- ✅ Fix package installation issues (especially DESeq2)
- ✅ Improve performance and stability
- ✅ Get latest Bioconductor support
- ✅ Access newer package versions

## Quick Install (2 minutes)

### Option 1: Install Downloaded Package (Easiest)

I've already downloaded R 4.4.2 for you. Just run:

```bash
sudo installer -pkg /tmp/R-latest.pkg -target /
```

**Or** double-click on the file in Finder:
```bash
open /tmp/R-latest.pkg
```

Then follow the on-screen installer.

### Option 2: Manual Download

If the above doesn't work:

1. Go to: https://cloud.r-project.org/bin/macosx/
2. Download: **R-4.4.2-arm64.pkg** (for Apple Silicon)
3. Double-click to install
4. Follow installer prompts

## After Installing R

### 1. Verify Installation

```bash
R --version
# Should show: R version 4.4.2 (2024-10-31)
```

### 2. Reinstall Packages

```R
# Start R
R

# Install BiocManager first
install.packages("BiocManager")

# Run the requirements script
source("/Users/jaisonjain/del/bulkseq/requirements.R")
```

**This will take 10-15 minutes** but should succeed without errors.

### 3. Run the Pipeline

```bash
cd /Users/jaisonjain/del/bulkseq
Rscript scripts/run_pipeline.R
```

## Expected Results

After updating R and reinstalling packages:

- ✅ DESeq2 will install successfully
- ✅ clusterProfiler will install successfully
- ✅ All pipeline scripts will run
- ✅ You'll get ~25 output files in `results/`

## Alternative: Quick Test Without Updating

If you want to test the pipeline logic without updating R, I can modify the scripts to use `edgeR` instead of `DESeq2` (both are already installed). However, updating R is the better long-term solution.

## Troubleshooting

### Issue: "R not found" after update

```bash
# Check where R is installed
which R
# Should be: /Library/Frameworks/R.framework/Resources/bin/R

# If not, add to PATH:
echo 'export PATH="/Library/Frameworks/R.framework/Resources/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

### Issue: Old R still showing

```bash
# Restart terminal or run:
hash -r
R --version
```

### Issue: Package installation still fails

```bash
# Install system dependencies first:
brew install harfbuzz fribidi freetype libpng cairo

# Then reinstall packages
Rscript requirements.R
```

## Quick Command Summary

```bash
# 1. Install R (choose one):
sudo installer -pkg /tmp/R-latest.pkg -target /
# OR
open /tmp/R-latest.pkg

# 2. Verify
R --version

# 3. Reinstall packages
cd /Users/jaisonjain/del/bulkseq
Rscript requirements.R

# 4. Run pipeline
Rscript scripts/run_pipeline.R
```

## Timeline

- **Update R:** 2-3 minutes
- **Reinstall packages:** 10-15 minutes  
- **Run pipeline:** 3-5 minutes
- **Total:** ~20 minutes to complete everything

---

## Ready to Update?

The installer is ready at `/tmp/R-latest.pkg`. Just run:

```bash
sudo installer -pkg /tmp/R-latest.pkg -target /
```

Enter your password when prompted, and you'll have the latest R in ~2 minutes!

