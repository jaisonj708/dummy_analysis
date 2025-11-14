# Quality Control Guide for Bulk RNA-seq

## Understanding Quality Control Metrics

Quality control is **critical** in RNA-seq analysis. Poor quality data will produce unreliable results, no matter how sophisticated your statistical methods are. This guide explains each QC metric, what to look for, and how to interpret the results.

---

## Overview: Why QC Matters

**The Problem:** RNA-seq experiments can fail at many stages:

- Poor RNA quality during extraction
- Library preparation issues
- Sequencing errors
- Sample mix-ups or contamination
- Batch effects
- Outlier samples

**The Solution:** QC checks help you:

- ✓ Identify problematic samples before analysis
- ✓ Detect technical issues that need correction
- ✓ Ensure biological signal exceeds technical noise
- ✓ Validate that samples group as expected
- ✓ Decide which genes to include in analysis

---

## QC Feature 1: Library Size Distribution

### What It Is

**Library size** = total number of reads (counts) mapped to genes for each sample.

The plot shows a bar chart with each sample's total read count.

### Why It's Necessary

- **Sequencing depth matters:** More reads = better detection of lowly expressed genes
- **Consistency is key:** Large differences between samples can bias results
- **Quality indicator:** Very low library sizes may indicate failed libraries or poor sequencing

### What to Look For

#### ✅ REASSURING (Good):

```
Library sizes are:
- Similar across all samples (within 2-3 fold of each other)
- At least 10-20 million reads per sample
- No obvious outliers
- Relatively uniform within each condition
```

**Example:** All samples between 15-30 million reads

#### ⚠️ CONCERNING (Bad):

```
Problems include:
- Huge variation (some samples <5M, others >50M)
- Very low counts (<5 million) in some samples
- One or two extreme outliers
- Systematic differences between conditions
```

**Example:** Control samples all ~30M, Treated samples all ~5M

### What It Means

**If library sizes vary greatly:**

- DESeq2 normalization will handle moderate differences
- But extreme differences suggest technical problems
- May need to sequence some samples deeper
- Or remove failed samples from analysis

**Action:** If a sample has <20% of the median library size, consider removing it.

---

## QC Feature 2: Gene Detection Rates

### What It Is

**Gene detection** = number of genes with at least 1 read (count > 0) in each sample.

The plot shows how many genes are detected per sample.

### Why It's Necessary

- **RNA quality indicator:** Degraded RNA shows fewer detected genes
- **Library complexity:** Good libraries capture diverse transcripts
- **Technical success:** Failed libraries detect fewer genes
- **Consistency check:** Samples should detect similar numbers of genes

### What to Look For

#### ✅ REASSURING:

```
All samples show:
- 15,000-20,000+ genes detected (for human/mouse)
- Similar numbers across all samples
- No dramatic outliers
- Consistent within experimental groups
```

**Example:** All samples detect 16,000-18,000 genes

#### ⚠️ CONCERNING:

```
Warning signs:
- Large variation between samples (e.g., 10,000 vs 20,000)
- Very low detection (<10,000 genes)
- One or two samples detect far fewer genes
- Systematic difference between conditions
```

**Example:** Control samples detect 18,000 genes, one treated sample detects only 8,000

### What It Means

**Low gene detection suggests:**

- RNA degradation (low RIN score)
- Failed library preparation
- Very low sequencing depth
- Sample contamination

**Action:** Remove samples detecting <50% of the median number of genes.

---

## QC Feature 3: Sample-to-Sample Correlation

### What It Is

**Pearson correlation** between all pairs of samples based on log-transformed counts.

The heatmap shows how similar each sample's expression profile is to every other sample.

### Why It's Necessary

- **Biological replicates should correlate:** Samples from the same condition should be similar
- **Detects outliers:** Samples that don't correlate with their group are problematic
- **Identifies sample swaps:** Wrong labels will show up here
- **Batch effects:** Can reveal systematic technical differences

### What to Look For

#### ✅ REASSURING:

```
Correlation matrix shows:
- High correlation within groups (r > 0.95)
- Lower correlation between groups (r = 0.85-0.95)
- Clear block structure (replicates cluster together)
- All correlations positive and strong (r > 0.80)
```

**Visual pattern:**

```
        C1   C2   C3   T1   T2   T3
C1    [1.00][0.98][0.97][0.90][0.89][0.91]
C2    [0.98][1.00][0.98][0.88][0.87][0.90]
C3    [0.97][0.98][1.00][0.89][0.88][0.91]
T1    [0.90][0.88][0.89][1.00][0.97][0.98]
T2    [0.89][0.87][0.88][0.97][1.00][0.97]
T3    [0.91][0.90][0.91][0.98][0.97][1.00]
```

Controls cluster, treated cluster, but groups are distinguishable.

#### ⚠️ CONCERNING:

```
Warning signs:
- One replicate correlates poorly with its group (r < 0.85)
- Sample correlates better with wrong group
- Very low overall correlations (r < 0.80)
- No clear block structure (no clustering by condition)
- Negative correlations (very bad!)
```

**Example:** Control_3 correlates better with Treated samples than with other Controls → likely sample swap

### What It Means

**Low correlation within replicates:**

- Biological variability too high
- Technical issues with specific samples
- Sample mislabeling
- Contamination

**High correlation between all samples:**

- Very small treatment effect
- Experimental design issues
- Or just samples are very similar (could be real)

**Action:**

- Investigate samples with r < 0.85 to their group
- Check for sample swaps if patterns look wrong
- Consider removing outliers

---

## QC Feature 4: PCA (Principal Component Analysis)

### What It Is

**PCA** reduces the thousands of gene expression values to 2-3 main axes of variation.

The plot shows each sample as a point in 2D space, positioned by gene expression similarity.

### Why It's Necessary

- **Most important QC check!** Reveals overall data structure
- **Validates experimental design:** Conditions should separate
- **Detects batch effects:** Samples may cluster by batch instead of biology
- **Identifies outliers:** Samples far from their group
- **Quality indicator:** Good data has clear patterns

### What to Look For

#### ✅ REASSURING:

```
PCA plot shows:
- Clear separation between experimental groups
- Replicates cluster tightly together
- PC1 explains the treatment effect (30-60% variance)
- No obvious outliers far from their group
- Batch effects minimal or on PC2/PC3
```

**Visual example:**

```
      PC1 (50% variance)
PC2    ●●●         ▲▲▲
(20%)  ●●●         ▲▲▲
       Controls   Treated
```

Clear left-right separation by treatment.

#### ⚠️ CONCERNING:

```
Warning signs:
- No separation between groups (overlap completely)
- One or more samples far from their group
- Separation by batch, not by treatment
- Very low variance explained (<10% on PC1)
- Unexpected clustering patterns
```

**Example bad PCA:**

```
      Batch1          Batch2
PC1    ●▲●▲           ●▲●▲
PC2    ●▲●▲           ●▲●▲
```

Samples cluster by batch, not by treatment!

### What It Means

**No separation between groups:**

- Treatment had no effect (or very small effect)
- Wrong samples were compared
- Biological variability >> treatment effect
- Might still find DE genes, but pattern unclear

**Batch effects dominate:**

- Technical variation > biological variation
- Need batch correction (include batch in model)
- May obscure real biological signal

**Outliers:**

- Sample-specific technical issues
- Contamination
- Mislabeling
- Consider removing

**Action:**

- Remove outliers >3 SD from group mean
- Include batch in design formula if batch effects present
- If no separation, consider if experiment worked

---

## QC Feature 5: Dispersion Estimates (DESeq2 Diagnostic)

### What It Is

**Dispersion** = gene-to-gene variability in expression levels.

The plot shows:

- X-axis: mean expression level
- Y-axis: dispersion (variance measure)
- Black points: gene-wise estimates
- Red points: fitted trend
- Blue circles: final shrunken estimates

### Why It's Necessary

- **Statistical modeling check:** DESeq2 estimates variance for each gene
- **Model fit quality:** Ensures statistical tests are valid
- **Outlier detection:** Genes with unusual variance
- **Low replicate adjustment:** Borrowing information across genes

### What to Look For

#### ✅ REASSURING:

```
Dispersion plot shows:
- Clear decreasing trend (high expression → low dispersion)
- Most genes follow the red trend line
- Blue circles shrink toward red line
- Dispersion values mostly between 0.01-0.5
- Smooth curve, not jagged
```

**Pattern:** Classic "dispersing cloud" that narrows with higher expression

#### ⚠️ CONCERNING:

```
Warning signs:
- No clear trend (cloud is uniform)
- Dispersions all very high (>1.0)
- Many genes far from trend line
- Weird patterns or bumps in the curve
- All dispersions nearly identical
```

### What It Means

**Good dispersion plot = valid statistical tests**

**High dispersions everywhere:**

- High biological variability
- Small sample size
- Technical variation
- Tests will be conservative (fewer significant genes)

**No trend:**

- Something wrong with the data
- Possible filtering issue
- May need different approach

**Action:** Usually just diagnostic - if it looks weird, check your data and filtering.

---

## QC Feature 6: P-value Distribution

### What It Is

**Histogram** of raw (unadjusted) p-values from the differential expression test.

### Why It's Necessary

- **Test validity check:** Ensures statistical tests worked properly
- **Signal detection:** Reveals if there are real differences
- **Multiple testing context:** Understand FDR correction impact

### What to Look For

#### ✅ REASSURING (Strong signal):

```
P-value histogram shows:
- Peak at very small p-values (0-0.05)
- Relatively flat for p-values 0.1-1.0
- More small p-values than expected by chance
```

**Shape:** Spike at left (many small p-values) + flat elsewhere

This means: **Many truly differentially expressed genes detected!**

#### ✅ REASSURING (No signal - but valid):

```
Alternative good pattern:
- Uniform distribution (flat across 0-1)
- No peak anywhere
```

**Shape:** Completely flat

This means: **No differential expression detected (could be real biology)**

#### ⚠️ CONCERNING:

```
Warning signs:
- Peak near p = 1.0 (many high p-values)
- U-shape (peaks at both ends)
- Peak in the middle (p = 0.5)
- Very few small p-values
```

**Bad shapes indicate:**

- Statistical test assumptions violated
- Model doesn't fit data
- Need different approach

### What It Means

**Ideal:** Spike at left = many DE genes, worth proceeding

**Flat:** No DE genes, but tests are valid (might be real)

**Weird shapes:** Something wrong - check data quality, model specification

**Action:** If distribution looks bad, revisit QC steps, consider different normalization or test.

---

## QC Feature 7: Fold Change Distribution

### What It Is

**Histogram** of log2 fold changes for all genes, colored by significance.

### Why It's Necessary

- **Effect size overview:** How big are the expression changes?
- **Validation check:** Expect both up and down regulation
- **Threshold setting:** Helps choose appropriate log2FC cutoff
- **Biological plausibility:** Very large FCs may indicate artifacts

### What to Look For

#### ✅ REASSURING:

```
Fold change distribution shows:
- Centered near zero (most genes unchanged)
- Symmetrical (similar up and down)
- Some genes with large fold changes (±2-4)
- Significant genes scattered across range
- Tail on both sides (not all changes in one direction)
```

**Shape:** Bell curve centered at 0, with significant genes in the tails

#### ⚠️ CONCERNING:

```
Warning signs:
- Heavily skewed (all changes in one direction)
- Very wide distribution (many huge fold changes)
- All significant genes have tiny fold changes
- Bimodal (two peaks)
```

### What It Means

**Symmetrical, centered at zero:** Normal, expected pattern

**All changes up (or down):**

- Global scaling issue
- RNA degradation in one group
- Check normalization

**Tiny fold changes:**

- Small biological effect
- Statistically significant but not biologically relevant
- Consider raising fold change threshold

**Huge fold changes:**

- Possibly real (very dramatic response)
- Or artifacts (check those genes specifically)

**Action:** Use this to set meaningful fold change thresholds (typically |log2FC| > 1).

---

## QC Feature 8: Mean-Variance Relationship

### What It Is

**Scatter plot** showing mean expression vs. variance for each gene.

- Kept genes (green): Pass filtering criteria
- Removed genes (red): Failed filtering

### Why It's Necessary

- **Filtering validation:** Ensures filtering was appropriate
- **Heteroscedasticity check:** Variance should increase with mean
- **Quality assessment:** Identifies genes with unusual properties

### What to Look For

#### ✅ REASSURING:

```
Mean-variance plot shows:
- Clear positive relationship (higher mean → higher variance)
- Red points (filtered) cluster at low mean/variance
- Green points (kept) span the range
- No weird outliers far off the trend
```

**Pattern:** Positive correlation on log-log scale, with filtered genes at bottom left

#### ⚠️ CONCERNING:

```
Warning signs:
- No relationship between mean and variance
- Many high-mean genes filtered out
- Huge outliers kept in analysis
- All variance similar regardless of expression
```

### What It Means

**Good pattern:** Filtering worked correctly, removing uninformative genes while keeping expressed genes

**Too stringent filtering:** Many mid-expression genes removed

**Too lenient:** Keeping lots of near-zero genes (will add noise)

**Action:** Adjust filtering parameters in config.R if needed.

---

## Comprehensive QC Checklist

Run through this checklist for every RNA-seq dataset:

### ✅ Sample Quality

- [ ] Library sizes > 10M reads per sample
- [ ] Library sizes within 3-fold of each other
- [ ] Gene detection > 15,000 genes per sample
- [ ] No samples with dramatically fewer detected genes

### ✅ Sample Relationships

- [ ] Within-group correlation > 0.95
- [ ] Between-group correlation 0.85-0.95
- [ ] No obvious outliers (r < 0.80)
- [ ] No evidence of sample swaps

### ✅ Experimental Structure

- [ ] PCA shows separation between groups
- [ ] Biological replicates cluster together
- [ ] PC1 or PC2 captures treatment effect
- [ ] No extreme outliers in PCA
- [ ] Batch effects minimal or controlled

### ✅ Statistical Modeling

- [ ] Dispersion plot shows decreasing trend
- [ ] Dispersions reasonable (mostly 0.01-0.5)
- [ ] P-value distribution looks valid
- [ ] Fold changes centered at zero
- [ ] Mean-variance relationship evident

### ✅ Results Validation

- [ ] Reasonable number of DE genes found (not 0, not 100%)
- [ ] Both up and down regulated genes
- [ ] Effect sizes (fold changes) make biological sense
- [ ] Known positive controls are significant

---

## Decision Tree: What to Do Based on QC Results

### Scenario 1: All QC Looks Good ✅

**Action:** Proceed with confidence!

- Continue to DE analysis
- Trust your results
- Generate publication figures

### Scenario 2: One Outlier Sample ⚠️

**Action:** Remove and reanalyze

- Drop the outlier
- Re-run from QC step
- Compare results with/without

### Scenario 3: Batch Effects Visible 🔧

**Action:** Add batch to model

- Include batch in design formula: `~ batch + condition`
- Or use batch correction (removeBatchEffect for visualization)
- Check if results change

### Scenario 4: No Separation on PCA ⚠️

**Action:** Investigate thoroughly

- Check sample labels (possible mix-up)
- Review experimental design
- Consider if treatment simply had no effect
- Look for DE genes anyway (may still find some)
- Plan validation experiments

### Scenario 5: Multiple Issues ❌

**Action:** Stop and troubleshoot

- Don't proceed to DE analysis yet
- Check raw data quality (FASTQ QC)
- Review lab notebooks for issues
- Consider re-sequencing
- Consult with bioinformatics expert

---

## Summary: The Big Picture

**Goal of QC:** Ensure your data is trustworthy before drawing biological conclusions

**Key Principle:** Technical variation should be much smaller than biological signal

**Red Flags:**

1. Samples don't correlate with their replicates
2. No separation between experimental groups on PCA
3. Outlier samples far from their group
4. Batch effects dominate the analysis
5. Statistical tests produce weird p-value distributions

**Green Lights:**

1. Replicates cluster tightly together
2. Clear separation between conditions
3. Consistent library sizes and gene detection
4. Proper statistical model fit
5. Reasonable number of DE genes with interpretable fold changes

**Remember:** It's better to remove a problematic sample than to include it and get misleading results!

---

## Further Reading

- **DESeq2 vignette:** Detailed explanation of dispersion estimation
- **PCA interpretation:** https://liorpachter.wordpress.com/2014/05/26/what-is-principal-component-analysis/
- **RNA-seq QC:** FastQC documentation for raw read QC
- **Batch effects:** Johnson et al. (2007) "Adjusting batch effects in microarray expression data using empirical Bayes methods"

---

**Bottom Line:** Good QC = Trustworthy Results. If your QC checks pass, you can be confident in your differential expression analysis! 🧬✅
