# Alzheimer's Disease Biomarker Discovery — Microarray Differential Expression & ROC Analysis

A reproducible R/Bioconductor pipeline that identifies candidate gene-expression biomarkers for Alzheimer's disease from a public brain microarray dataset, using differential expression analysis and ROC/AUC-based biomarker evaluation.

> **Note:** this is a **microarray** study (Affymetrix platform, analyzed with `limma`), not RNA-seq. See "Before You Push" at the bottom for a repo-naming fix.

## Dataset

- **Source:** NCBI GEO accession [GSE5281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281)
- **Platform:** Affymetrix Human Genome U133 Plus 2.0 Array
- **Subset analyzed:** Entorhinal Cortex (EC) samples only, Control vs. Alzheimer's Disease (n = 21 samples after filtering for region and disease-state annotation)

## Methodology

1. **Data retrieval & QC** — expression matrix and sample metadata pulled via `GEOquery`; missing values checked before proceeding.
2. **Normalization** — quantile normalization (`normalizeBetweenArrays`), with before/after boxplots saved to confirm the correction.
3. **Low-expression filtering** — genes in the bottom 25% by mean expression are removed.
4. **Region/condition subsetting** — samples restricted to the Entorhinal Cortex region, split into Control vs. AD groups.
5. **Differential expression** — `limma` linear model (`~0 + group` design), empirical Bayes moderation (`eBayes`), contrast `AD - Control`.
6. **Significance filtering** — exploratory threshold (adj. P < 0.05, |log2FC| > 1) and a stricter biomarker-candidate threshold (adj. P < 0.01, |log2FC| > 2).
7. **Visualization** — volcano plot of all tested genes (colored by significance) and a heatmap of the top 50 differentially expressed genes (row-scaled, `pheatmap`).
8. **Probe-to-gene mapping** — significant probe IDs mapped to gene symbols using `hgu133plus2.db` / `AnnotationDbi`.
9. **Biomarker evaluation** — top candidate probes evaluated individually as classifiers (Control vs. AD) via ROC curves and AUC, using `pROC`; optimal sensitivity/specificity threshold reported for the best-performing probe.

## Repository Structure

```
├── Alzheimer-Microarray-DEG-analysis.R   # full analysis script
└── README.md
```

## Key Results

- 21 Entorhinal Cortex samples analyzed (Control vs. Alzheimer's Disease)
- Strict biomarker-candidate filter (adj. P < 0.01, |log2FC| > 2) applied on top of the exploratory DEG list
- Three candidate biomarker probes evaluated by ROC/AUC: `205914_s_at`, `214333_x_at`, `209364_at`
- _[Add your top gene symbols and the actual AUC value for your best-performing probe here, once you re-run the script — see "Before You Push" below]_

## Tools & Libraries

| Category | Tools |
|---|---|
| Data retrieval | `GEOquery` |
| Differential expression | `limma` |
| Visualization | `ggplot2` (volcano plot), `pheatmap` (heatmap) |
| Probe annotation | `hgu133plus2.db`, `AnnotationDbi` |
| Biomarker evaluation | `pROC` (ROC curves, AUC) |

## How to Reproduce

```r
install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "hgu133plus2.db", "AnnotationDbi"))
install.packages(c("dplyr", "ggplot2", "pheatmap", "pROC"))

source("Alzheimer-Microarray-DEG-analysis.R")
```

## Limitations & Future Work

- Single brain region (Entorhinal Cortex) and a small sample size (n = 21); findings are exploratory and not adjusted for multiple-region or multi-cohort validation.
- Candidate biomarkers are evaluated on the same dataset used to discover them — independent validation on a separate GEO cohort (e.g., GSE1297 or GSE36980) would be the natural next step.
- No correction for potential confounders (age, sex) was applied in the current design matrix.

---

