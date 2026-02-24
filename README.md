# TCGA-KIRC Hallmark-Based Survival Modeling Pipeline
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18761184.svg)](https://doi.org/10.5281/zenodo.18761184)
![R](https://img.shields.io/badge/R-%3E%3D4.2-blue)
![Status](https://img.shields.io/badge/status-reproducible-success)
![Validation](https://img.shields.io/badge/validation-bootstrap%20%2B%20train%2Ftest-orange)

------------------------------------------------------------------------

## Overview

This repository implements a fully reproducible end-to-end
bioinformatics pipeline for **TCGA Kidney Renal Clear Cell Carcinoma
(TCGA-KIRC)** integrating:

-   RNA-seq data processing
-   Hallmark pathway activity scoring (ssGSEA)
-   Differential expression analysis (DESeq2)
-   Hallmark GSEA enrichment
-   LASSO-based pathway selection
-   Multivariable Cox regression modeling
-   Bootstrap validation (B = 500)
-   Time-dependent ROC analysis
-   RMS calibration
-   Train/Test validation without data leakage

The goal is to demonstrate a rigorous translational bioinformatics
workflow from raw transcriptomic data to validated survival modeling.

------------------------------------------------------------------------

## Biological Rationale

Clear cell renal carcinoma (KIRC) is characterized by:

-   Hypoxia signaling activation
-   Metabolic reprogramming
-   Inflammatory pathway enrichment
-   EMT-related progression mechanisms

This pipeline models survival using **Hallmark pathway-level activity**
rather than individual genes, improving interpretability and biological
coherence.

------------------------------------------------------------------------

## Complete Workflow

### 1. Data Acquisition

-   TCGA-KIRC RNA-seq (Primary Tumor)
-   Downloaded via `TCGAbiolinks`
-   STAR count workflow

### 2. Preprocessing

-   Clinical integration
-   Survival variable construction
-   AJCC stage binarization (Early vs Late)
-   Tumor grade binarization (Low vs High, if available)

### 3. Hallmark Pathway Scoring

-   Variance Stabilizing Transformation (VST)
-   Hallmark gene sets (MSigDB)
-   ssGSEA pathway activity scores per sample

### 4. Differential Expression (DESeq2)

Differential expression between molecular risk groups:

-   `DESeq2_mgroup_adjusted_fullResults.csv`
-   `DESeq2_mgroup_significant_genes.csv`

These results support biological differences between High and Low
molecular risk groups.

Volcano plot: - `PANEL_volcano_labeled.png`

### 5. Hallmark Gene Set Enrichment

GSEA performed on DESeq2-ranked genes:

-   `GSEA_Hallmark_mgroup.csv`
-   `PANEL_gsea_colored.png`

Enriched pathways include hypoxia, EMT, inflammatory response, and
metabolic signaling --- consistent with KIRC biology.

### 6. Feature Selection & Risk Score

-   Univariate Cox screening: `Cox_univariate_Hallmarks.csv`
-   Top 30 hallmarks: `VALID_top30_univariate_hallmarks.csv`
-   LASSO selection: `SelectedPathways_LASSO_beta.csv`
-   Stability assessment: `VALID_lasso_selection_stability.csv`

A molecular risk score was constructed from selected pathways.

### 7. Multivariable Survival Modeling

Adjusted Cox model including:

-   Molecular risk score
-   Age
-   AJCC stage
-   Tumor grade (if available)

Outputs: - `Cox_final_multivariable_grade_adjusted.csv` -
`Model_summaries.txt`

Kaplan--Meier curve: - `KM_mgroup.png`

### 8. Model Validation

#### Time-dependent ROC

-   `timeROC_AUC.csv`
-   `PANEL_timeROC.png`
-   `TEST_timeROC.png`

Observed AUC range: \~0.65--0.75 depending on time horizon.

#### RMS Calibration (Bootstrap B = 500)

-   `PANEL_calibration_RMS_1y3y5y_legendBottom.png`
-   `VALID_bootstrap_grade_adjusted.csv`
-   `VALID_bootstrap_grade_adjusted.txt`

Calibration demonstrates reasonable agreement between predicted and
observed survival probabilities.

#### Proportional Hazards Diagnostics

-   `VALID_cox_zph.txt`
-   `VALID_cox_zph.png`

#### Train/Test Validation

-   `TRAIN_TEST_summary.csv`
-   `TEST_KM_mgroup.png`

Feature selection was performed exclusively on the training set to
prevent data leakage.

------------------------------------------------------------------------

## Summary Panel

Integrated molecular summary figure:

-   `FIGURE_SummaryPanel_Molecular.png`

Includes: - PCA visualization - Volcano plot - Hallmark GSEA barplot -
Kaplan--Meier survival curve - Time-dependent ROC curves

------------------------------------------------------------------------

## Reproducibility

The pipeline is modular and cache-aware:

-   Intermediate objects stored as `.rds`
-   Downstream steps load only required inputs
-   Steps can be forced to rebuild with `force=TRUE`
-   `sessionInfo.txt` exported for environment tracking

### Run From Scratch

``` r
source("run_all.R")

main(
  force = TRUE,
  use_pec_calibration = TRUE,
  use_rms_calibration = TRUE,
  run_train_test = TRUE
)
```

### Fast Re-run

``` r
source("run_all.R")
main()
```

------------------------------------------------------------------------

## Key Strengths

-   End-to-end RNA-seq to survival modeling pipeline
-   Integration of molecular and clinical covariates
-   Penalized regression (glmnet)
-   Bootstrap validation (B = 500)
-   Calibration assessment
-   Honest train/test validation
-   Publication-ready visual outputs
-   Fully reproducible modular architecture

------------------------------------------------------------------------

## Scientific Interpretation

The molecular risk score shows strong association with overall survival
in the full cohort.

Discrimination decreases in the held-out test set, emphasizing the need
for external validation --- reflecting transparent and honest modeling
practices rather than overfitting.

------------------------------------------------------------------------

## Requirements

R ≥ 4.2

Core packages:

-   TCGAbiolinks
-   DESeq2
-   GSVA
-   glmnet
-   survival
-   survminer
-   timeROC
-   rms
-   pec
-   ggplot2
-   dplyr

------------------------------------------------------------------------

## Methodological Strengths

- No data leakage (feature selection restricted to training set)
- Bootstrap validation (B = 500)
- Calibration assessment (RMS)
- Proportional hazards diagnostics
- Penalized regression (glmnet)
- Modular cache-aware design

------------------------------------------------------------------------
## Limitations

- Internal validation only (no external cohort)
- Based on TCGA bulk RNA-seq
- Pathway-level abstraction may mask gene-level effects

------------------------------------------------------------------------
## Author
Somayeh Sarirchi\
Bioinformatics & Translational Cancer Research\
RNA-seq \| Survival Modeling \| Pathway Analysis \| Reproducible
Pipelines
------------------------------------------------------------------------

If using this workflow structure academically, please cite TCGA and
MSigDB resources accordingly.
