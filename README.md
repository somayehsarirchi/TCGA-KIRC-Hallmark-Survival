# TCGA-KIRC Hallmark-Based Survival Modeling Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18761184.svg)](https://doi.org/10.5281/zenodo.18761184)\
![R](https://img.shields.io/badge/R-%3E%3D4.2-blue)\
![Status](https://img.shields.io/badge/status-reproducible-success)\
![Validation](https://img.shields.io/badge/validation-bootstrap%20%2B%20train%2Ftest-orange)

------------------------------------------------------------------------

## Executive Summary

This repository implements a fully reproducible, modular, and
validation-aware translational bioinformatics pipeline for **TCGA Kidney
Renal Clear Cell Carcinoma (TCGA-KIRC)**.

The workflow progresses from raw RNA-seq counts to internally validated
multivariable survival modeling using Hallmark pathway activity scores.

Importantly, the pipeline distinguishes between:

-   Biological pathway enrichment
-   Independent prognostic contribution
-   Penalized multivariable selection

This enables mechanistic interpretation while preserving statistical
rigor.

------------------------------------------------------------------------

## Conceptual Framework

This pipeline integrates three analytical layers:

1.  **Biological Signal Discovery**
    -   Differential expression (DESeq2)
    -   Hallmark GSEA enrichment
2.  **Prognostic Signal Identification**
    -   Univariate Cox screening
    -   LASSO penalized regression
    -   Stability selection
3.  **Clinical Risk Modeling & Validation**
    -   Multivariable Cox modeling
    -   Bootstrap internal validation (B = 500)
    -   Time-dependent ROC
    -   RMS calibration
    -   Honest Train/Test validation

This separation prevents conflating enrichment significance with
independent predictive power.

------------------------------------------------------------------------

## Biological Rationale

Clear cell renal carcinoma (KIRC) is characterized by:

-   Hypoxia signaling activation
-   Metabolic reprogramming
-   Inflammatory pathway enrichment
-   EMT-related progression mechanisms

While multiple Hallmark pathways are enriched between molecular risk
groups, penalized modeling identifies a partially distinct subset as the
most stable independent prognostic features.

This highlights the difference between:

> Biological enrichment and independent multivariable prognostic
> contribution.

------------------------------------------------------------------------

## Complete Workflow

### 1. Data Acquisition

-   TCGA-KIRC RNA-seq (Primary Tumor)
-   Downloaded via `TCGAbiolinks`
-   STAR count workflow

### 2. Preprocessing

-   Clinical integration
-   Survival time construction
-   AJCC stage binarization (Early vs Late)
-   Tumor grade binarization (Low vs High, if available)
-   Variance Stabilizing Transformation (DESeq2 VST)

### 3. Hallmark Pathway Activity (ssGSEA)

-   MSigDB Hallmark gene sets
-   Sample-level pathway activity scoring
-   Pathway matrix used for downstream modeling

### 4. Differential Expression & GSEA

DESeq2 comparison between molecular risk groups:

-   `DESeq2_mgroup_adjusted_fullResults.csv`
-   `DESeq2_mgroup_significant_genes.csv`

GSEA on ranked genes:

-   `GSEA_Hallmark_mgroup.csv`
-   `PANEL_gsea_colored.png`

These analyses characterize tumor biology but are not directly used for
feature selection.

------------------------------------------------------------------------

### 5. Feature Selection Strategy

Stepwise selection process:

1.  Univariate Cox screening
2.  Top 30 Hallmarks retained
3.  LASSO penalized Cox regression (`glmnet`)
4.  Stability selection across resampling

Outputs:

-   `Cox_univariate_Hallmarks.csv`
-   `SelectedPathways_LASSO_beta.csv`
-   `VALID_lasso_selection_stability.csv`

A molecular risk score was constructed from LASSO-selected pathways.

------------------------------------------------------------------------

### 6. Multivariable Survival Model

Adjusted Cox model including:

-   Molecular risk score
-   Age
-   AJCC stage
-   Tumor grade (if available)

Outputs:

-   `Cox_final_multivariable_grade_adjusted.csv`
-   `Model_summaries.txt`
-   `KM_mgroup.png`

------------------------------------------------------------------------

### 7. Internal Validation

#### Time-Dependent ROC

-   `timeROC_AUC.csv`
-   `PANEL_timeROC.png`
-   `TEST_timeROC.png`

Observed AUC range: \~0.65--0.75

#### Bootstrap Validation (B = 500)

-   `VALID_bootstrap_grade_adjusted.csv`
-   `VALID_bootstrap_grade_adjusted.txt`

#### RMS Calibration

-   `PANEL_calibration_RMS_1y3y5y_legendBottom.png`

#### Proportional Hazards Diagnostics

-   `VALID_cox_zph.txt`
-   `VALID_cox_zph.png`

#### Honest Train/Test Validation

-   `TRAIN_TEST_summary.csv`
-   `TEST_KM_mgroup.png`

Feature selection was performed exclusively on the training set to
prevent data leakage.

------------------------------------------------------------------------

## Summary Panel

Integrated multi-layer visualization:

`FIGURE_SummaryPanel_Molecular.png`

Includes:

-   PCA visualization
-   Volcano plot
-   Hallmark GSEA barplot
-   Kaplan--Meier curve
-   Time-dependent ROC curves

------------------------------------------------------------------------

## Reproducibility & Architecture

-   Modular script structure
-   Cache-aware `.rds` storage
-   Deterministic seed control
-   sessionInfo exported
-   Force rebuild option

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

------------------------------------------------------------------------

## Key Strengths

-   Clear separation of enrichment vs prognostic modeling
-   Penalized regression (glmnet)
-   Stability selection
-   Bootstrap validation (B = 500)
-   Calibration assessment
-   Honest train/test split
-   Modular reproducible design
-   DOI archived release

------------------------------------------------------------------------

## Scientific Interpretation

The molecular risk score demonstrates strong association with overall
survival in the full cohort.

Performance attenuation in the held-out test set reflects transparent
modeling and avoids over-optimistic bias.

The divergence between enriched pathways and LASSO-selected pathways
underscores the distinction between tumor biology and independent
prognostic contribution.

------------------------------------------------------------------------

## Limitations

-   Internal validation only (no external cohort)
-   Based on TCGA bulk RNA-seq
-   Pathway-level abstraction may mask gene-level heterogeneity

------------------------------------------------------------------------

## Requirements

R ≥ 4.2

Core packages:

TCGAbiolinks\
DESeq2\
GSVA\
glmnet\
survival\
timeROC\
rms\
pec\
ggplot2\
dplyr\

------------------------------------------------------------------------

## Author

Somayeh Sarirchi\
Bioinformatics & Translational Cancer Research\
RNA-seq \| Survival Modeling \| Pathway Analysis \| Reproducible
Pipelines

------------------------------------------------------------------------

If using this workflow academically, please cite TCGA and MSigDB
resources accordingly.
