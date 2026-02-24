# TCGA-KIRC Hallmark-Based Survival Modeling Pipeline

![R](https://img.shields.io/badge/R-%3E%3D4.2-blue)
![Status](https://img.shields.io/badge/status-reproducible-success)
![Pipeline](https://img.shields.io/badge/workflow-modular-green)
![Validation](https://img.shields.io/badge/validation-bootstrap%20%2B%20train%2Ftest-orange)

------------------------------------------------------------------------

## Overview

This repository implements a fully reproducible, modular end-to-end
bioinformatics pipeline for survival modeling in **TCGA Kidney Renal
Clear Cell Carcinoma (TCGA-KIRC)** using Hallmark pathway activity
scores.

The project demonstrates:

-   Robust survival modeling
-   Internal bootstrap validation (B = 500)
-   Time-dependent ROC analysis
-   RMS calibration curves
-   Train/Test validation without data leakage
-   LASSO feature selection stability assessment
-   Publication-ready visualization
-   Cache-aware modular workflow design

This pipeline prioritizes **methodological rigor and transparent
reporting** over over-optimized predictive performance.

------------------------------------------------------------------------

## Why This Project Matters

Kidney Renal Clear Cell Carcinoma (KIRC) is characterized by:

-   Hypoxia signaling
-   EMT activation
-   Inflammatory pathways
-   Metabolic reprogramming

Instead of gene-level modeling, this pipeline leverages **Hallmark
pathway-level scoring (ssGSEA)** to build interpretable and biologically
meaningful survival models.

The project demonstrates how to:

-   Integrate transcriptomic data with clinical variables
-   Apply penalized Cox regression for feature selection
-   Perform rigorous internal validation
-   Assess calibration and discrimination
-   Avoid data leakage in model evaluation

This is a translational bioinformatics-oriented modeling framework
suitable for real-world research settings.

------------------------------------------------------------------------

## Project Structure

    TCGA-KIRC-Hallmark-Survival/
    │
    ├── R/
    │   ├── 00_utils.R
    │   ├── 01_gdc.R
    │   ├── 02_preprocess.R
    │   ├── 03_vst_symbol.R
    │   ├── 04_ssgsea.R
    │   ├── 05_uni_cox.R
    │   ├── 06_lasso_risk.R
    │   ├── 07_deseq2_gsea.R
    │   ├── 08_eval_panels.R
    │   └── 09_train_test.R
    │
    ├── cache/              # Auto-generated intermediate RDS objects
    ├── results/
    │   ├── figures/        # Publication-ready figures
    │   ├── tables/         # CSV summaries
    │   └── objects/        # Final model objects
    │
    ├── run_all.R
    └── README.md

------------------------------------------------------------------------

## Pipeline Workflow

### 1. Data Acquisition

-   TCGA-KIRC RNA-seq (Primary Tumor)
-   Retrieved using TCGAbiolinks
-   STAR count workflow

### 2. Preprocessing

-   Clinical filtering
-   Survival variable construction
-   AJCC stage binarization (Early vs Late)
-   Optional tumor grade binarization (Low vs High)

### 3. Pathway Scoring

-   VST normalization
-   Hallmark gene sets (MSigDB)
-   ssGSEA pathway activity scores

### 4. Feature Selection

-   Univariate Cox screening
-   Top 5 pathways by FDR
-   LASSO Cox regression
-   Molecular risk score construction

### 5. Multivariable Modeling

Cox model including: - Molecular risk score - Age - AJCC stage - Tumor
grade (if available)

### 6. Validation Strategy

-   Bootstrap internal validation (B = 500)
-   RMS calibration curves (1, 3, 5 years)
-   Time-dependent ROC curves
-   Proportional hazards diagnostics
-   LASSO selection stability
-   70/30 Train/Test split (no leakage)

------------------------------------------------------------------------

## Key Outputs

### Figures

-   Workflow diagram
-   Summary molecular panel
-   Kaplan-Meier curves
-   Time-dependent ROC curves
-   RMS calibration plots
-   Train/Test evaluation
-   PH diagnostics

### Tables

-   Multivariable hazard ratios
-   Bootstrap validation metrics
-   ROC AUC values
-   LASSO stability frequencies
-   Top univariate pathways
-   sessionInfo for reproducibility

------------------------------------------------------------------------

## How to Run

### First Run (Full rebuild)

``` r
source("run_all.R")

main(
  force = TRUE,
  use_pec_calibration = TRUE,
  use_rms_calibration = TRUE,
  run_train_test = TRUE
)
```

### Fast Re-run (Uses cache)

``` r
source("run_all.R")
main()
```

------------------------------------------------------------------------

## Reproducibility Features

-   Modular step-wise execution
-   Cache-aware design
-   Explicit RDS-based dependencies
-   Bootstrap validation
-   Transparent reporting of limitations
-   sessionInfo export

------------------------------------------------------------------------

## Validation Summary

-   Strong association in full cohort KM analysis
-   Moderate discrimination (AUC \~0.65--0.75)
-   Calibration curves demonstrate reasonable agreement
-   Test-set discrimination decreases (expected in honest validation)

This emphasizes robustness and transparency rather than overfitting.

------------------------------------------------------------------------

## For Principal Investigators

This repository demonstrates:

-   Advanced survival modeling in R
-   Penalized regression expertise (glmnet)
-   Pathway-level bioinformatics analysis
-   Internal validation best practices
-   Reproducible research engineering
-   Clean modular project architecture

This project is intended as a demonstration of translational
bioinformatics competency suitable for research engineer or
computational oncology positions.

------------------------------------------------------------------------

## Requirements

R ≥ 4.2

Key packages: - TCGAbiolinks - DESeq2 - GSVA - glmnet - survival -
survminer - timeROC - rms - pec - ggplot2 - dplyr

------------------------------------------------------------------------

## Author
Somayeh Sarirchi
Bioinformatics & Translational Cancer Research\
RNA-seq \| Survival Modeling \| Pathway Analysis \| Reproducible
Pipelines

------------------------------------------------------------------------

If you use this workflow structure in academic work, please cite TCGA
and MSigDB resources accordingly.
