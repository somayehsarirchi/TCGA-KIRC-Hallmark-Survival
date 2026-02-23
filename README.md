# TCGA-KIRC Hallmark Survival Modeling

[![R](https://img.shields.io/badge/language-R-blue.svg)]() [![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)]()

------------------------------------------------------------------------

## Overview

This repository provides a fully reproducible computational oncology
pipeline for molecular risk stratification in **TCGA Kidney Renal Clear
Cell Carcinoma (KIRC)**.

The workflow integrates:

-   Hallmark pathway activity scoring (ssGSEA)
-   Penalized Cox modeling (LASSO)
-   Multivariable clinical adjustment
-   Internal bootstrap validation
-   Time-dependent discrimination
-   Calibration analysis (rms)
-   Differential expression (DESeq2)
-   Hallmark GSEA
-   Automated multi-panel figure generation

The objective is to evaluate whether a Hallmark-derived molecular risk
score provides independent prognostic value beyond established
clinicopathological variables.

------------------------------------------------------------------------

## Objectives

-   Quantify Hallmark pathway activity using ssGSEA
-   Identify prognostic pathways via univariable Cox regression
-   Construct a molecular risk score using LASSO-regularized Cox
    modeling
-   Adjust for:
    -   Age at diagnosis
    -   AJCC stage
    -   Tumor grade
-   Evaluate:
    -   Discrimination (C-index, timeROC)
    -   Calibration (bootstrap-corrected)
    -   PH assumption
    -   LASSO stability
-   Characterize downstream tumor biology (DESeq2 + GSEA)

------------------------------------------------------------------------

## Dataset

-   Cohort: **TCGA-KIRC (Primary Tumor)**
-   RNA-Seq: **STAR -- Counts**
-   Clinical variables:
    -   Age
    -   AJCC stage
    -   Tumor grade
-   Data accessed via **TCGAbiolinks**

Final cohort size, number of events, and follow-up statistics are
automatically extracted into:

    results/README_numbers.txt

------------------------------------------------------------------------

## Repository Structure

    scripts/
      00_config.R
      01_download_preprocess.R
      02_ssgsea_scoring.R
      03_survival_modeling.R
      04_DEG_GSEA.R
      05_validation.R
      06_summary_panel.R
      07_rms_calibration.R
      07_extract_readme_numbers.R
      run_all.R

    data_cache/
    results/
      figures/
      tables/
      sessionInfo.txt
      README_numbers.txt

------------------------------------------------------------------------

## How to Run the Full Pipeline

### Option 1 --- From R

From project root:

``` r
source("scripts/run_all.R")
```

### Option 2 --- From command line

From project root:

``` bash
Rscript scripts/run_all.R
```

This will automatically generate:

-   All intermediate objects
-   All tables
-   All figures
-   Calibration analysis
-   Final multi-panel summary figure
-   Session information

------------------------------------------------------------------------

## Reproducibility

-   Fixed random seeds
-   Bootstrap validation implemented
-   All intermediate objects cached
-   sessionInfo saved
-   Deterministic pipeline
-   MIT License

------------------------------------------------------------------------

## Author

**Somayeh Sarirchi**\
Computational Oncology / Bioinformatics
