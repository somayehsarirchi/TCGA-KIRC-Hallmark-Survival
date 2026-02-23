# TCGA-KIRC Hallmark Survival Modeling

[![R](https://img.shields.io/badge/language-R-blue.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)]()

---

## Overview

This repository provides a fully reproducible computational oncology pipeline for molecular risk stratification in **TCGA Kidney Renal Clear Cell Carcinoma (KIRC)**.

The workflow integrates:

- Hallmark pathway activity scoring (ssGSEA)
- Penalized Cox modeling (LASSO)
- Multivariable clinical adjustment
- Internal bootstrap validation
- Time-dependent discrimination
- Calibration analysis (rms)
- Differential expression (DESeq2)
- Hallmark GSEA
- Automated figure panel generation

The objective is to evaluate whether a Hallmark-derived molecular risk score provides independent prognostic value beyond established clinicopathological variables.

---

## Objectives

- Quantify Hallmark pathway activity using ssGSEA
- Identify prognostic pathways via univariable Cox regression
- Construct a molecular risk score using LASSO-regularized Cox modeling
- Adjust for:
  - Age at diagnosis
  - AJCC stage
  - Tumor grade
- Evaluate:
  - Discrimination (C-index, timeROC)
  - Calibration (bootstrap-corrected)
  - PH assumption
  - LASSO stability
- Characterize downstream tumor biology (DESeq2 + GSEA)

---

## Dataset

- Cohort: **TCGA-KIRC (Primary Tumor)**
- RNA-Seq: **STAR – Counts**
- Clinical variables:
  - Age
  - AJCC stage
  - Tumor grade
- Data accessed via **TCGAbiolinks**

Final cohort size, events, and follow-up statistics are automatically extracted into:
results/README_numbers.txt

---

## Methods Summary

### 1. Pathway Scoring
- Hallmark gene sets (MSigDB category H)
- ssGSEA scoring per sample

### 2. Survival Modeling
- Univariable Cox across all Hallmarks
- LASSO Cox on top-ranked pathways
- Risk score = linear predictor of selected pathways
- Risk groups = median split (High vs Low)

### 3. Multivariable Adjustment
Two model types (when available):

- **Binary model**
  - Stage (Early vs Late)
  - Grade (Low vs High)

- **Multi-level model**
  - AJCC 4-level stage
  - Grade G1–G4

### 4. Internal Validation
- Bootstrap validation (rms::validate)
- Optimism-corrected discrimination
- Calibration at 3 years
- PH assumption test (cox.zph)
- LASSO feature stability via bootstrap

### 5. Biological Characterization
- DESeq2 (High vs Low risk)
- Adjusted for age/stage/grade
- Hallmark GSEA
- Multi-panel summary figure

---

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

---

## How to Run the Full Pipeline

### Option 1 — From R

From project root:

```r
source("scripts/run_all.R")
Option 2 — From command line
Rscript scripts/run_all.R
This will automatically generate:

* All intermediate objects

* All tables

* All figures

* Calibration analysis

* Final multi-panel summary figure

* Session information

## Key Outputs
## Survival Modeling

* Cox_univariate_Hallmarks.csv

* SelectedPathways_LASSO_beta.csv

* SurvivalTable_molecularRisk.csv

* Model_compare_AIC_Cindex.csv

* Model_summaries.txt

## Validation

* VALID_rms_validate.txt

* VALID_cox_zph.txt

* VALID_lasso_selection_stability.csv

* timeROC_AUC.csv

## Biology

DESeq2_mgroup_adjusted_fullResults.csv

DESeq2_mgroup_significant_genes.csv

GSEA_Hallmark_mgroup.csv

Figures

PCA

Volcano

Hallmark GSEA

Kaplan–Meier

Time-dependent ROC

Calibration plot

FIGURE_SummaryPanel_Molecular.png

Reproducibility

Fixed random seeds

Bootstrap validation implemented

All intermediate objects cached

sessionInfo saved

Deterministic pipeline

MIT License

Author

Somayeh Sarirchi
Computational Oncology / Bioinformatics
