# TCGA-KIRC Hallmark Survival Modeling

[![R](https://img.shields.io/badge/language-R-blue.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)]()

---

## Overview

This repository presents a fully reproducible computational oncology pipeline for molecular risk stratification in TCGA Kidney Renal Clear Cell Carcinoma (KIRC).

The framework integrates pathway-level activity scoring, penalized survival modeling, internal validation, and biological characterization.

The primary objective is to evaluate whether a Hallmark-derived molecular risk score provides independent prognostic value beyond established clinicopathological variables.

---

## Objectives

- Quantify Hallmark pathway activity using ssGSEA
- Identify prognostic pathways via univariable Cox regression
- Construct a molecular risk score using LASSO-regularized Cox modeling
- Evaluate independent prognostic value adjusting for:
  - Age
  - AJCC stage
  - Tumor grade
- Assess model discrimination and calibration
- Characterize downstream tumor biology using DESeq2 and GSEA

---

## Dataset

| Variable | Value |
|----------|-------|
| Cohort | TCGA-KIRC (Primary Tumor) |
| Initial patients | 437 |
| Final analyzable cohort | 428 |
| Events | 125 |
| Median follow-up | ~1140 days |
| RNA-Seq | STAR counts |
| Clinical covariates | Age, AJCC stage, Tumor grade |

Data were accessed using TCGAbiolinks.

---

## Multivariable Survival Model

Final adjusted Cox model:

Surv(time, event) ~ pw_risk + age_years + stage_bin + grade_bin

### Results

| Variable | Hazard Ratio | 95% CI |
|----------|--------------|--------|
| Molecular risk score | 1.65 | 1.39 – 1.96 |
| Age | 1.93 | 1.43 – 2.60 |
| Stage (Late vs Early) | 2.40 | 1.64 – 3.53 |
| Grade (High vs Low) | 1.50* | borderline |

*Tumor grade showed borderline significance after adjustment for stage.

---

## Model Performance

- Apparent C-index: ~0.77
- Bootstrap optimism-corrected C-index: ~0.76
- Calibration slope: 0.96
- Minimal optimism (Δ ≈ 0.013)

The molecular risk score remains independently associated with survival after full clinical adjustment.

---

## Time-Dependent Discrimination

| Time | AUC |
|------|------|
| 1 year | 0.85 |
| 3 years | 0.78 |
| 5 years | 0.75 |

---

## Internal Validation Strategy

- 500-bootstrap internal validation
- Optimism-corrected discrimination (Dxy → C-index)
- Calibration at 3 years
- Proportional hazards testing
- LASSO feature stability assessment

---

## Feature Stability

Pathways consistently selected across bootstrap resamples:

- HALLMARK_BILE_ACID_METABOLISM
- HALLMARK_UNFOLDED_PROTEIN_RESPONSE

---

## Differential Expression Analysis

High vs Low molecular risk comparison:

- Significant DEGs: 1965
- Adjusted for age and stage
- Downstream Hallmark GSEA performed

---

## Pipeline Structure

data/
scripts/
results/

To run the full analysis:

source("scripts/run_all.R")

---

## Reproducibility

- Random seeds fixed
- Bootstrap validation implemented
- sessionInfo saved
- No data leakage
- MIT License

---

## Author

Somayeh Sarirchi  
Computational Oncology / Bioinformatics
