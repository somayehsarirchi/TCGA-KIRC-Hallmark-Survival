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
