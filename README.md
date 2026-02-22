
Molecular Hallmark–Based Survival Modeling in TCGA-KIRC
Overview

This repository implements a fully reproducible computational oncology pipeline for molecular risk stratification in TCGA Kidney Renal Clear Cell Carcinoma (KIRC).

The goal was to:

Quantify Hallmark pathway activity using ssGSEA

Identify prognostic pathways via Cox regression

Build a regularized (LASSO) survival model

Construct a molecular risk score

Validate discrimination, calibration, and model assumptions

Characterize downstream tumor biology using DESeq2 and GSEA

Dataset

TCGA-KIRC (Primary Tumor samples)

n = 437 patients

Events = 125

Median follow-up ≈ 1140 days

RNA-Seq: STAR counts

Clinical covariates: age, AJCC stage

Molecular Risk Model

Adjusted Cox model:

Surv(time, event) ~ pw_risk + age + stage

Results:

Molecular risk HR = 2.61 (p = 2.1e-10)

Age HR = 1.038 (p < 1e-5)

Stage (Late vs Early) HR = 2.72 (p < 1e-7)

Concordance index (C-index) = 0.775

The model demonstrates strong discrimination.

Validation Strategy

✔ Train/Test split (70/30, no data leakage)
✔ Time-dependent ROC (1y / 3y / 5y)
✔ Bootstrap calibration (rms + PEC)
✔ Proportional hazards testing (Schoenfeld residuals)
✔ LASSO feature selection stability (50 bootstrap resamples)

Feature Stability

Consistently selected pathways (100% frequency):

HALLMARK_BILE_ACID_METABOLISM

HALLMARK_UNFOLDED_PROTEIN_RESPONSE

Differential Expression Analysis

High vs Low molecular risk comparison:

Significant DEGs = 1965

Adjusted for age and stage

Downstream Hallmark GSEA performed

Biological Insights

Top prognostic hallmarks:

Unfolded Protein Response

G2M Checkpoint

MYC Targets

E2F Targets

IL6–JAK–STAT3 signaling

These pathways reflect proliferative stress, oncogenic activation, and inflammatory signaling consistent with renal carcinoma biology.

Reproducibility

All random seeds fixed

sessionInfo saved

Cached GDC query objects

Fully scripted pipeline

Bootstrap validation included

Summary Panel

<img width="2400" height="1600" alt="FIGURE_SummaryPanel_Molecular" src="https://github.com/user-attachments/assets/801be9e1-3156-479a-b829-26ca58ec423a" />


Author

Somayeh Sarirchi
Computational Oncology / Bioinformatics
