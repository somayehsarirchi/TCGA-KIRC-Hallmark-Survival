# =====================================================
# 03_survival_modeling.R
# Cox screening + LASSO, risk score, KM, adjusted Cox
# Outputs:
#   results/tables/Cox_univariate_Hallmarks.csv
#   results/tables/SelectedPathways_LASSO_beta.csv
#   results/tables/SurvivalTable_molecularRisk.csv
#   results/tables/Model_compare_AIC_Cindex.csv
#   results/figures/KM_mgroup.png
#   data_cache/cox_final_models.rds
# =====================================================

source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
  library(glmnet)
})

set.seed(1)

pal_group <- c("Low"="#2C7BB6", "High"="#D7191C")

meta2   <- readRDS(file.path(DIR$cache, "meta2_survival.rds"))
hall_df <- readRDS(file.path(DIR$cache, "hall_df_scores.rds"))

# ------------------------------
# 0) Safety checks + harmonize
# ------------------------------

# Ensure rownames exist
if (is.null(rownames(meta2)) || any(rownames(meta2) == "")) {
  stop("[03] meta2 has no rownames. Please set rownames to sample IDs before saving meta2_survival.rds")
}
if (is.null(rownames(hall_df)) || any(rownames(hall_df) == "")) {
  stop("[03] hall_df has no rownames. Please set rownames to sample IDs before saving hall_df_scores.rds")
}

# Align samples by common IDs
common_ids <- intersect(rownames(meta2), rownames(hall_df))
if (length(common_ids) < 50) {
  stop(sprintf("[03] Too few common samples between meta2 and hall_df: %d", length(common_ids)))
}
meta2   <- meta2[common_ids, , drop=FALSE]
hall_df <- hall_df[common_ids, , drop=FALSE]

# Basic survival sanity
meta2 <- meta2 %>%
  mutate(
    time  = as.numeric(time),
    event = as.numeric(event)
  )

meta2 <- meta2 %>%
  filter(!is.na(time), !is.na(event), time > 0)

if (sum(meta2$event == 1, na.rm = TRUE) < 5) {
  stop(sprintf("[03] Non-positive/too few events after filtering. Events=1 count: %d",
               sum(meta2$event == 1, na.rm = TRUE)))
}

# After filtering meta2, keep hall_df aligned
hall_df <- hall_df[rownames(meta2), , drop=FALSE]

# ------------------------------
# 1) Stage + Grade engineering
# ------------------------------

# stage_bin: enforce factor levels if present
if ("stage_bin" %in% colnames(meta2)) {
  meta2$stage_bin <- factor(meta2$stage_bin, levels = c("Early", "Late"))
}

# stage_4level optional (only if you already have it)
# If your meta2 has a 4-level stage column (e.g., stage, ajcc_stage, stage4),
# you can map it here. Otherwise it will be ignored safely.
stage4_candidates <- c("stage_4", "stage4", "ajcc_stage", "stage")
stage4_name <- stage4_candidates[stage4_candidates %in% colnames(meta2)][1]
if (!is.na(stage4_name) && length(stage4_name) == 1) {
  meta2$stage_4 <- as.character(meta2[[stage4_name]])
  meta2$stage_4 <- factor(meta2$stage_4)  # keep as-is unless you want explicit ordering
} else {
  meta2$stage_4 <- NULL
}

# tumor grade: try to find grade column
grade_candidates <- c("tumor_grade", "grade", "histological_grade", "tumorgrade")
grade_name <- grade_candidates[grade_candidates %in% colnames(meta2)][1]

if (!is.na(grade_name) && length(grade_name) == 1) {

  g_raw <- meta2[[grade_name]]

  # Standardize to character, trim
  g_chr <- trimws(as.character(g_raw))

  # Normalize common encodings to G1..G4
  # Examples: "1","2","3","4"  OR "G1","G2"... OR "Grade 1"...
  g_chr <- gsub("^Grade\\s*", "", g_chr, ignore.case = TRUE)
  g_chr <- gsub("^G\\s*", "", g_chr, ignore.case = TRUE)

  # Keep only 1-4
  g_chr[!g_chr %in% c("1","2","3","4")] <- NA

  meta2$tumor_grade4 <- factor(
    paste0("G", g_chr),
    levels = c("G1","G2","G3","G4")
  )

  # Bin version (G1-2 vs G3-4)
  meta2$tumor_grade_bin <- dplyr::case_when(
    meta2$tumor_grade4 %in% c("G1","G2") ~ "Low",
    meta2$tumor_grade4 %in% c("G3","G4") ~ "High",
    TRUE ~ NA_character_
  )
  meta2$tumor_grade_bin <- factor(meta2$tumor_grade_bin, levels = c("Low","High"))

} else {
  meta2$tumor_grade4 <- NULL
  meta2$tumor_grade_bin <- NULL
  message("[03] WARNING: No tumor grade column found in meta2. Final models will run without grade.")
}

# ------------------------------
# 2) Univariate Cox (all hallmarks)
# ------------------------------
# Use safe fitting: skip pathways with zero variance / too many NA
safe_uni_cox <- function(pw_vec, time, event) {
  # remove NA pairs
  ok <- is.finite(pw_vec) & is.finite(time) & is.finite(event)
  pw_vec <- pw_vec[ok]; time <- time[ok]; event <- event[ok]
  if (length(pw_vec) < 30) return(NULL)
  if (sd(pw_vec) == 0) return(NULL)
  fit <- coxph(Surv(time, event) ~ pw_vec)
  s <- summary(fit)
  data.frame(
    HR = s$coef[1, "exp(coef)"],
    p  = s$coef[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

uni_list <- lapply(colnames(hall_df), function(pw){
  res <- safe_uni_cox(hall_df[[pw]], meta2$time, meta2$event)
  if (is.null(res)) return(NULL)
  cbind(data.frame(pathway=pw, stringsAsFactors = FALSE), res)
})
uni_res <- bind_rows(uni_list)

if (nrow(uni_res) == 0) stop("[03] Univariate Cox produced 0 results (check hall_df content / filtering).")

uni_res$padj <- p.adjust(uni_res$p, method="BH")
uni_res <- uni_res %>% arrange(padj)

write.csv(uni_res, file.path(DIR$tab, "Cox_univariate_Hallmarks.csv"), row.names = FALSE)

topN <- 5
top_pw <- head(uni_res$pathway, topN)
message("[03] Top hallmarks:")
print(top_pw)

# ------------------------------
# 3) LASSO Cox on top hallmarks
# ------------------------------
X_pw   <- as.matrix(hall_df[, top_pw, drop=FALSE])
y_surv <- Surv(meta2$time, meta2$event)

# If any NA in X, glmnet will fail
if (any(!is.finite(X_pw))) {
  stop("[03] NA/Inf detected in X_pw. Please check hallmark scores (hall_df).")
}

cv_pw  <- cv.glmnet(X_pw, y_surv, family="cox", alpha=1, nfolds=10)

# Try lambda.min; if selects 0, fallback to lambda.1se; if still 0, stop with guidance
fit_min <- glmnet(X_pw, y_surv, family="cox", alpha=1, lambda=cv_pw$lambda.min)
coef_min <- as.matrix(coef(fit_min))
sel_min  <- rownames(coef_min)[coef_min[,1] != 0]

chosen_lambda <- "lambda.min"
fit_pw <- fit_min
coef_pw <- coef_min
sel_pw <- sel_min

if (length(sel_pw) == 0) {
  fit_1se <- glmnet(X_pw, y_surv, family="cox", alpha=1, lambda=cv_pw$lambda.1se)
  coef_1se <- as.matrix(coef(fit_1se))
  sel_1se  <- rownames(coef_1se)[coef_1se[,1] != 0]
  if (length(sel_1se) == 0) {
    stop("[03] LASSO selected 0 pathways with both lambda.min and lambda.1se. Expand features (top10/20) or use ridge/elastic-net.")
  } else {
    chosen_lambda <- "lambda.1se"
    fit_pw <- fit_1se
    coef_pw <- coef_1se
    sel_pw <- sel_1se
  }
}

beta_pw <- coef_pw[sel_pw, 1, drop=TRUE]

write.csv(
  data.frame(pathway=sel_pw, beta=as.numeric(beta_pw)),
  file.path(DIR$tab, "SelectedPathways_LASSO_beta.csv"),
  row.names=FALSE
)

message(sprintf("[03] LASSO chosen: %s | selected: %d", chosen_lambda, length(sel_pw)))

# Risk score
meta2$pw_risk <- as.numeric(X_pw[, sel_pw, drop=FALSE] %*% beta_pw)
meta2$mgroup  <- factor(ifelse(meta2$pw_risk >= median(meta2$pw_risk, na.rm=TRUE), "High", "Low"),
                        levels=c("Low","High"))

write.csv(
  meta2[, intersect(c("time","event","age_years","stage_bin","stage_4","tumor_grade4","tumor_grade_bin","pw_risk","mgroup"),
                    colnames(meta2)), drop=FALSE],
  file.path(DIR$tab, "SurvivalTable_molecularRisk.csv"),
  row.names=TRUE
)

# ------------------------------
# 4) KM plot
# ------------------------------
fit_km <- survfit(Surv(time, event) ~ mgroup, data=meta2)

p_km <- ggsurvplot(
  fit_km, data=meta2,
  pval=TRUE, risk.table=TRUE,
  palette=c(pal_group["Low"], pal_group["High"]),
  legend.labs=c("Low","High"),
  legend.title="Risk group",
  title="Molecular group (Hallmark-pathway risk): Overall survival",
  xlab="Time (days)",
  ylab="Overall survival probability",
  risk.table.height = 0.25
)

ggsave(file.path(DIR$fig, "KM_mgroup.png"), p_km$plot, width=7, height=6, dpi=300)

# ------------------------------
# 5) Adjusted Cox models (with grade)
#    - Model A: stage_bin + grade_bin
#    - Model B: stage_4 + grade4  (if available)
# ------------------------------

# helper: C-index from survConcordance (base survival)
get_cindex <- function(fit) {
  sc <- survConcordance(fit$y ~ predict(fit))
  as.numeric(sc$concordance)
}

models <- list()

# Model A (bin/bin)
vars_A <- c("pw_risk", "age_years")
if ("stage_bin" %in% colnames(meta2)) vars_A <- c(vars_A, "stage_bin")
if ("tumor_grade_bin" %in% colnames(meta2)) vars_A <- c(vars_A, "tumor_grade_bin")

dfA <- meta2 %>% dplyr::select(all_of(c("time","event", vars_A))) %>% na.omit()
if (nrow(dfA) >= 50 && sum(dfA$event==1) >= 5) {
  fA <- as.formula(paste0("Surv(time, event) ~ ", paste(vars_A, collapse=" + ")))
  models$cox_A_bin <- coxph(fA, data=dfA, x=TRUE, y=TRUE)
} else {
  message("[03] WARNING: Model A could not run (too few rows/events after NA removal).")
}

# Model B (4-level, optional)
vars_B <- c("pw_risk", "age_years")
if ("stage_4" %in% colnames(meta2)) vars_B <- c(vars_B, "stage_4")
if ("tumor_grade4" %in% colnames(meta2)) vars_B <- c(vars_B, "tumor_grade4")

dfB <- meta2 %>% dplyr::select(all_of(c("time","event", vars_B))) %>% na.omit()
if (all(c("stage_4","tumor_grade4") %in% colnames(meta2)) &&
    nrow(dfB) >= 50 && sum(dfB$event==1) >= 5) {
  fB <- as.formula(paste0("Surv(time, event) ~ ", paste(vars_B, collapse=" + ")))
  models$cox_B_4level <- coxph(fB, data=dfB, x=TRUE, y=TRUE)
} else {
  message("[03] INFO: Model B (4-level) not available (missing stage_4 or grade4, or too few complete rows).")
}

# Keep a "final" pointer: prefer Model A if exists else any model
cox_final <- if (!is.null(models$cox_A_bin)) models$cox_A_bin else models[[1]]

# Save models
saveRDS(
  list(
    models=models,
    cox_final=cox_final,
    lasso=list(
      top_features=top_pw,
      selected=sel_pw,
      beta=beta_pw,
      chosen_lambda=chosen_lambda,
      lambda.min=cv_pw$lambda.min,
      lambda.1se=cv_pw$lambda.1se
    )
  ),
  file.path(DIR$cache, "cox_final_models.rds")
)

# Compare AIC + C-index
cmp <- data.frame(
  model=character(),
  n=integer(),
  events=integer(),
  AIC=numeric(),
  C_index=numeric(),
  stringsAsFactors=FALSE
)

if (!is.null(models$cox_A_bin)) {
  cmp <- rbind(cmp, data.frame(
    model="cox_A_bin (stage_bin + grade_bin)",
    n=nrow(dfA),
    events=sum(dfA$event==1),
    AIC=as.numeric(AIC(models$cox_A_bin)),
    C_index=get_cindex(models$cox_A_bin)
  ))
}
if (!is.null(models$cox_B_4level)) {
  cmp <- rbind(cmp, data.frame(
    model="cox_B_4level (stage_4 + grade4)",
    n=nrow(dfB),
    events=sum(dfB$event==1),
    AIC=as.numeric(AIC(models$cox_B_4level)),
    C_index=get_cindex(models$cox_B_4level)
  ))
}

write.csv(cmp, file.path(DIR$tab, "Model_compare_AIC_Cindex.csv"), row.names=FALSE)

# Write summaries
sink(file.path(DIR$res, "Model_summaries.txt"))
cat("===== Univariate hallmarks: saved to Cox_univariate_Hallmarks.csv =====\n\n")

cat("===== LASSO Cox =====\n")
cat("TopN candidates:", paste(top_pw, collapse=", "), "\n")
cat("Chosen lambda:", chosen_lambda, "\n")
cat("lambda.min:", cv_pw$lambda.min, "\n")
cat("lambda.1se:", cv_pw$lambda.1se, "\n\n")
cat("Selected pathways:\n")
print(data.frame(pathway=sel_pw, beta=as.numeric(beta_pw)))

cat("\n===== KM: mgroup =====\n")
cat("KM plot saved: results/figures/KM_mgroup.png\n\n")

cat("===== Adjusted Cox models =====\n\n")
if (!is.null(models$cox_A_bin)) {
  cat("---- Model A (bin/bin) ----\n")
  print(summary(models$cox_A_bin))
  cat("\n")
}
if (!is.null(models$cox_B_4level)) {
  cat("---- Model B (4-level) ----\n")
  print(summary(models$cox_B_4level))
  cat("\n")
}

cat("===== Model comparison (AIC + C-index) =====\n")
print(cmp)
sink()

message("[03] DONE: Cox + LASSO + KM + adjusted models (with grade if available) saved.")
