# =====================================================
# 07_extract_readme_numbers.R
# Extract key numbers for README from pipeline outputs
# Output:
#   results/README_numbers.txt
# =====================================================

source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(readr)
})

assert_dir(DIR$res)

out_txt <- file.path(DIR$res, "README_numbers.txt")

# -----------------------------
# Helpers
# -----------------------------
fmt <- function(x, digits=2) {
  if (length(x)==0 || all(is.na(x))) return(NA_character_)
  format(round(as.numeric(x), digits), nsmall = digits, trim = TRUE)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

extract_cox_table <- function(fit) {
  s <- summary(fit)
  co <- as.data.frame(s$coefficients)
  ci <- as.data.frame(s$conf.int)

  out <- data.frame(
    term = rownames(co),
    HR   = ci$`exp(coef)`,
    CI_low = ci$`lower .95`,
    CI_high= ci$`upper .95`,
    p = co$`Pr(>|z|)`,
    row.names = NULL,
    check.names = FALSE
  )
  out
}

# -----------------------------
# 1) Cohort stats
# -----------------------------
meta_path <- file.path(DIR$cache, "meta2_survival.rds")
if (!file.exists(meta_path)) stop("[07] meta2_survival.rds not found.")
meta2 <- readRDS(meta_path)

meta2$time  <- as.numeric(meta2$time)
meta2$event <- as.numeric(meta2$event)

# use the same survival filters as modeling
meta_s <- meta2 %>% filter(is.finite(time), is.finite(event), time > 0)

N_total <- nrow(meta_s)
N_events <- sum(meta_s$event == 1, na.rm = TRUE)
median_follow <- median(meta_s$time, na.rm = TRUE)

# -----------------------------
# 2) Cox models summary
# -----------------------------
cox_obj_path <- file.path(DIR$cache, "cox_final_models.rds")
cox_old_path <- file.path(DIR$cache, "cox_final.rds")

cox_obj <- NULL
if (file.exists(cox_obj_path)) {
  cox_obj <- readRDS(cox_obj_path)
} else if (file.exists(cox_old_path)) {
  cox_obj <- list(models=list(cox_final=readRDS(cox_old_path)), cox_final=readRDS(cox_old_path))
} else {
  stop("[07] No Cox model object found (cox_final_models.rds / cox_final.rds).")
}

models <- cox_obj$models
cox_final <- cox_obj$cox_final

cox_final_tab <- extract_cox_table(cox_final)

# try to also extract A/B if exist
coxA_tab <- NULL; coxB_tab <- NULL
if (!is.null(models$cox_A_bin)) coxA_tab <- extract_cox_table(models$cox_A_bin)
if (!is.null(models$cox_B_4level)) coxB_tab <- extract_cox_table(models$cox_B_4level)

# -----------------------------
# 3) Model compare table (AIC + C-index)
# -----------------------------
cmp_path <- file.path(DIR$tab, "Model_compare_AIC_Cindex.csv")
cmp <- safe_read_csv(cmp_path)

# -----------------------------
# 4) timeROC AUC
# -----------------------------
auc_path <- file.path(DIR$tab, "timeROC_AUC.csv")
auc <- safe_read_csv(auc_path)

# -----------------------------
# 5) DEG counts
# -----------------------------
deg_path <- file.path(DIR$tab, "DESeq2_mgroup_significant_genes.csv")
deg <- safe_read_csv(deg_path)
deg_n <- if (is.null(deg)) NA_integer_ else nrow(deg)

# -----------------------------
# 6) LASSO stability (top selected)
# -----------------------------
stab_path <- file.path(DIR$tab, "VALID_lasso_selection_stability.csv")
stab <- safe_read_csv(stab_path)

top_stable <- NULL
if (!is.null(stab) && all(c("pathway","freq","prop") %in% colnames(stab))) {
  top_stable <- stab %>%
    arrange(desc(freq), desc(prop)) %>%
    slice_head(n = 5)
}

# -----------------------------
# Write report
# -----------------------------
sink(out_txt)
cat("===== README numbers (auto-extracted) =====\n")
cat("Generated at: ", format(Sys.time()), "\n\n")

cat("## Cohort\n")
cat("Final analyzable cohort (time>0, non-missing): ", N_total, "\n", sep="")
cat("Events (event==1): ", N_events, "\n", sep="")
cat("Median follow-up (days): ", round(median_follow, 1), "\n", sep="")
cat("\n")

cat("## Cox final model (as saved in cox_final_models.rds)\n")
print(cox_final_tab)
cat("\n")

if (!is.null(coxA_tab)) {
  cat("## Model A (cox_A_bin)\n")
  print(coxA_tab)
  cat("\n")
}
if (!is.null(coxB_tab)) {
  cat("## Model B (cox_B_4level)\n")
  print(coxB_tab)
  cat("\n")
}

cat("## Model comparison (AIC + C-index)\n")
if (is.null(cmp)) cat("Model_compare_AIC_Cindex.csv NOT FOUND\n\n") else print(cmp)
cat("\n")

cat("## timeROC AUC\n")
if (is.null(auc)) cat("timeROC_AUC.csv NOT FOUND (maybe timeROC skipped)\n\n") else print(auc)
cat("\n")

cat("## Differential expression\n")
cat("Significant DEGs (padj<0.05 & |log2FC|>=1): ", deg_n, "\n", sep="")
cat("\n")

cat("## LASSO stability (top 5)\n")
if (is.null(top_stable)) cat("VALID_lasso_selection_stability.csv NOT FOUND\n") else print(top_stable)
cat("\n")

cat("===== End =====\n")
sink()

message("[07] DONE: README numbers saved -> ", out_txt)
