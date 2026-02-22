# =====================================================
# 05_validation.R
# timeROC + PH test + LASSO stability + sessionInfo
# Outputs:
#   results/figures/PANEL_timeROC.png
#   results/VALID_cox_zph.txt + VALID_cox_zph.png
#   results/tables/VALID_lasso_selection_stability.csv
#   results/tables/timeROC_AUC.csv
#   results/sessionInfo.txt
#   data_cache/timeROC_object.rds
# =====================================================

source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(timeROC)
  library(dplyr)
})

set.seed(1)

meta2   <- readRDS(file.path(DIR$cache, "meta2_survival.rds"))
hall_df <- readRDS(file.path(DIR$cache, "hall_df_scores.rds"))

# ------------------------------
# 0) Load final Cox model robustly
# ------------------------------
cox_final <- NULL

new_model_path <- file.path(DIR$cache, "cox_final_models.rds")
old_model_path <- file.path(DIR$cache, "cox_final.rds")

if (file.exists(new_model_path)) {
  obj <- readRDS(new_model_path)
  if (!is.null(obj$cox_final)) cox_final <- obj$cox_final
} else if (file.exists(old_model_path)) {
  cox_final <- readRDS(old_model_path)
} else {
  stop("[05] No Cox model found. Run step03 first (cox_final_models.rds or cox_final.rds).")
}

if (is.null(cox_final) || !inherits(cox_final, "coxph")) {
  stop("[05] Loaded cox_final is not a coxph object. Check step03 outputs.")
}

# ------------------------------
# 1) Bring pw_risk from survival table & harmonize samples
# ------------------------------
survtab_path <- file.path(DIR$tab, "SurvivalTable_molecularRisk.csv")
if (!file.exists(survtab_path)) stop("[05] SurvivalTable_molecularRisk.csv not found. Run step03 first.")

survtab <- read.csv(survtab_path, row.names=1, check.names=FALSE)

# Ensure rownames exist
if (is.null(rownames(meta2))) stop("[05] meta2 has no rownames (sample IDs).")
if (is.null(rownames(hall_df))) stop("[05] hall_df has no rownames (sample IDs).")

common_ids <- Reduce(intersect, list(rownames(meta2), rownames(hall_df), rownames(survtab)))
if (length(common_ids) < 50) {
  stop(sprintf("[05] Too few common samples across meta2/hall_df/survtab: %d", length(common_ids)))
}

meta2   <- meta2[common_ids, , drop=FALSE]
hall_df <- hall_df[common_ids, , drop=FALSE]
survtab <- survtab[common_ids, , drop=FALSE]

meta2$pw_risk <- as.numeric(survtab[, "pw_risk"])

# Survival sanity
meta2$time  <- as.numeric(meta2$time)
meta2$event <- as.numeric(meta2$event)

keep <- is.finite(meta2$time) & is.finite(meta2$event) & meta2$time > 0 & is.finite(meta2$pw_risk)
meta2   <- meta2[keep, , drop=FALSE]
hall_df <- hall_df[rownames(meta2), , drop=FALSE]

if (sum(meta2$event == 1, na.rm=TRUE) < 5) {
  stop(sprintf("[05] Too few events after filtering: %d", sum(meta2$event == 1, na.rm=TRUE)))
}

# ------------------------------
# 2) timeROC
# ------------------------------
# Desired times (days)
times0 <- c(365, 1095, 1825)

# Ensure requested times are within follow-up range (timeROC can misbehave otherwise)
tmax <- max(meta2$time, na.rm=TRUE)
times <- times0[times0 < tmax]

if (length(times) == 0) {
  message(sprintf("[05] WARNING: All requested times exceed max follow-up (tmax=%.1f). Skipping timeROC.", tmax))
} else {

  roc_pw <- timeROC(
    T     = meta2$time,
    delta = meta2$event,
    marker= meta2$pw_risk,
    cause = 1,
    times = times,
    iid   = TRUE
  )

  write.csv(
    data.frame(time=times, AUC=roc_pw$AUC),
    file.path(DIR$tab, "timeROC_AUC.csv"),
    row.names=FALSE
  )

  png(file.path(DIR$fig, "PANEL_timeROC.png"), width=900, height=750, res=150)

  # Plot each time with different line types (avoid relying on color only)
  ltys <- c(1,2,3,4,5)
  ltys <- ltys[seq_along(times)]

  plot(roc_pw, time=times[1], lwd=2, lty=ltys[1])
  if (length(times) >= 2) {
    for (i in 2:length(times)) plot(roc_pw, time=times[i], lwd=2, lty=ltys[i], add=TRUE)
  }

  leg <- paste0(round(times/365, 1), "y AUC=", sprintf("%.2f", roc_pw$AUC))
  legend("bottomright", legend=leg, lty=ltys, lwd=2, bty="n")

  dev.off()

  saveRDS(roc_pw, file.path(DIR$cache, "timeROC_object.rds"))
  message("[05] timeROC saved for times: ", paste(times, collapse=", "))
}

# ------------------------------
# 3) PH assumption test (cox.zph)
# ------------------------------
zph_out_txt <- file.path(DIR$res, "VALID_cox_zph.txt")
zph_out_png <- file.path(DIR$res, "VALID_cox_zph.png")

zph <- NULL
zph_ok <- TRUE
tryCatch({
  zph <- cox.zph(cox_final)
}, error=function(e){
  zph_ok <<- FALSE
  capture.output(
    cat("[05] cox.zph failed:\n", conditionMessage(e), "\n"),
    file = zph_out_txt
  )
})

if (zph_ok) {
  capture.output(print(zph), file=zph_out_txt)

  png(zph_out_png, width=1100, height=800, res=150)
  plot(zph)
  dev.off()

  message("[05] PH test saved: VALID_cox_zph.txt/png")
}

# ------------------------------
# 4) LASSO stability (bootstrap selection frequency)
# ------------------------------
uni_path <- file.path(DIR$tab, "Cox_univariate_Hallmarks.csv")
if (!file.exists(uni_path)) stop("[05] Cox_univariate_Hallmarks.csv not found. Run step03 first.")
uni_res <- read.csv(uni_path, stringsAsFactors=FALSE)

topN <- 5
top_pw <- head(uni_res$pathway, topN)

# safety: keep only those that exist in hall_df
top_pw <- top_pw[top_pw %in% colnames(hall_df)]
if (length(top_pw) < 2) stop("[05] Too few top pathways found in hall_df to run LASSO stability.")

X_pw   <- as.matrix(hall_df[, top_pw, drop=FALSE])
y_surv <- Surv(meta2$time, meta2$event)

if (any(!is.finite(X_pw))) stop("[05] NA/Inf detected in hallmark matrix for top pathways.")

# CV once on full data to choose lambda
cv_pw <- cv.glmnet(X_pw, y_surv, family="cox", alpha=1, nfolds=10)

B <- 50
sel_list <- vector("list", B)
kept <- 0

for (b in seq_len(B)) {
  idx <- sample(seq_len(nrow(meta2)), replace=TRUE)

  # need enough events in bootstrap sample
  if (sum(meta2$event[idx] == 1, na.rm=TRUE) < 5) next

  Xb <- X_pw[idx, , drop=FALSE]
  yb <- Surv(meta2$time[idx], meta2$event[idx])

  fit <- glmnet(Xb, yb, family="cox", alpha=1, lambda=cv_pw$lambda.min)
  co  <- as.matrix(coef(fit))
  sel <- rownames(co)[co[,1] != 0]

  kept <- kept + 1
  sel_list[[kept]] <- sel
}

sel_list <- sel_list[seq_len(kept)]
if (kept < 10) {
  message(sprintf("[05] WARNING: Only %d/%d bootstraps usable (low events). Stability may be noisy.", kept, B))
}

tab <- sort(table(unlist(sel_list)), decreasing=TRUE)

stab_df <- data.frame(
  pathway = names(tab),
  freq    = as.integer(tab),
  prop    = as.numeric(tab) / max(1, kept),
  stringsAsFactors=FALSE
)

# ensure all top features appear even if never selected
miss_pw <- setdiff(top_pw, stab_df$pathway)
if (length(miss_pw) > 0) {
  stab_df <- rbind(stab_df, data.frame(pathway=miss_pw, freq=0L, prop=0, stringsAsFactors=FALSE))
  stab_df <- stab_df %>% arrange(desc(freq), desc(prop), pathway)
}

write.csv(
  stab_df,
  file.path(DIR$tab, "VALID_lasso_selection_stability.csv"),
  row.names=FALSE
)

# ------------------------------
# 5) sessionInfo
# ------------------------------
sink(file.path(DIR$res, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("[05] DONE: validation outputs saved.")
