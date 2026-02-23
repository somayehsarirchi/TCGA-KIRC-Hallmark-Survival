# =====================================================
# 07_rms_calibration.R
# rms validate + calibrate for final adjusted Cox model
# Outputs:
#   results/VALID_rms_validate.txt
#   results/figures/PANEL_calibration.png
#   data_cache/rms_cph_final.rds
# =====================================================

source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(survival)
  library(rms)
})

set.seed(1)

assert_dir(DIR$res); assert_dir(DIR$fig); assert_dir(DIR$cache)

# ------------------------------
# 0) Load final Cox model (from step03)
# ------------------------------
new_model_path <- file.path(DIR$cache, "cox_final_models.rds")
old_model_path <- file.path(DIR$cache, "cox_final.rds")

cox_final <- NULL
if (file.exists(new_model_path)) {
  obj <- readRDS(new_model_path)
  if (!is.null(obj$cox_final)) cox_final <- obj$cox_final
} else if (file.exists(old_model_path)) {
  cox_final <- readRDS(old_model_path)
} else {
  stop("[07] No Cox model found. Run step03 first.")
}
if (is.null(cox_final) || !inherits(cox_final, "coxph")) {
  stop("[07] cox_final is not a coxph object.")
}

# Extract the exact data used by coxph (after na.omit inside step03)
mf <- model.frame(cox_final)
if (!all(c("time","event") %in% colnames(mf))) {
  # model.frame(cph) usually contains Surv() as response, not time/event columns
  # We rebuild time/event by parsing response
  y <- model.response(mf)
  if (!inherits(y, "Surv")) stop("[07] Could not recover Surv response from model.frame(cox_final).")
  mf$time  <- as.numeric(y[, "time"])
  mf$event <- as.numeric(y[, "status"])
}

# Reconstruct formula as Surv(time,event) ~ predictors
rhs <- attr(terms(cox_final), "term.labels")
if (length(rhs) == 0) stop("[07] No predictors found in cox_final terms.")
f_rms <- as.formula(paste0("Surv(time, event) ~ ", paste(rhs, collapse=" + ")))

# ------------------------------
# 1) Fit rms::cph with same formula and data
# ------------------------------
dd <- datadist(mf)
options(datadist = "dd")

cph_fit <- cph(
  formula = f_rms,
  data    = mf,
  x = TRUE, y = TRUE,
  surv = TRUE, time.inc = 1095
)

saveRDS(cph_fit, file.path(DIR$cache, "rms_cph_final.rds"))

# ------------------------------
# 2) Validate (bootstrap)
# ------------------------------
B <- 500
val <- validate(cph_fit, method = "boot", B = B)

# Save validate output
val_txt <- file.path(DIR$res, "VALID_rms_validate.txt")
capture.output(
  cat("===== rms::validate (boot) =====\n"),
  cat("B =", B, "\n"),
  cat("Model:", deparse(f_rms), "\n\n"),
  print(val),
  file = val_txt
)

# ------------------------------
# 3) Calibrate (bootstrap) at 3 years (1095 days)
# ------------------------------
u_time <- 1095

cal <- calibrate(
  cph_fit,
  method = "boot",
  B      = 200,     # شروع سبک‌تر برای کالیبراسیون (می‌توانی 500 هم کنی)
  u      = u_time
)

# Plot calibration (clean, balanced font)
out_png <- file.path(DIR$fig, "PANEL_calibration.png")
png(out_png, width = 900, height = 750, res = 150)

par(
  mar = c(4.2, 4.2, 2.2, 1.2),
  cex.axis = 0.95,
  cex.lab  = 1.0,
  cex.main = 1.05
)

plot(
  cal,
  xlab = "Predicted survival probability",
  ylab = "Observed survival probability",
  subtitles = FALSE,
  lwd = 2,
  col = 1
)

abline(0, 1, lty = 2)
title(main = paste0("Calibration at ", round(u_time/365,1), " years (bootstrap)"))

legend(
  "topleft",
  legend = c("Bootstrap-corrected", "Ideal (45°)"),
  lty = c(1, 2),
  lwd = c(2, 1),
  bty = "n"
)

dev.off()

message("[07] DONE: rms validate + calibration saved.")
message("[07] Validate text: ", val_txt)
message("[07] Calibration plot: ", out_png)
