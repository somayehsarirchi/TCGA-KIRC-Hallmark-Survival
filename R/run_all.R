# =====================================================
# run_all.R
# Run the full pipeline end-to-end (robust + logging)
# =====================================================

source("scripts/00_config.R")

# ---- Logging ----
log_file <- file.path(DIR$res, "pipeline_log.txt")
assert_dir(DIR$res)

sink(log_file, split = TRUE)
on.exit({ sink() }, add = TRUE)

cat("===== PIPELINE START =====\n")
cat("Project root:", DIR$root, "\n")
cat("Start time:", format(Sys.time()), "\n\n")

run_step <- function(script_path) {
  cat("--------------------------------------------------\n")
  cat("Running:", script_path, "\n")
  cat("Start:", format(Sys.time()), "\n")
  cat("--------------------------------------------------\n")

  if (!file.exists(script_path)) stop("Script not found: ", script_path)

  t0 <- Sys.time()

  tryCatch(
    {
      # isolate each step in its own environment
      sys.source(script_path, envir = new.env())
      dt <- difftime(Sys.time(), t0, units = "secs")
      cat("✅ SUCCESS:", script_path, "\n")
      cat("Duration (sec):", round(as.numeric(dt), 2), "\n\n")
    },
    error = function(e) {
      dt <- difftime(Sys.time(), t0, units = "secs")
      cat("❌ FAILED:", script_path, "\n")
      cat("Duration (sec):", round(as.numeric(dt), 2), "\n")
      cat("Error message:\n", conditionMessage(e), "\n\n")
      stop(e)
    }
  )
}

# ---- Steps ----
run_step("scripts/01_download_preprocess.R")
run_step("scripts/02_ssgsea_scoring.R")
run_step("scripts/03_survival_modeling.R")
run_step("scripts/04_DEG_GSEA.R")
run_step("scripts/05_validation.R")
run_step("scripts/07_rms_calibration.R")
run_step("scripts/06_summary_panel.R")

cat("===== PIPELINE END =====\n")
cat("End time:", format(Sys.time()), "\n\n")

cat("===== sessionInfo() =====\n")
print(sessionInfo())
cat("\n")

message("ALL DONE ✅  (log saved to: ", log_file, ")")
