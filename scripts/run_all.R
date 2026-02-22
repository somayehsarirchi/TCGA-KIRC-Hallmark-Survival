# =====================================================
# run_all.R
# Run the full pipeline end-to-end (robust + logging)
# =====================================================

source("scripts/00_config.R")

# ---- Logging ----
log_file <- file.path(DIR$res, "pipeline_log.txt")
dir.create(DIR$res, recursive = TRUE, showWarnings = FALSE)

sink(log_file, split = TRUE)
cat("===== PIPELINE START =====\n")
cat("Start time:", format(Sys.time()), "\n\n")

run_step <- function(script_path) {
  cat("--------------------------------------------------\n")
  cat("Running:", script_path, "\n")
  cat("Time:", format(Sys.time()), "\n")
  cat("--------------------------------------------------\n")

  if (!file.exists(script_path)) {
    stop("Script not found: ", script_path)
  }

  tryCatch(
    {
      source(script_path, local = new.env())
      cat("✅ SUCCESS:", script_path, "\n\n")
    },
    error = function(e) {
      cat("❌ FAILED:", script_path, "\n")
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
run_step("scripts/06_summary_panel.R")

cat("===== PIPELINE END =====\n")
cat("End time:", format(Sys.time()), "\n\n")

# Extra: sessionInfo at end (even if step05 is skipped in future edits)
cat("===== sessionInfo() =====\n")
print(sessionInfo())

sink()

message("ALL DONE ✅  (log saved to: ", log_file, ")")
