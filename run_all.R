suppressPackageStartupMessages({
  library(here)
})

source(here("R","00_utils.R"))
source(here("R","01_gdc.R"))
source(here("R","02_preprocess.R"))
source(here("R","03_vst_symbol.R"))
source(here("R","04_ssgsea.R"))
source(here("R","05_uni_cox.R"))
source(here("R","06_lasso_risk.R"))
source(here("R","07_deseq2_gsea.R"))
source(here("R","08_eval_panels.R"))
source(here("R","09_train_test.R"))

main <- function(force=FALSE,
                 use_pec_calibration=TRUE,
                 use_rms_calibration=TRUE,
                 run_train_test=TRUE) {
  
  dir_create(here("logs"))
  dir_create(here("results","figures"))
  dir_create(here("results","tables"))
  dir_create(here("results","objects"))
  
  log_file <- here("logs", paste0("pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
  log_msg("Project root: ", here::here(), log_file=log_file)
  log_msg("force=", force,
          " use_pec_calibration=", use_pec_calibration,
          " use_rms_calibration=", use_rms_calibration,
          " run_train_test=", run_train_test, log_file=log_file)
  
  tryCatch({
    step_01_gdc(log_file=log_file, force=force)
    step_02_preprocess(log_file=log_file, force=force)
    step_03_vst_symbol(log_file=log_file, force=force)
    step_04_ssgsea(log_file=log_file, force=force)
    step_05_uni_cox(log_file=log_file, force=force)
    step_06_lasso_risk(log_file=log_file, force=force)
    step_07_deseq2_gsea(log_file=log_file, force=force)
    step_08_eval_panels(log_file=log_file, force=force,
                        use_pec_calibration=use_pec_calibration,
                        use_rms_calibration=use_rms_calibration)
    
    if (isTRUE(run_train_test)) {
      step_09_train_test(log_file=log_file, force=force)
    }
    
    log_msg("✅ Pipeline finished successfully.", log_file=log_file)
  }, error=function(e) {
    log_msg("❌ Pipeline failed: ", conditionMessage(e), level="ERROR", log_file=log_file)
    stop(e)
  })
}

main(force=FALSE, use_pec_calibration=TRUE, use_rms_calibration=TRUE, run_train_test=TRUE)