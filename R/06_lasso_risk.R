source(here::here("R","00_utils.R"))

step_06_lasso_risk <- function(cache_in_meta=here::here("cache","02_preprocess"),
                               cache_in_scores=here::here("cache","04_scores"),
                               cache_in_uni=here::here("cache","05_unicox"),
                               cache_out=here::here("cache","06_lasso"),
                               out_tables=here::here("results","tables"),
                               log_file=NULL,
                               force=FALSE) {
  
  suppressPackageStartupMessages({
    library(glmnet)
    library(survival)
  })
  
  dir_create(cache_out); dir_create(out_tables)
  
  in_meta2 <- file.path(cache_in_meta, "meta2.rds")
  in_halldf <- file.path(cache_in_scores, "hall_df.rds")
  in_top5 <- file.path(cache_in_uni, "top5.rds")
  
  out_cv <- file.path(cache_out, "cv_pw.rds")
  out_fit <- file.path(cache_out, "fit_pw.rds")
  out_beta <- file.path(cache_out, "beta_pw.rds")
  out_sel <- file.path(cache_out, "sel_pw.rds")
  out_meta_risk <- file.path(cache_out, "meta2_risk.rds")
  
  out_surv_csv <- file.path(out_tables, "SurvivalTable_molecularRisk.csv")
  out_beta_csv <- file.path(out_tables, "SelectedPathways_LASSO_beta.csv")
  
  run_step("06_lasso_risk",
           c(out_cv, out_fit, out_beta, out_sel, out_meta_risk, out_surv_csv, out_beta_csv),
           force=force, log_file=log_file,
           fn_build=function() {
             
             meta2 <- read_rds(in_meta2, log_file)
             hall_df <- read_rds(in_halldf, log_file)
             top5 <- read_rds(in_top5, log_file)
             
             X_pw <- as.matrix(hall_df[, top5, drop=FALSE])
             y_surv <- Surv(meta2$time, meta2$event)
             
             cv_pw <- cv.glmnet(X_pw, y_surv, family="cox", alpha=1, nfolds=10)
             fit_pw <- glmnet(X_pw, y_surv, family="cox", alpha=1, lambda=cv_pw$lambda.min)
             
             coef_pw <- as.matrix(coef(fit_pw))
             sel_pw <- rownames(coef_pw)[coef_pw[,1] != 0]
             if (length(sel_pw) == 0) stopf("LASSO selected 0 pathways. Consider lambda.1se or more pathways.", log_file=log_file)
             
             beta_pw <- coef_pw[sel_pw, 1, drop=TRUE]
             
             meta2$pw_risk <- as.numeric(X_pw[, sel_pw, drop=FALSE] %*% beta_pw)
             meta2$mgroup <- factor(ifelse(meta2$pw_risk >= median(meta2$pw_risk, na.rm=TRUE), "High", "Low"),
                                    levels=c("Low","High"))
             
             write.csv(meta2[, c("time","event","age_years","gender","stage_bin","grade_bin","pw_risk","mgroup")],
                       out_surv_csv, row.names=TRUE)
             write.csv(data.frame(pathway=sel_pw, beta=as.numeric(beta_pw)),
                       out_beta_csv, row.names=FALSE)
             
             save_rds(cv_pw, out_cv, log_file)
             save_rds(fit_pw, out_fit, log_file)
             save_rds(beta_pw, out_beta, log_file)
             save_rds(sel_pw, out_sel, log_file)
             save_rds(meta2, out_meta_risk, log_file)
           })
  
  invisible(TRUE)
}