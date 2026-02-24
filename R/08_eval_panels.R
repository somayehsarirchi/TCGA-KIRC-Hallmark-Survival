# (KM/ROC/Calibration/Workflow/SummaryPanel/VALID)
source(here::here("R","00_utils.R"))

step_08_eval_panels <- function(cache_in_lasso=here::here("cache","06_lasso"),
                                cache_in_scores=here::here("cache","04_scores"),
                                cache_in_uni=here::here("cache","05_unicox"),
                                fig_dir=here::here("results","figures"),
                                out_tables=here::here("results","tables"),
                                out_objects=here::here("results","objects"),
                                cache_out=here::here("cache","08_eval"),
                                use_pec_calibration=TRUE,
                                use_rms_calibration=TRUE,
                                log_file=NULL,
                                force=FALSE) {
  
  suppressPackageStartupMessages({
    library(survival)
    library(survminer)
    library(timeROC)
    library(pec)
    library(rms)
    library(ggplot2)
    library(grid)
    library(gridExtra)
    library(png)
    library(dplyr)
    library(glmnet)
  })
  
  dir_create(fig_dir); dir_create(out_tables); dir_create(out_objects); dir_create(cache_out)
  
  in_meta <- file.path(cache_in_lasso, "meta2_risk.rds")
  in_halldf <- file.path(cache_in_scores, "hall_df.rds")
  in_top5 <- file.path(cache_in_uni, "top5.rds")
  in_cv <- file.path(cache_in_lasso, "cv_pw.rds")
  
  out_cox_final <- file.path(out_objects, "Cox_final_pwRisk.rds")
  out_roc <- file.path(out_objects, "timeROC_object.rds")
  out_pec <- file.path(out_objects, "Calibration_object_pec.rds")
  out_cox_rms <- file.path(out_objects, "Cox_rms_grade_adjusted.rds")
  
  km_png <- file.path(fig_dir, "KM_mgroup.png")
  roc_png <- file.path(fig_dir, "PANEL_timeROC.png")
  cal_pec_png <- file.path(fig_dir, "PANEL_calibration_PEC_1y3y5y.png")
  cal_rms_png <- file.path(fig_dir, "PANEL_calibration_RMS_1y3y5y_legendBottom.png")
  workflow_png <- file.path(fig_dir, "PANEL_workflow.png")
  summary_png <- file.path(fig_dir, "FIGURE_SummaryPanel_Molecular.png")
  
  out_auc_csv <- file.path(out_tables, "timeROC_AUC.csv")
  out_zph_txt <- file.path(out_tables, "VALID_cox_zph.txt")
  out_zph_png <- file.path(fig_dir, "VALID_cox_zph.png")
  out_lasso_stab <- file.path(out_tables, "VALID_lasso_selection_stability.csv")
  out_top30 <- file.path(out_tables, "VALID_top30_univariate_hallmarks.csv")
  out_session <- file.path(out_tables, "sessionInfo.txt")
  out_model_txt <- file.path(out_tables, "Model_summaries.txt")
  out_hr_csv <- file.path(out_tables, "Cox_final_multivariable_grade_adjusted.csv")
  out_boot_txt <- file.path(out_tables, "VALID_bootstrap_grade_adjusted.txt")
  out_boot_csv <- file.path(out_tables, "VALID_bootstrap_grade_adjusted.csv")
  
  run_step("08_eval_panels",
           c(out_cox_final, out_roc, out_auc_csv, km_png, roc_png, workflow_png, summary_png,
             out_model_txt, out_zph_txt, out_zph_png, out_lasso_stab, out_top30, out_session,
             out_hr_csv, out_boot_txt, out_boot_csv),
           force=force, log_file=log_file,
           fn_build=function() {
             
             meta2 <- read_rds(in_meta, log_file)
             hall_df <- read_rds(in_halldf, log_file)
             top5 <- read_rds(in_top5, log_file)
             cv_pw <- read_rds(in_cv, log_file)
             
             # --- KM (mgroup) ---
             pal_group <- c("Low"="#2C7BB6","High"="#D7191C")
             fit_km <- survfit(Surv(time, event) ~ mgroup, data=meta2)
             p_km <- ggsurvplot(
               fit_km, data=meta2, pval=TRUE, risk.table=TRUE,
               palette=c(pal_group["Low"], pal_group["High"]),
               legend.labs=c("Low","High"),
               title="Molecular group (Hallmark-pathway risk): Overall survival"
             )
             ggsave(km_png, p_km$plot, width=7, height=6, dpi=300)
             
             # --- Cox adjusted (grade-aware if available) ---
             # prefer grade_bin if non-missing
             has_grade <- "grade_bin" %in% colnames(meta2) && any(!is.na(meta2$grade_bin))
             if (has_grade) {
               cox_final <- coxph(Surv(time, event) ~ pw_risk + age_years + stage_bin + grade_bin, data=meta2)
             } else {
               cox_final <- coxph(Surv(time, event) ~ pw_risk + age_years + stage_bin, data=meta2)
             }
             saveRDS(cox_final, out_cox_final)
             
             # --- timeROC ---
             times <- c(365,1095,1825)
             roc_pw <- timeROC(T=meta2$time, delta=meta2$event, marker=meta2$pw_risk, cause=1, times=times, iid=TRUE)
             write.csv(data.frame(time=times, AUC=roc_pw$AUC), out_auc_csv, row.names=FALSE)
             saveRDS(roc_pw, out_roc)
             
             png(roc_png, width=900, height=750, res=150)
             plot(roc_pw, time=365,  col=1, lwd=2, lty=1)
             plot(roc_pw, time=1095, col=2, lwd=2, lty=2, add=TRUE)
             plot(roc_pw, time=1825, col=3, lwd=2, lty=3, add=TRUE)
             legend("bottomright",
                    legend=c(paste0("1y AUC=", round(roc_pw$AUC[1],2)),
                             paste0("3y AUC=", round(roc_pw$AUC[2],2)),
                             paste0("5y AUC=", round(roc_pw$AUC[3],2))),
                    col=c(1,2,3), lty=c(1,2,3), lwd=2, bty="n")
             dev.off()
             
    
             # --- PEC calibration (FIXED: avoid 'fml' symbol in model call) ---
             times <- c(365, 1095, 1825)
             
             cal_cols <- c("time","event","pw_risk","age_years","stage_bin")
             if (has_grade) cal_cols <- c(cal_cols, "grade_bin")
             
             cal_dat <- meta2[, cal_cols]
             cal_dat <- cal_dat[complete.cases(cal_dat) & cal_dat$time > 0, , drop=FALSE]
             
             fml <- if (has_grade)
               Surv(time,event) ~ pw_risk + age_years + stage_bin + grade_bin
             else
               Surv(time,event) ~ pw_risk + age_years + stage_bin
             
             # ✅ کلید اصلی: do.call باعث می‌شود فرمول به صورت مقدار داخل call ذخیره شود
             coxph_fit <- do.call(survival::coxph, list(
               formula = fml,
               data    = cal_dat,
               x       = TRUE,
               y       = TRUE,
               model   = TRUE
             ))
             
             pec_fit <- pec::pec(
               object      = list("Cox" = coxph_fit),
               formula     = Surv(time,event) ~ 1,
               data        = cal_dat,
               times       = times,
               cens.model  = "cox",
               splitMethod = "Boot632",
               B           = 500   # یا 200، پایین توضیح دادم
             )
             
             # --- RMS calibration + bootstrap (optional) ---
             if (use_rms_calibration) {
               meta_rms <- meta2 %>%
                 dplyr::select(time, event, pw_risk, age_years, stage_bin, dplyr::any_of("grade_bin")) %>%
                 dplyr::filter(complete.cases(.), time > 0)
               
               meta_rms$stage_bin <- factor(meta_rms$stage_bin, levels=c("Early","Late"))
               if ("grade_bin" %in% colnames(meta_rms)) meta_rms$grade_bin <- factor(meta_rms$grade_bin, levels=c("Low","High"))
               
               # --- RMS calibration + bootstrap (fixed dd scope) ---
               meta_rms <- meta2 %>%
                 dplyr::select(time, event, pw_risk, age_years, stage_bin, dplyr::any_of("grade_bin")) %>%
                 dplyr::filter(complete.cases(.), time > 0)
               
               meta_rms$stage_bin <- factor(meta_rms$stage_bin, levels=c("Early","Late"))
               if ("grade_bin" %in% colnames(meta_rms)) meta_rms$grade_bin <- factor(meta_rms$grade_bin, levels=c("Low","High"))
               
               dd_local <- rms::datadist(meta_rms)
               
               # ✅ کلید: dd باید با همین نام برای rms قابل پیدا کردن باشد
               assign("dd", dd_local, envir = .GlobalEnv)
               
               old_opt <- options(datadist = "dd")
               on.exit({
                 options(old_opt)
                 # اگر دوست داری dd بعد از اتمام پاک شود:
                 # if (exists("dd", envir=.GlobalEnv)) rm("dd", envir=.GlobalEnv)
               }, add = TRUE)
               
               fml_rms <- if ("grade_bin" %in% colnames(meta_rms))
                 Surv(time,event) ~ pw_risk + age_years + stage_bin + grade_bin
               else
                 Surv(time,event) ~ pw_risk + age_years + stage_bin
               
               cox_rms <- do.call(rms::cph, list(
                 formula = fml_rms,
                 data = meta_rms,
                 x = TRUE, y = TRUE, surv = TRUE
               ))
               
               set.seed(1)
               val <- rms::validate(cox_rms, method="boot", B=500)
               
               set.seed(1)
               cal1 <- rms::calibrate(cox_rms, method="boot", u=365,  B=500)
               cal3 <- rms::calibrate(cox_rms, method="boot", u=1095, B=500)
               cal5 <- rms::calibrate(cox_rms, method="boot", u=1825, B=500)
               saveRDS(cox_rms, out_cox_rms)
               
               set.seed(1)
               val <- validate(cox_rms, method="boot", B=500)
               
               capture.output(val, file=out_boot_txt)
               val_mat <- as.matrix(val)
               out_mat <- cbind(Metric=rownames(val_mat), val_mat)
               write.table(out_mat, file=out_boot_csv, sep=",", row.names=FALSE, col.names=TRUE, quote=TRUE)
               
               set.seed(1)
               cal1 <- calibrate(cox_rms, method="boot", u=365,  B=500)
               cal3 <- calibrate(cox_rms, method="boot", u=1095, B=500)
               cal5 <- calibrate(cox_rms, method="boot", u=1825, B=500)
               
               png(cal_rms_png, width=1800, height=750, res=300)
               layout(matrix(c(1,2,3,4,4,4), nrow=2, byrow=TRUE), heights=c(4.5,1.2))
               par(oma=c(3.2,5.0,0.8,0.6))
               
               plot_cal_panel <- function(cal_obj, main_title) {
                 par(mar=c(4.0,2.2,3.0,1.0))
                 plot(cal_obj, xlim=c(0,1), ylim=c(0,1),
                      xlab="", ylab="", main=main_title, subtitles=FALSE,
                      lwd=2, cex=0.90, cex.axis=0.85, cex.main=1.0, cex.lab=1.0)
                 abline(0,1,lty=2,lwd=1.5)
               }
               plot_cal_panel(cal1, "1-Year")
               plot_cal_panel(cal3, "3-Year")
               plot_cal_panel(cal5, "5-Year")
               mtext("Predicted survival probability", side=1, outer=TRUE, line=1.4, cex=1.0)
               mtext("Observed survival probability",  side=2, outer=TRUE, line=3.0, cex=1.0)
               
               par(mar=c(0,0,0,0)); plot.new()
               legend("center",
                      legend=c("Apparent","Bias-corrected","Bootstrap samples","Ideal"),
                      col=c("black","blue","grey70","black"),
                      lty=c(1,1,1,2), lwd=c(2,2,2,1.5),
                      bty="n", horiz=TRUE, cex=1.0)
               dev.off()
               layout(1); par(mfrow=c(1,1)); par(oma=c(0,0,0,0))
             }
             
             # --- Model summary + HR table ---
             sink(out_model_txt)
             cat("===== Cox adjusted model =====\n")
             print(summary(cox_final))
             sink()
             
             summ <- summary(cox_final)
             hr_table <- data.frame(
               Variable = rownames(summ$coefficients),
               HR = summ$coefficients[, "exp(coef)"],
               Lower95 = summ$conf.int[, "lower .95"],
               Upper95 = summ$conf.int[, "upper .95"],
               pvalue = summ$coefficients[, "Pr(>|z|)"]
             )
             write.csv(hr_table, out_hr_csv, row.names=FALSE)
             
             # --- zph diagnostics ---
             zph <- cox.zph(cox_final)
             capture.output(print(zph), file=out_zph_txt)
             png(out_zph_png, width=1100, height=800, res=150); plot(zph); dev.off()
             
             # --- LASSO stability (bootstrapped resampling on top5) ---
             set.seed(1)
             sel_list <- replicate(50, {
               idx <- sample(seq_len(nrow(meta2)), replace=TRUE)
               Xb <- as.matrix(hall_df[idx, top5, drop=FALSE])
               yb <- Surv(meta2$time[idx], meta2$event[idx])
               fit <- glmnet(Xb, yb, family="cox", alpha=1, lambda=cv_pw$lambda.min)
               cf <- as.matrix(coef(fit))
               rownames(cf)[cf[,1] != 0]
             }, simplify=FALSE)
             tab <- sort(table(unlist(sel_list)), decreasing=TRUE)
             write.csv(data.frame(pathway=names(tab), freq=as.integer(tab)), out_lasso_stab, row.names=FALSE)
             
             # --- top30 uni (for README) ---
             # (اگر uni_res قبلاً ذخیره شده، بهتر است از cache/05_unicox بخوانیم؛ اینجا دوباره نمی‌سازیم)
             # برای اینکه فایل همیشه ساخته شود، از cache می‌خوانیم:
             uni_res <- read_rds(file.path(cache_in_uni, "uni_res.rds"), log_file)
             write.csv(head(uni_res, 30), out_top30, row.names=FALSE)
             
             # --- sessionInfo ---
             sink(out_session); print(sessionInfo()); sink()
             
             # --- Workflow panel (clean & consistent path) ---
             p_workflow <- ggplot() + theme_void() +
               annotate("text", x=0, y=6.5, label="TCGA-KIRC (Primary Tumor)", size=6, fontface="bold") +
               annotate("label", x=0, y=5.6, label="Data download & QC\n(TCGAbiolinks)", size=4) +
               annotate("segment", x=0, xend=0, y=5.2, yend=4.8, arrow=arrow(length=unit(0.12,"in")), linewidth=0.6) +
               annotate("label", x=0, y=4.4, label="ssGSEA Hallmark scoring", size=4) +
               annotate("segment", x=0, xend=0, y=4.0, yend=3.6, arrow=arrow(length=unit(0.12,"in")), linewidth=0.6) +
               annotate("label", x=0, y=3.2, label="Cox screening + LASSO", size=4) +
               annotate("segment", x=0, xend=0, y=2.8, yend=2.4, arrow=arrow(length=unit(0.12,"in")), linewidth=0.6) +
               annotate("label", x=0, y=2.0, label="Molecular risk score\n(High / Low)", size=4) +
               annotate("segment", x=0, xend=-2, y=1.6, yend=1.0, arrow=arrow(length=unit(0.12,"in")), linewidth=0.6) +
               annotate("segment", x=0, xend=0,  y=1.6, yend=1.0, arrow=arrow(length=unit(0.12,"in")), linewidth=0.6) +
               annotate("segment", x=0, xend=2,  y=1.6, yend=1.0, arrow=arrow(length=unit(0.12,"in")), linewidth=0.6) +
               annotate("label", x=-2, y=0.6, label="DESeq2\n(adjusted)", size=3.8) +
               annotate("label", x=0,  y=0.6, label="Hallmark GSEA", size=3.8) +
               annotate("label", x=2,  y=0.6, label="Adjusted Cox model", size=3.8) +
               coord_cartesian(xlim=c(-3,3), ylim=c(0,7)) +
               theme(plot.margin=margin(20,20,20,20))
             ggsave(workflow_png, p_workflow, width=6, height=7, dpi=300)
             
             # --- Summary panel (robust: arrangeGrob + png THEN draw) ---
             req_imgs <- c("PANEL_pca.png","PANEL_volcano_labeled.png","PANEL_gsea_colored.png","KM_mgroup.png","PANEL_timeROC.png")
             miss <- req_imgs[!file.exists(file.path(fig_dir, req_imgs))]
             if (length(miss) > 0) stopf("Missing figure(s) for summary panel:\n", paste(miss, collapse="\n"), log_file=log_file)
             
             g_pca     <- rasterGrob(readPNG(file.path(fig_dir,"PANEL_pca.png")), interpolate=TRUE)
             g_volcano <- rasterGrob(readPNG(file.path(fig_dir,"PANEL_volcano_labeled.png")), interpolate=TRUE)
             g_gsea    <- rasterGrob(readPNG(file.path(fig_dir,"PANEL_gsea_colored.png")), interpolate=TRUE)
             g_km      <- rasterGrob(readPNG(file.path(fig_dir,"KM_mgroup.png")), interpolate=TRUE)
             g_roc     <- rasterGrob(readPNG(file.path(fig_dir,"PANEL_timeROC.png")), interpolate=TRUE)
             
             panel_grob <- arrangeGrob(
               g_pca,
               g_volcano, g_gsea,
               g_km,      g_roc,
               ncol=2, heights=c(1,1,1)
             )
             
             png(summary_png, width=2000, height=2600, res=300)
             grid.draw(panel_grob)
             dev.off()
           })
  
  invisible(TRUE)
}