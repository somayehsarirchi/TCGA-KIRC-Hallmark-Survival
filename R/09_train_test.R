source(here::here("R","00_utils.R"))

step_09_train_test <- function(cache_in_lasso=here::here("cache","06_lasso"),
                               cache_in_scores=here::here("cache","04_scores"),
                               fig_dir=here::here("results","figures"),
                               out_tables=here::here("results","tables"),
                               cache_out=here::here("cache","09_train_test"),
                               log_file=NULL,
                               force=FALSE) {
  
  suppressPackageStartupMessages({
    library(survival)
    library(glmnet)
    library(survminer)
    library(timeROC)
  })
  
  dir_create(fig_dir); dir_create(out_tables); dir_create(cache_out)
  
  in_meta <- file.path(cache_in_lasso, "meta2_risk.rds")
  in_halldf <- file.path(cache_in_scores, "hall_df.rds")
  
  out_perf <- file.path(out_tables, "TRAIN_TEST_summary.csv")
  out_km   <- file.path(fig_dir, "TEST_KM_mgroup.png")
  out_roc  <- file.path(fig_dir, "TEST_timeROC.png")
  out_auc  <- file.path(out_tables, "TEST_timeROC_AUC.csv")
  out_train_uni <- file.path(out_tables, "TRAIN_univariate_Hallmarks.csv")
  out_train_beta <- file.path(out_tables, "TRAIN_LASSO_selectedPathways_beta.csv")
  
  run_step("09_train_test",
           c(out_perf, out_km, out_roc, out_auc, out_train_uni, out_train_beta),
           force=force, log_file=log_file,
           fn_build=function() {
             
             meta2 <- read_rds(in_meta, log_file)
             hall_df <- read_rds(in_halldf, log_file)
             stopifnot(all(rownames(hall_df) == rownames(meta2)))
             
             set.seed(1)
             n <- nrow(meta2)
             train_idx <- sample(seq_len(n), size=floor(0.7*n))
             meta_tr <- meta2[train_idx, , drop=FALSE]
             meta_te <- meta2[-train_idx, , drop=FALSE]
             hall_tr <- hall_df[rownames(meta_tr), , drop=FALSE]
             hall_te <- hall_df[rownames(meta_te), , drop=FALSE]
             
             # Univariate Cox TRAIN only -> top5
             uni_tr <- lapply(colnames(hall_tr), function(pw){
               fit <- coxph(Surv(time,event) ~ hall_tr[[pw]], data=meta_tr)
               s <- summary(fit)
               data.frame(pathway=pw, HR=s$coef[1,"exp(coef)"], p=s$coef[1,"Pr(>|z|)"], stringsAsFactors=FALSE)
             })
             uni_tr <- do.call(rbind, uni_tr)
             uni_tr$padj <- p.adjust(uni_tr$p, method="BH")
             uni_tr <- uni_tr[order(uni_tr$padj), ]
             write.csv(uni_tr, out_train_uni, row.names=FALSE)
             top5_tr <- head(uni_tr$pathway, 5)
             
             # LASSO TRAIN only
             X_tr <- as.matrix(hall_tr[, top5_tr, drop=FALSE])
             y_tr <- Surv(meta_tr$time, meta_tr$event)
             set.seed(1)
             cv_tr <- cv.glmnet(X_tr, y_tr, family="cox", alpha=1, nfolds=10)
             fit_tr <- glmnet(X_tr, y_tr, family="cox", alpha=1, lambda=cv_tr$lambda.min)
             
             coef_tr <- as.matrix(coef(fit_tr))
             sel_tr <- rownames(coef_tr)[coef_tr[,1] != 0]
             if (length(sel_tr) == 0) stopf("TRAIN LASSO selected 0 pathways.", log_file=log_file)
             
             beta_tr <- coef_tr[sel_tr, 1, drop=TRUE]
             write.csv(data.frame(pathway=sel_tr, beta=as.numeric(beta_tr)), out_train_beta, row.names=FALSE)
             
             # risk score for train/test using TRAIN betas
             meta_tr$pw_risk_tr <- as.numeric(as.matrix(hall_tr[, sel_tr, drop=FALSE]) %*% beta_tr)
             meta_te$pw_risk_tr <- as.numeric(as.matrix(hall_te[, sel_tr, drop=FALSE]) %*% beta_tr)
             
             cut_tr <- median(meta_tr$pw_risk_tr, na.rm=TRUE)
             meta_tr$mgroup_tr <- factor(ifelse(meta_tr$pw_risk_tr >= cut_tr, "High","Low"), levels=c("Low","High"))
             meta_te$mgroup_tr <- factor(ifelse(meta_te$pw_risk_tr >= cut_tr, "High","Low"), levels=c("Low","High"))
             
             # KM on TEST
             pal_group <- c("Low"="#2C7BB6","High"="#D7191C")
             fit_km_te <- survfit(Surv(time,event) ~ mgroup_tr, data=meta_te)
             p_km_te <- ggsurvplot(fit_km_te, data=meta_te, pval=TRUE, risk.table=TRUE,
                                   palette=c(pal_group["Low"], pal_group["High"]),
                                   legend.labs=c("Low","High"),
                                   title="TEST set: Molecular group (trained threshold)")
             ggsave(out_km, p_km_te$plot, width=7, height=6, dpi=300)
             
             # timeROC on TEST
             times <- c(365,1095,1825)
             roc_te <- timeROC(T=meta_te$time, delta=meta_te$event, marker=meta_te$pw_risk_tr, cause=1, times=times, iid=TRUE)
             write.csv(data.frame(time=times, AUC=roc_te$AUC), out_auc, row.names=FALSE)
             
             png(out_roc, width=900, height=750, res=150)
             plot(roc_te, time=365,  col=1, lwd=2, lty=1)
             plot(roc_te, time=1095, col=2, lwd=2, lty=2, add=TRUE)
             plot(roc_te, time=1825, col=3, lwd=2, lty=3, add=TRUE)
             legend("bottomright",
                    legend=c(paste0("1y AUC=", round(roc_te$AUC[1],2)),
                             paste0("3y AUC=", round(roc_te$AUC[2],2)),
                             paste0("5y AUC=", round(roc_te$AUC[3],2))),
                    col=c(1,2,3), lty=c(1,2,3), lwd=2, bty="n")
             dev.off()
             
             perf <- data.frame(
               set=c("train","test"),
               n=c(nrow(meta_tr), nrow(meta_te)),
               events=c(sum(meta_tr$event), sum(meta_te$event)),
               event_rate=c(mean(meta_tr$event), mean(meta_te$event)),
               n_low=c(sum(meta_tr$mgroup_tr=="Low"), sum(meta_te$mgroup_tr=="Low")),
               n_high=c(sum(meta_tr$mgroup_tr=="High"), sum(meta_te$mgroup_tr=="High"))
             )
             write.csv(perf, out_perf, row.names=FALSE)
           })
  
  invisible(TRUE)
}