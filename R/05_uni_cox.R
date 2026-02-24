source(here::here("R","00_utils.R"))

step_05_uni_cox <- function(cache_in_meta=here::here("cache","02_preprocess"),
                            cache_in_scores=here::here("cache","04_scores"),
                            cache_out=here::here("cache","05_unicox"),
                            out_tables=here::here("results","tables"),
                            log_file=NULL,
                            force=FALSE) {
  
  suppressPackageStartupMessages({
    library(survival)
  })
  
  dir_create(cache_out); dir_create(out_tables)
  
  in_meta2 <- file.path(cache_in_meta, "meta2.rds")
  in_halldf <- file.path(cache_in_scores, "hall_df.rds")
  
  out_uni <- file.path(cache_out, "uni_res.rds")
  out_top5 <- file.path(cache_out, "top5.rds")
  out_csv <- file.path(out_tables, "Cox_univariate_Hallmarks.csv")
  
  run_step("05_uni_cox", c(out_uni, out_top5, out_csv), force=force, log_file=log_file, fn_build=function() {
    
    meta2 <- read_rds(in_meta2, log_file)
    hall_df <- read_rds(in_halldf, log_file)
    
    if (!all(rownames(hall_df) == rownames(meta2))) stopf("hall_df rownames not aligned to meta2", log_file=log_file)
    
    uni_res <- lapply(colnames(hall_df), function(pw){
      fit <- coxph(Surv(time, event) ~ hall_df[[pw]], data=meta2)
      s <- summary(fit)
      data.frame(pathway=pw,
                 HR=s$coef[1,"exp(coef)"],
                 p=s$coef[1,"Pr(>|z|)"],
                 stringsAsFactors=FALSE)
    })
    uni_res <- do.call(rbind, uni_res)
    uni_res$padj <- p.adjust(uni_res$p, method="BH")
    uni_res <- uni_res[order(uni_res$padj), ]
    top5 <- head(uni_res$pathway, 5)
    
    save_rds(uni_res, out_uni, log_file)
    save_rds(top5, out_top5, log_file)
    write.csv(uni_res, out_csv, row.names=FALSE)
  })
  
  invisible(TRUE)
}