source(here::here("R","00_utils.R"))

step_02_preprocess <- function(cache_in=here::here("cache","01_gdc"),
                               cache_out=here::here("cache","02_preprocess"),
                               fig_dir=here::here("results","figures"),
                               log_file=NULL,
                               force=FALSE) {
  
  suppressPackageStartupMessages({
    library(TCGAbiolinks)
    library(SummarizedExperiment)
    library(dplyr)
  })
  
  dir_create(cache_out)
  dir_create(fig_dir)
  
  g_path <- file.path(cache_in, "G_TCGA-KIRC_STARcounts.rds")
  
  out_counts_mat <- file.path(cache_out, "counts_mat.rds")
  out_meta2      <- file.path(cache_out, "meta2.rds")
  out_counts2    <- file.path(cache_out, "counts2.rds")
  out_gene_annot <- file.path(cache_out, "gene_annot.rds")
  out_qc_png     <- file.path(cache_out, "qc_png_path.rds")
  
  run_step("02_preprocess",
           c(out_counts_mat, out_meta2, out_counts2, out_gene_annot, out_qc_png),
           force=force, log_file=log_file,
           fn_build=function() {
             
             G <- read_rds(g_path, log_file)
             
             # remove NA stage
             if (!("ajcc_pathologic_stage" %in% colnames(colData(G)))) {
               stopf("colData(G) missing ajcc_pathologic_stage", log_file=log_file)
             }
             G <- G[, !is.na(colData(G)$ajcc_pathologic_stage)]
             
             qc_png <- file.path(fig_dir, "QC_G.png")
             Gpre <- TCGAanalyze_Preprocessing(
               G, cor.cut=0.75, datatype="unstranded", filename=qc_png
             )
             
             counts_mat <- as.matrix(Gpre)
             assert_nonempty(counts_mat, "counts_mat", log_file)
             
             meta <- as.data.frame(colData(G))
             meta <- meta[colnames(counts_mat), , drop=FALSE]
             if (!all(rownames(meta) == colnames(counts_mat))) stopf("Meta not aligned to counts_mat", log_file=log_file)
             
             gene_annot <- as.data.frame(rowData(G))
             gene_annot$gene_id <- rownames(gene_annot)
             gene_annot <- gene_annot[match(rownames(counts_mat), gene_annot$gene_id), , drop=FALSE]
             
             if (!("gene_name" %in% colnames(gene_annot))) {
               stopf("rowData(G) does not contain 'gene_name' column", log_file=log_file)
             }
             
             # survival + covariates
             meta$time  <- ifelse(!is.na(meta$days_to_death), meta$days_to_death, meta$days_to_last_follow_up)
             meta$event <- ifelse(meta$vital_status == "Dead", 1, 0)
             
             meta$age_at_diagnosis <- suppressWarnings(as.numeric(meta$age_at_diagnosis))
             meta$age_years <- meta$age_at_diagnosis / 365.25
             
             stage_raw <- gsub("^Stage\\s+", "", as.character(meta$ajcc_pathologic_stage))
             stage_raw <- toupper(stage_raw)
             meta$stage_bin <- NA_character_
             meta$stage_bin[stage_raw %in% c("I","IA","IB","II","IIA","IIB","IIC")] <- "Early"
             meta$stage_bin[stage_raw %in% c("III","IIIA","IIIB","IIIC","IV","IVA","IVB","IVC")] <- "Late"
             meta$stage_bin <- factor(meta$stage_bin, levels=c("Early","Late"))
             
             # grade_bin robust
             meta <- make_grade_bin(meta, log_file=log_file)
             
             keep_surv <- !is.na(meta$time) & meta$time > 0 &
               !is.na(meta$event) &
               !is.na(meta$age_years) &
               !is.na(meta$stage_bin)
             
             meta2   <- meta[keep_surv, , drop=FALSE]
             counts2 <- counts_mat[, rownames(meta2), drop=FALSE]
             if (!all(colnames(counts2) == rownames(meta2))) stopf("counts2 not aligned to meta2", log_file=log_file)
             
             assert_nonempty(meta2, "meta2", log_file)
             assert_nonempty(counts2, "counts2", log_file)
             
             save_rds(counts_mat, out_counts_mat, log_file)
             save_rds(meta2,      out_meta2,      log_file)
             save_rds(counts2,    out_counts2,    log_file)
             save_rds(gene_annot, out_gene_annot, log_file)
             save_rds(qc_png,     out_qc_png,     log_file)
           })
  
  invisible(TRUE)
}