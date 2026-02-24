source(here::here("R","00_utils.R"))

step_03_vst_symbol <- function(cache_in=here::here("cache","02_preprocess"),
                               cache_out=here::here("cache","03_vst"),
                               log_file=NULL,
                               force=FALSE) {
  
  suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
  })
  
  dir_create(cache_out)
  
  in_counts2 <- file.path(cache_in, "counts2.rds")
  in_meta2   <- file.path(cache_in, "meta2.rds")
  in_annot   <- file.path(cache_in, "gene_annot.rds")
  
  out_matgn  <- file.path(cache_out, "mat_gn.rds")
  
  run_step("03_vst_symbol", c(out_matgn), force=force, log_file=log_file, fn_build=function() {
    
    counts2    <- read_rds(in_counts2, log_file)
    meta2      <- read_rds(in_meta2,   log_file)
    gene_annot <- read_rds(in_annot,   log_file)
    
    dds0 <- DESeqDataSetFromMatrix(countData=round(counts2), colData=meta2, design=~1)
    dds0 <- dds0[rowSums(counts(dds0) >= 10) >= 10, ]
    vsd0 <- vst(dds0, blind=TRUE)
    mat_vst <- assay(vsd0)
    
    gene_annot2 <- gene_annot[match(rownames(mat_vst), gene_annot$gene_id), , drop=FALSE]
    gene_names <- gene_annot2$gene_name
    
    mat_gn <- mat_vst
    rownames(mat_gn) <- gene_names
    mat_gn <- mat_gn[!is.na(rownames(mat_gn)) & rownames(mat_gn)!="", , drop=FALSE]
    
    if (any(duplicated(rownames(mat_gn)))) {
      tmp <- as.data.frame(mat_gn)
      tmp$gene <- rownames(tmp)
      tmp <- tmp %>% group_by(gene) %>% summarise(across(where(is.numeric), mean), .groups="drop")
      mat_gn <- as.matrix(tmp[, -1])
      rownames(mat_gn) <- tmp$gene
    }
    
    if (any(duplicated(rownames(mat_gn)))) stopf("Duplicated gene symbols remain after aggregation", log_file=log_file)
    
    save_rds(mat_gn, out_matgn, log_file)
  })
  
  invisible(TRUE)
}