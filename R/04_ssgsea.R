source(here::here("R","00_utils.R"))

step_04_ssgsea <- function(cache_in=here::here("cache","03_vst"),
                           cache_out=here::here("cache","04_scores"),
                           log_file=NULL,
                           force=FALSE) {
  
  suppressPackageStartupMessages({
    library(msigdbr)
    library(dplyr)
    library(GSVA)
  })
  
  dir_create(cache_out)
  
  in_matgn <- file.path(cache_in, "mat_gn.rds")
  out_path <- file.path(cache_out, "pathways_hallmark.rds")
  out_halldf <- file.path(cache_out, "hall_df.rds")
  
  run_step("04_ssgsea", c(out_path, out_halldf), force=force, log_file=log_file, fn_build=function() {
    
    mat_gn <- read_rds(in_matgn, log_file)
    hallmark <- msigdbr(species="Homo sapiens", category="H") %>%
      dplyr::select(gs_name, gene_symbol)
    pathways <- split(hallmark$gene_symbol, hallmark$gs_name)
    
    hall_scores <- tryCatch({
      par <- GSVA::ssgseaParam(exprData=as.matrix(mat_gn), geneSets=pathways)
      GSVA::gsva(par)
    }, error=function(e){
      GSVA::gsva(as.matrix(mat_gn), pathways, method="ssgsea", kcdf="Gaussian", abs.ranking=TRUE)
    })
    
    hall_df <- as.data.frame(t(hall_scores)) # samples x pathways
    assert_nonempty(hall_df, "hall_df", log_file)
    
    save_rds(pathways, out_path, log_file)
    save_rds(hall_df,  out_halldf, log_file)
  })
  
  invisible(TRUE)
}