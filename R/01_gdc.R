source(here::here("R","00_utils.R"))

step_01_gdc <- function(project="TCGA-KIRC",
                        input_tsv=here::here("data","KidneyRNA-Seq.tsv"),
                        cache_dir=here::here("cache","01_gdc"),
                        log_file=NULL,
                        force=FALSE) {
  
  suppressPackageStartupMessages({
    library(data.table)
    library(TCGAbiolinks)
    library(SummarizedExperiment)
  })
  
  dir_create(cache_dir)
  
  qg_path <- file.path(cache_dir, "QG_TCGA-KIRC_STARcounts.rds")
  g_path  <- file.path(cache_dir, "G_TCGA-KIRC_STARcounts.rds")
  id_path <- file.path(cache_dir, "myKidneyDataID_TCGA-KIRC.rds")
  
  run_step("01_gdc", c(qg_path, g_path, id_path), force=force, log_file=log_file, fn_build=function() {
    
    if (!file.exists(input_tsv)) stopf("Missing TSV: ", input_tsv, log_file=log_file)
    
    myKidneyData <- fread(input_tsv, data.table=FALSE)
    myKidneyData <- myKidneyData[myKidneyData$project.project_id == project, ]
    myKidneyDataID <- unique(myKidneyData$cases.submitter_id)
    
    QG <- GDCquery(
      project = project,
      barcode = myKidneyDataID,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type = "Primary Tumor"
    )
    GDCdownload(QG)
    G <- GDCprepare(QG, summarizedExperiment=TRUE)
    
    assert_nonempty(myKidneyDataID, "myKidneyDataID", log_file)
    assert_nonempty(G, "G SummarizedExperiment", log_file)
    
    save_rds(QG, qg_path, log_file)
    save_rds(G,  g_path,  log_file)
    save_rds(myKidneyDataID, id_path, log_file)
  })
  
  invisible(TRUE)
}