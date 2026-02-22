# =====================================================
# 01_download_preprocess.R
# Download TCGA-KIRC STAR counts + basic preprocessing
# Outputs:
#   data_cache/G_raw.rds
#   data_cache/counts_mat.rds
#   data_cache/meta_raw.rds
#   data_cache/gene_annot.rds
#   results/figures/QC_G.png
# =====================================================

source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(ggplot2)
})

QG_path <- file.path(DIR$cache, "QG_TCGA-KIRC_STARcounts.rds")
G_path  <- file.path(DIR$cache, "G_TCGA-KIRC_STARcounts.rds")

if (!file.exists(QG_path)) {
  message("[01] Querying GDC...")
  QG <- GDCquery(
    project = "TCGA-KIRC",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor"
  )
  saveRDS(QG, QG_path)
} else {
  message("[01] Loading cached query object...")
  QG <- readRDS(QG_path)
}

if (!file.exists(G_path)) {
  message("[01] Downloading + preparing SummarizedExperiment...")
  GDCdownload(QG)
  G <- GDCprepare(QG, summarizedExperiment = TRUE)
  saveRDS(G, G_path)
} else {
  message("[01] Loading cached SummarizedExperiment...")
  G <- readRDS(G_path)
}

# remove NA stage
G <- G[, !is.na(colData(G)$ajcc_pathologic_stage)]
saveRDS(G, file.path(DIR$cache, "G_raw.rds"))

# preprocessing
qc_png <- file.path(DIR$fig, "QC_G.png")
message("[01] Running TCGAanalyze_Preprocessing...")
Gpre <- TCGAanalyze_Preprocessing(
  G, cor.cut = 0.75,
  datatype = "unstranded",
  filename = qc_png
)

counts_mat <- as.matrix(Gpre)

meta <- as.data.frame(colData(G))
meta <- meta[colnames(counts_mat), , drop=FALSE]
stopifnot(all(rownames(meta) == colnames(counts_mat)))

gene_annot <- as.data.frame(rowData(G))
gene_annot$gene_id <- rownames(gene_annot)
gene_annot <- gene_annot[match(rownames(counts_mat), gene_annot$gene_id), , drop=FALSE]

if (!("gene_name" %in% colnames(gene_annot))) {
  stop("[01] rowData(G) does not contain 'gene_name'. Please check correct column name.")
}

saveRDS(counts_mat, file.path(DIR$cache, "counts_mat.rds"))
saveRDS(meta,       file.path(DIR$cache, "meta_raw.rds"))
saveRDS(gene_annot, file.path(DIR$cache, "gene_annot.rds"))

message("[01] DONE: counts_mat/meta/gene_annot saved in data_cache/")
