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
})

QG_path <- file.path(DIR$cache, "QG_TCGA-KIRC_STARcounts.rds")
G_path  <- file.path(DIR$cache, "G_TCGA-KIRC_STARcounts.rds")

# -------------------------
# 1) Query / Download / Prepare (cached)
# -------------------------
if (!file.exists(QG_path)) {
  message("[01] Querying GDC (TCGA-KIRC, Primary Tumor, STAR-Counts)...")
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

# -------------------------
# 2) Basic filters
# -------------------------
# remove NA stage (needed downstream)
G <- G[, !is.na(colData(G)$ajcc_pathologic_stage)]
saveRDS(G, file.path(DIR$cache, "G_raw.rds"))

# -------------------------
# 3) Preprocessing (TCGAbiolinks QC + correlation filter)
# -------------------------
qc_png <- file.path(DIR$fig, "QC_G.png")
message("[01] Running TCGAanalyze_Preprocessing (cor.cut=0.75)...")
Gpre <- TCGAanalyze_Preprocessing(
  G, cor.cut = 0.75,
  datatype = "unstranded",
  filename = qc_png
)

counts_mat <- as.matrix(Gpre)

# -------------------------
# 4) meta + gene annotation (aligned)
# -------------------------
meta <- as.data.frame(colData(G))
meta <- meta[colnames(counts_mat), , drop = FALSE]
stopifnot(all(rownames(meta) == colnames(counts_mat)))

gene_annot <- as.data.frame(rowData(G))
gene_annot$gene_id <- rownames(gene_annot)
gene_annot <- gene_annot[match(rownames(counts_mat), gene_annot$gene_id), , drop = FALSE]

# ---- robust gene symbol column detection ----
candidate_cols <- c("gene_name", "external_gene_name", "gene", "gene_symbol", "Symbol", "symbol")
hit <- candidate_cols[candidate_cols %in% colnames(gene_annot)]

if (length(hit) == 0) {
  stop(paste0(
    "[01] Cannot find a gene symbol column in rowData(G). Available columns:\n",
    paste(colnames(gene_annot), collapse = ", ")
  ))
} else if (!("gene_name" %in% colnames(gene_annot))) {
  gene_annot$gene_name <- gene_annot[[hit[1]]]
  message("[01] NOTE: gene_name was not present; using '", hit[1], "' as gene_name.")
}

# -------------------------
# 5) Sanity checks
# -------------------------
stopifnot(!any(duplicated(colnames(counts_mat))))
stopifnot(!any(duplicated(rownames(counts_mat))))
stopifnot(nrow(meta) == ncol(counts_mat))
stopifnot(nrow(gene_annot) == nrow(counts_mat))

# -------------------------
# 6) Save
# -------------------------
saveRDS(counts_mat, file.path(DIR$cache, "counts_mat.rds"))
saveRDS(meta,       file.path(DIR$cache, "meta_raw.rds"))
saveRDS(gene_annot, file.path(DIR$cache, "gene_annot.rds"))

message("[01] DONE: counts_mat/meta/gene_annot saved in data_cache/")
message("[01] Samples: ", ncol(counts_mat), " | Genes: ", nrow(counts_mat))
