# =====================================================
# 02_ssgsea_scoring.R
# Build survival-ready cohort + ssGSEA Hallmark scores
# Outputs:
#   data_cache/meta2_survival.rds
#   data_cache/counts2_survival.rds
#   data_cache/mat_gn_vst.rds
#   data_cache/hall_df_scores.rds
#   data_cache/hallmark_pathways.rds
# =====================================================

source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(msigdbr)
  library(GSVA)
  library(survival)
})

counts_mat <- readRDS(file.path(DIR$cache, "counts_mat.rds"))
meta       <- readRDS(file.path(DIR$cache, "meta_raw.rds"))
gene_annot <- readRDS(file.path(DIR$cache, "gene_annot.rds"))

# -------------------------
# 1) Build survival fields
# -------------------------
meta$time  <- ifelse(!is.na(meta$days_to_death), meta$days_to_death, meta$days_to_last_follow_up)
meta$event <- ifelse(meta$vital_status == "Dead", 1, 0)

meta$age_at_diagnosis <- suppressWarnings(as.numeric(meta$age_at_diagnosis))
meta$age_years <- meta$age_at_diagnosis / 365.25

# -------------------------
# 2) Stage binning (robust)
# -------------------------
stage_raw <- as.character(meta$ajcc_pathologic_stage)
stage_raw <- toupper(trimws(gsub("^STAGE\\s+", "", stage_raw)))

# Keep only roman prefix: I, II, III, IV (ignore suffix A/B/C)
stage_prefix <- gsub("^(I{1,3}V?|IV).*$", "\\1", stage_raw)  # e.g. "IIIA" -> "III"
stage_prefix[stage_prefix == stage_raw & !stage_raw %in% c("I","II","III","IV")] <- NA

meta$stage_bin <- NA_character_
meta$stage_bin[stage_prefix %in% c("I","II")] <- "Early"
meta$stage_bin[stage_prefix %in% c("III","IV")] <- "Late"
meta$stage_bin <- factor(meta$stage_bin, levels = c("Early","Late"))

# -------------------------
# 3) Tumor grade binning (G1/G2 vs G3/G4)
# -------------------------
# NOTE: some TCGA fields can be empty or "GX"
if ("tumor_grade" %in% colnames(meta)) {
  tg <- toupper(trimws(as.character(meta$tumor_grade)))
  tg[tg %in% c("", "NA", "N/A")] <- NA
  meta$tumor_grade <- tg

  meta$grade_bin <- dplyr::case_when(
    meta$tumor_grade %in% c("G1","G2") ~ "Low",
    meta$tumor_grade %in% c("G3","G4") ~ "High",
    TRUE ~ NA_character_
  )
  meta$grade_bin <- factor(meta$grade_bin, levels = c("Low","High"))
} else {
  meta$grade_bin <- factor(NA_character_, levels = c("Low","High"))
  warning("[02] Column 'tumor_grade' not found in meta; grade_bin will be NA.")
}

# -------------------------
# 4) Survival-ready cohort
# -------------------------
keep_surv <- !is.na(meta$time) & meta$time > 0 &
  !is.na(meta$event) &
  !is.na(meta$age_years) &
  !is.na(meta$stage_bin)

meta2   <- meta[keep_surv, , drop = FALSE]
counts2 <- counts_mat[, rownames(meta2), drop = FALSE]
stopifnot(all(colnames(counts2) == rownames(meta2)))

message("[02] Survival-ready samples: ", ncol(counts2))
message("[02] Events: ", sum(meta2$event))
message("[02] Stage table:")
print(table(meta2$stage_bin, useNA = "ifany"))
message("[02] grade_bin NA count (not filtered here): ", sum(is.na(meta2$grade_bin)))

# -------------------------
# 5) VST for ssGSEA (gene-level)
# -------------------------
dds0 <- DESeqDataSetFromMatrix(
  countData = round(counts2),
  colData   = meta2,
  design    = ~ 1
)

# keep genes with >=10 counts in >=10 samples
dds0 <- dds0[rowSums(counts(dds0) >= 10) >= 10, ]
message("[02] Genes after count filter: ", nrow(dds0))

vsd0 <- vst(dds0, blind = TRUE)
mat_vst <- assay(vsd0)

# -------------------------
# 6) Map gene symbols + de-duplicate
# -------------------------
# align annotation to current gene IDs
gene_annot2 <- gene_annot[match(rownames(mat_vst), gene_annot$gene_id), , drop = FALSE]

mat_gn <- mat_vst
mat_gn_sym <- gene_annot2$gene_name
mat_gn_sym <- as.character(mat_gn_sym)

# remove empty/NA symbols
keep_sym <- !is.na(mat_gn_sym) & mat_gn_sym != ""
mat_gn <- mat_gn[keep_sym, , drop = FALSE]
mat_gn_sym <- mat_gn_sym[keep_sym]
rownames(mat_gn) <- mat_gn_sym

# de-duplicate gene symbols by mean
n_dup <- sum(duplicated(rownames(mat_gn)))
if (n_dup > 0) {
  message("[02] Duplicated gene symbols detected: ", n_dup, " -> aggregating by mean.")
  tmp <- as.data.frame(mat_gn)
  tmp$gene <- rownames(tmp)
  tmp <- tmp %>% group_by(gene) %>% summarise(across(where(is.numeric), mean), .groups = "drop")
  mat_gn <- as.matrix(tmp[, -1])
  rownames(mat_gn) <- tmp$gene
}
stopifnot(!any(duplicated(rownames(mat_gn))))
message("[02] Final gene symbols (unique): ", nrow(mat_gn))

# -------------------------
# 7) Hallmark gene sets
# -------------------------
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

# -------------------------
# 8) ssGSEA (GSVA compatibility)
# -------------------------
hall_scores <- tryCatch({
  # GSVA >= 2 style
  par <- GSVA::ssgseaParam(exprData = as.matrix(mat_gn), geneSets = pathways)
  GSVA::gsva(par)  # pathways x samples
}, error = function(e) {
  # older GSVA
  GSVA::gsva(as.matrix(mat_gn), pathways, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)
})

hall_df <- as.data.frame(t(hall_scores))  # samples x pathways
stopifnot(all(rownames(hall_df) == rownames(meta2)))

# -------------------------
# 9) Save
# -------------------------
saveRDS(meta2,    file.path(DIR$cache, "meta2_survival.rds"))
saveRDS(counts2,  file.path(DIR$cache, "counts2_survival.rds"))
saveRDS(mat_gn,   file.path(DIR$cache, "mat_gn_vst.rds"))
saveRDS(hall_df,  file.path(DIR$cache, "hall_df_scores.rds"))
saveRDS(pathways, file.path(DIR$cache, "hallmark_pathways.rds"))

message("[02] DONE: meta2/counts2/mat_gn/hall_df saved in data_cache/")
