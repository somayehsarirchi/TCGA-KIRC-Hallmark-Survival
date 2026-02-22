# =====================================================
# 04_DEG_GSEA.R
# DESeq2 High vs Low (adjusted), PCA, Volcano, Hallmark GSEA
# Outputs:
#   results/tables/DESeq2_mgroup_adjusted_fullResults.csv
#   results/tables/DESeq2_mgroup_significant_genes.csv
#   results/tables/GSEA_Hallmark_mgroup.csv
#   results/figures/PANEL_pca.png
#   results/figures/PANEL_volcano_labeled.png
#   results/figures/PANEL_gsea_colored.png
# =====================================================

source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(fgsea)
  library(msigdbr)
})

set.seed(1)

meta2      <- readRDS(file.path(DIR$cache, "meta2_survival.rds"))
counts2    <- readRDS(file.path(DIR$cache, "counts2_survival.rds"))
gene_annot <- readRDS(file.path(DIR$cache, "gene_annot.rds"))

# hallmarks pathways (prefer cached, else build)
pathway_file <- file.path(DIR$cache, "hallmark_pathways.rds")
if (file.exists(pathway_file)) {
  pathways <- readRDS(pathway_file)
} else {
  message("[04] hallmark_pathways.rds not found. Building from msigdbr...")
  msig <- msigdbr(species = "Homo sapiens", category = "H")
  pathways <- split(msig$gene_symbol, msig$gs_name)
  saveRDS(pathways, pathway_file)
}

# ------------------------------
# 0) Harmonize samples meta2 <-> counts2
# ------------------------------
if (is.null(colnames(counts2))) stop("[04] counts2 has no colnames (sample IDs).")
if (is.null(rownames(meta2)))   stop("[04] meta2 has no rownames (sample IDs).")

common_ids <- intersect(rownames(meta2), colnames(counts2))
if (length(common_ids) < 50) {
  stop(sprintf("[04] Too few common samples between meta2 and counts2: %d", length(common_ids)))
}

meta2   <- meta2[common_ids, , drop=FALSE]
counts2 <- counts2[, common_ids, drop=FALSE]

# Survival sanity not strictly required for DESeq2, but mgroup is
# ------------------------------
# 1) Ensure mgroup exists (from step03)
# ------------------------------
if (!("mgroup" %in% colnames(meta2))) {
  survtab_path <- file.path(DIR$tab, "SurvivalTable_molecularRisk.csv")
  if (!file.exists(survtab_path)) {
    stop("[04] meta2 has no mgroup and SurvivalTable_molecularRisk.csv not found. Run step03 first.")
  }
  survtab <- read.csv(survtab_path, row.names=1, check.names=FALSE)
  # align via rownames
  miss <- setdiff(rownames(meta2), rownames(survtab))
  if (length(miss) > 0) {
    stop(sprintf("[04] Some meta2 samples missing in SurvivalTable_molecularRisk.csv: %d", length(miss)))
  }
  meta2$mgroup <- survtab[rownames(meta2), "mgroup"]
}

meta2$mgroup <- factor(as.character(meta2$mgroup), levels=c("Low","High"))
if (any(is.na(meta2$mgroup))) stop("[04] mgroup contains NA after alignment. Check step03 outputs.")

# ------------------------------
# 2) Build covariates robustly (stage + grade)
# ------------------------------
# stage: prefer stage_bin; else stage_4 if present; else none
stage_var <- NULL
if ("stage_bin" %in% colnames(meta2)) {
  stage_var <- "stage_bin"
  meta2$stage_bin <- factor(meta2$stage_bin, levels=c("Early","Late"))
} else if ("stage_4" %in% colnames(meta2)) {
  stage_var <- "stage_4"
  meta2$stage_4 <- factor(meta2$stage_4)
} else {
  message("[04] WARNING: No stage_bin/stage_4 found. DESeq2 will run without stage adjustment.")
}

# grade: prefer tumor_grade_bin; else tumor_grade4; else none
grade_var <- NULL
if ("tumor_grade_bin" %in% colnames(meta2)) {
  grade_var <- "tumor_grade_bin"
  meta2$tumor_grade_bin <- factor(meta2$tumor_grade_bin, levels=c("Low","High"))
} else if ("tumor_grade4" %in% colnames(meta2)) {
  grade_var <- "tumor_grade4"
  meta2$tumor_grade4 <- factor(meta2$tumor_grade4, levels=c("G1","G2","G3","G4"))
} else {
  message("[04] INFO: No tumor_grade_bin/tumor_grade4 found. DESeq2 will run without grade adjustment.")
}

# age: ensure numeric if exists
if ("age_years" %in% colnames(meta2)) meta2$age_years <- as.numeric(meta2$age_years)

# Build model formula:
# Always include mgroup; include age/stage/grade if available
covars <- c()
if ("age_years" %in% colnames(meta2)) covars <- c(covars, "age_years")
if (!is.null(stage_var))             covars <- c(covars, stage_var)
if (!is.null(grade_var))             covars <- c(covars, grade_var)
covars <- c(covars, "mgroup")

design_formula <- as.formula(paste("~", paste(covars, collapse=" + ")))
message("[04] DESeq2 design: ", deparse(design_formula))

# ------------------------------
# 3) DESeq2 setup
# ------------------------------
# Make sure counts are integer-like
count_mat <- round(as.matrix(counts2))
storage.mode(count_mat) <- "integer"

# Remove samples with NA in any design covariate
design_cols <- unique(covars)
meta_use <- meta2[, design_cols, drop=FALSE]
keep_samples <- complete.cases(meta_use)

if (sum(keep_samples) < 50) {
  stop(sprintf("[04] Too few samples after removing NA in covariates: %d", sum(keep_samples)))
}

meta2   <- meta2[keep_samples, , drop=FALSE]
count_mat <- count_mat[, rownames(meta2), drop=FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta2,
  design    = design_formula
)

# Filter genes: at least 10 counts in at least 10 samples (adjustable)
dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]
dds <- DESeq(dds)

# Main contrast: High vs Low
res <- results(dds, contrast=c("mgroup","High","Low"))
res <- res[order(res$padj), ]

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# robust join (gene_annot may have gene_id/gene_name)
if (!all(c("gene_id","gene_name") %in% colnames(gene_annot))) {
  message("[04] WARNING: gene_annot missing gene_id/gene_name columns. Writing results without gene_name.")
  res_annot <- res_df
} else {
  res_annot <- res_df %>%
    left_join(gene_annot %>% dplyr::select(gene_id, gene_name),
              by = c("gene_id"="gene_id")) %>%
    arrange(padj)
}

write.csv(res_annot, file.path(DIR$tab, "DESeq2_mgroup_adjusted_fullResults.csv"), row.names=FALSE)

# Significant genes (thresholds adjustable)
padj_cut <- 0.05
lfc_cut  <- 1

sig_annot <- res_annot %>%
  filter(!is.na(padj) & padj < padj_cut & !is.na(log2FoldChange) & abs(log2FoldChange) >= lfc_cut)

write.csv(sig_annot, file.path(DIR$tab, "DESeq2_mgroup_significant_genes.csv"), row.names=FALSE)

message("[04] Significant DEGs: ", nrow(sig_annot))

# ------------------------------
# 4) PCA
# ------------------------------
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup="mgroup", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData, aes(PC1, PC2, color=mgroup)) +
  geom_point(size=3, alpha=0.85) +
  theme_classic() +
  labs(
    title="PCA - Molecular group (High vs Low)",
    x=paste0("PC1: ", percentVar[1], "%"),
    y=paste0("PC2: ", percentVar[2], "%")
  ) +
  theme(plot.title = element_text(hjust=0.5))

# if save_plot exists in 00_config.R we use it; else fallback to ggsave
if (exists("save_plot")) {
  save_plot(p_pca, "PANEL_pca.png", w=7, h=5)
} else {
  ggsave(file.path(DIR$fig, "PANEL_pca.png"), p_pca, width=7, height=5, dpi=300)
}

# ------------------------------
# 5) Volcano (labeled)
# ------------------------------
if (!requireNamespace("ggrepel", quietly=TRUE)) {
  stop("Install ggrepel: install.packages('ggrepel')")
}

vol <- res_annot

# Safe -log10(padj): handle padj==0
vol$padj2 <- vol$padj
vol$padj2[!is.na(vol$padj2) & vol$padj2 == 0] <- .Machine$double.xmin

vol$neglog10padj <- -log10(vol$padj2)

vol$cls <- "NotSig"
vol$cls[!is.na(vol$padj) & vol$padj < padj_cut & vol$log2FoldChange >  lfc_cut] <- "Up"
vol$cls[!is.na(vol$padj) & vol$padj < padj_cut & vol$log2FoldChange < -lfc_cut] <- "Down"

top_up <- vol %>% filter(cls=="Up", !is.na(gene_name), gene_name!="") %>% arrange(padj) %>% slice_head(n=5)
top_dn <- vol %>% filter(cls=="Down", !is.na(gene_name), gene_name!="") %>% arrange(padj) %>% slice_head(n=5)
lab_df <- bind_rows(top_up, top_dn)

p_volcano <- ggplot(vol, aes(log2FoldChange, neglog10padj, color=cls)) +
  geom_point(alpha=0.6, size=1.1) +
  ggrepel::geom_text_repel(
    data=lab_df, aes(label=gene_name),
    size=4, max.overlaps=Inf, box.padding=0.4, point.padding=0.2
  ) +
  theme_classic() +
  labs(
    title="Volcano: High vs Low (adjusted)",
    x="log2FC (High - Low)",
    y="-log10(padj)"
  ) +
  scale_color_manual(values=c(Down="blue", NotSig="grey70", Up="red")) +
  theme(plot.title = element_text(hjust=0.5),
        legend.title = element_blank())

ggsave(file.path(DIR$fig, "PANEL_volcano_labeled.png"),
       p_volcano, width=7, height=5.2, dpi=300)

# ------------------------------
# 6) GSEA (Hallmark)
# ------------------------------
# Use DESeq2 "stat" (Wald statistic) as ranking metric
tmp <- vol[, intersect(c("gene_name","stat"), colnames(vol)), drop=FALSE]
if (!all(c("gene_name","stat") %in% colnames(tmp))) {
  stop("[04] Missing gene_name or stat for GSEA. Check gene_annot join and DESeq2 results.")
}

tmp <- tmp[!is.na(tmp$gene_name) & tmp$gene_name!="" & !is.na(tmp$stat), , drop=FALSE]

# Remove duplicated gene symbols (keep the one with strongest absolute stat)
tmp <- tmp %>%
  arrange(desc(abs(stat))) %>%
  distinct(gene_name, .keep_all=TRUE)

ranks <- tmp$stat
names(ranks) <- tmp$gene_name
ranks <- sort(ranks, decreasing=TRUE)

fg <- fgseaMultilevel(pathways=pathways, stats=ranks, minSize=15, maxSize=500)
fg <- fg[order(fg$padj), ]
fg_df <- as.data.frame(fg)

if ("leadingEdge" %in% names(fg_df)) {
  fg_df$leadingEdge <- vapply(fg_df$leadingEdge, function(x) paste(x, collapse=";"), FUN.VALUE=character(1))
}

write.csv(fg_df, file.path(DIR$tab, "GSEA_Hallmark_mgroup.csv"), row.names=FALSE)

# Plot top 10 enriched
topn <- fg_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n=10) %>%
  mutate(
    pathway2 = gsub("^HALLMARK_", "", pathway),
    dir = ifelse(NES >= 0, "High", "Low")
  )

# Preserve order for plotting
topn$pathway2 <- factor(topn$pathway2, levels=rev(topn$pathway2))

p_gsea <- ggplot(topn, aes(x=pathway2, y=NES, fill=dir)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(title="GSEA Hallmark (High vs Low)", x="", y="NES") +
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust=0.5))

ggsave(file.path(DIR$fig, "PANEL_gsea_colored.png"),
       p_gsea, width=7, height=5.2, dpi=300)

message("[04] DONE: DESeq2 + PCA/Volcano/GSEA saved.")
