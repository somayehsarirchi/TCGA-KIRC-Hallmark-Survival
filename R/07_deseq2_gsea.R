source(here::here("R","00_utils.R"))

step_07_deseq2_gsea <- function(cache_in_pre=here::here("cache","02_preprocess"),
                                cache_in_lasso=here::here("cache","06_lasso"),
                                cache_in_scores=here::here("cache","04_scores"),
                                cache_out=here::here("cache","07_deg_gsea"),
                                fig_dir=here::here("results","figures"),
                                out_tables=here::here("results","tables"),
                                log_file=NULL,
                                force=FALSE) {
  
  suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
    library(ggplot2)
    library(fgsea)
  })
  
  dir_create(cache_out); dir_create(fig_dir); dir_create(out_tables)
  
  in_counts2 <- file.path(cache_in_pre, "counts2.rds")
  in_meta_risk <- file.path(cache_in_lasso, "meta2_risk.rds")
  in_gene_annot <- file.path(cache_in_pre, "gene_annot.rds")
  in_pathways <- file.path(cache_in_scores, "pathways_hallmark.rds")
  
  out_dds <- file.path(cache_out, "dds.rds")
  out_res <- file.path(cache_out, "res_annot.rds")
  out_sig <- file.path(cache_out, "sig_annot.rds")
  out_fg  <- file.path(cache_out, "fg_df.rds")
  
  out_res_csv <- file.path(out_tables, "DESeq2_mgroup_adjusted_fullResults.csv")
  out_sig_csv <- file.path(out_tables, "DESeq2_mgroup_significant_genes.csv")
  out_gsea_csv <- file.path(out_tables, "GSEA_Hallmark_mgroup.csv")
  
  pca_png <- file.path(fig_dir, "PANEL_pca.png")
  volc_png <- file.path(fig_dir, "PANEL_volcano_labeled.png")
  gsea_png <- file.path(fig_dir, "PANEL_gsea_colored.png")
  
  run_step("07_deseq2_gsea",
           c(out_dds, out_res, out_sig, out_fg, out_res_csv, out_sig_csv, out_gsea_csv,
             pca_png, volc_png, gsea_png),
           force=force, log_file=log_file,
           fn_build=function() {
             
             counts2 <- read_rds(in_counts2, log_file)
             meta2   <- read_rds(in_meta_risk, log_file)
             gene_annot <- read_rds(in_gene_annot, log_file)
             pathways <- read_rds(in_pathways, log_file)
             
             dds <- DESeqDataSetFromMatrix(
               countData=round(counts2),
               colData=meta2,
               design=~ age_years + stage_bin + mgroup
             )
             dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]
             dds <- DESeq(dds)
             
             res <- results(dds, contrast=c("mgroup","High","Low"))
             res <- res[order(res$padj), ]
             res_df <- as.data.frame(res); res_df$gene_id <- rownames(res_df)
             
             res_annot <- res_df %>%
               left_join(gene_annot %>% dplyr::select(gene_id, gene_name), by="gene_id") %>%
               arrange(padj)
             
             sig_annot <- subset(res_annot, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
             
             write.csv(res_annot, out_res_csv, row.names=FALSE)
             write.csv(sig_annot, out_sig_csv, row.names=FALSE)
             
             # PCA
             vsd <- vst(dds, blind=FALSE)
             pcaData <- plotPCA(vsd, intgroup="mgroup", returnData=TRUE)
             pal_group <- c("Low"="#2C7BB6","High"="#D7191C")
             p_pca <- ggplot(pcaData, aes(PC1, PC2, color=mgroup)) +
               geom_point(size=3, alpha=0.85) +
               scale_color_manual(values=pal_group) +
               theme_classic() +
               ggtitle("PCA - Molecular group (High vs Low)") +
               theme(plot.title=element_text(hjust=0.5))
             ggsave(pca_png, p_pca, width=7, height=5, dpi=300)
             
             # Volcano (requires ggrepel)
             if (!requireNamespace("ggrepel", quietly=TRUE)) stopf("Install ggrepel", log_file=log_file)
             vol <- res_annot
             vol$neglog10padj <- -log10(vol$padj)
             vol$cls <- "NotSig"
             vol$cls[!is.na(vol$padj) & vol$padj < 0.05 & vol$log2FoldChange >  1] <- "Up"
             vol$cls[!is.na(vol$padj) & vol$padj < 0.05 & vol$log2FoldChange < -1] <- "Down"
             
             top_up <- vol %>% filter(cls=="Up", !is.na(gene_name), gene_name!="") %>% arrange(padj) %>% slice_head(n=5)
             top_dn <- vol %>% filter(cls=="Down", !is.na(gene_name), gene_name!="") %>% arrange(padj) %>% slice_head(n=5)
             lab_df <- bind_rows(top_up, top_dn)
             
             p_volc <- ggplot(vol, aes(log2FoldChange, neglog10padj, color=cls)) +
               geom_point(alpha=0.6, size=1.1) +
               ggrepel::geom_text_repel(data=lab_df, aes(label=gene_name), size=4,
                                        box.padding=0.35, point.padding=0.25, max.overlaps=Inf) +
               theme_classic() + labs(x="log2FC (High - Low)", y="-log10(padj)") +
               scale_color_manual(values=c(Down="blue", NotSig="grey70", Up="red"))
             ggsave(volc_png, p_volc, width=7, height=5.2, dpi=300)
             
             # GSEA ranks
             tmp <- vol[, c("gene_name","stat")]
             tmp <- tmp[!is.na(tmp$gene_name) & tmp$gene_name!="" & !is.na(tmp$stat), ]
             tmp <- tmp[order(abs(tmp$stat), decreasing=TRUE), ]
             tmp <- tmp[!duplicated(tmp$gene_name), ]
             ranks <- tmp$stat; names(ranks) <- tmp$gene_name
             ranks <- sort(ranks, decreasing=TRUE)
             
             fg <- fgseaMultilevel(pathways=pathways, stats=ranks, minSize=15, maxSize=500)
             fg <- fg[order(fg$padj), ]
             fg_df <- as.data.frame(fg)
             if ("leadingEdge" %in% names(fg_df)) fg_df$leadingEdge <- sapply(fg_df$leadingEdge, paste, collapse=";")
             write.csv(fg_df, out_gsea_csv, row.names=FALSE)
             
             topn <- fg_df %>%
               arrange(padj) %>% slice_head(n=10) %>%
               mutate(pathway2=gsub("^HALLMARK_","", pathway),
                      dir=ifelse(NES >= 0, "High","Low"))
             topn$pathway2 <- factor(topn$pathway2, levels=rev(topn$pathway2))
             
             p_gsea <- ggplot(topn, aes(pathway2, NES, fill=dir)) +
               geom_col() + coord_flip() + theme_classic() +
               labs(x="", y="NES") +
               scale_fill_manual(values=c("Low"="#2C7BB6","High"="#D7191C")) +
               theme(legend.title=element_blank())
             ggsave(gsea_png, p_gsea, width=7, height=5.2, dpi=300)
             
             save_rds(dds, out_dds, log_file)
             save_rds(res_annot, out_res, log_file)
             save_rds(sig_annot, out_sig, log_file)
             save_rds(fg_df, out_fg, log_file)
           })
  
  invisible(TRUE)
}