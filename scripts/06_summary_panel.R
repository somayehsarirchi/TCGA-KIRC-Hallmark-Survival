# =====================================================
# 06_summary_panel.R
# SUMMARY PANEL (without workflow) - Clean, compact, labeled
# Output:
#   results/FIGURE_SummaryPanel_Molecular.png
# =====================================================

source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(grid)
  library(gridExtra)
  library(png)
})

# ------------------------------
# Helpers
# ------------------------------
read_png_grob_safe <- function(path, label = NULL) {
  if (!file.exists(path)) {
    # placeholder if file missing
    msg <- paste0("Missing:\n", basename(path))
    if (!is.null(label)) msg <- paste0(label, "  ", msg)
    return(
      grobTree(
        rectGrob(gp = gpar(fill = "grey95", col = "grey70")),
        textGrob(msg, gp = gpar(col = "grey30", fontsize = 12))
      )
    )
  }
  rasterGrob(readPNG(path), interpolate = TRUE)
}

panel_with_label <- function(img_grob, label) {
  # label in top-left corner
  grobTree(
    img_grob,
    textGrob(
      label, x = unit(0.01, "npc"), y = unit(0.99, "npc"),
      just = c("left", "top"),
      gp = gpar(fontface = "bold", fontsize = 16, col = "black")
    )
  )
}

# ------------------------------
# Inputs (from previous steps)
# ------------------------------
fig_dir <- DIR$fig  # from 00_config.R (results/figures)

pca_path     <- file.path(fig_dir, "PANEL_pca.png")
volcano_path <- file.path(fig_dir, "PANEL_volcano_labeled.png")
gsea_path    <- file.path(fig_dir, "PANEL_gsea_colored.png")
km_path      <- file.path(fig_dir, "KM_mgroup.png")
roc_path     <- file.path(fig_dir, "PANEL_timeROC.png")

g_pca     <- read_png_grob_safe(pca_path)
g_volcano <- read_png_grob_safe(volcano_path)
g_gsea    <- read_png_grob_safe(gsea_path)
g_km      <- read_png_grob_safe(km_path)
g_roc     <- read_png_grob_safe(roc_path)

# Add panel letters
g_pca     <- panel_with_label(g_pca,     "A")
g_volcano <- panel_with_label(g_volcano, "B")
g_gsea    <- panel_with_label(g_gsea,    "C")
g_km      <- panel_with_label(g_km,      "D")
g_roc     <- panel_with_label(g_roc,     "E")

# empty cell (keeps symmetry)
g_empty <- nullGrob()

# ------------------------------
# Layout: 2 columns x 3 rows
# ------------------------------
panel <- arrangeGrob(
  g_pca,     g_volcano,
  g_gsea,    g_km,
  g_roc,     g_empty,
  ncol = 2,
  heights = c(1, 1, 1),
  widths  = c(1, 1),
  padding = unit(0.6, "lines")
)

# Optional outer border/background for clean export
panel_final <- grobTree(
  rectGrob(gp = gpar(fill = "white", col = NA)),
  panel
)

# ------------------------------
# Save
# ------------------------------
out_path <- file.path(DIR$res, "FIGURE_SummaryPanel_Molecular.png")

png(out_path, width = 2400, height = 3000, res = 300)
grid.draw(panel_final)
dev.off()

message("[06] DONE: Summary panel saved -> ", out_path)
