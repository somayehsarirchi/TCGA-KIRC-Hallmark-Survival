# =====================================================
# SUMMARY PANEL (without workflow) - Clean & compact
# =====================================================

library(gridExtra)
library(png)
library(grid)
library(here)

fig_dir <- here("results", "figures")

g_pca      <- rasterGrob(readPNG(file.path(fig_dir, "PANEL_pca.png")))
g_volcano  <- rasterGrob(readPNG(file.path(fig_dir, "PANEL_volcano_labeled.png")))
g_gsea     <- rasterGrob(readPNG(file.path(fig_dir, "PANEL_gsea_colored.png")))
g_km       <- rasterGrob(readPNG(file.path(fig_dir, "KM_mgroup.png")))
g_roc      <- rasterGrob(readPNG(file.path(fig_dir, "PANEL_timeROC.png")))

# 2 columns x 3 rows
panel_no_workflow <- arrangeGrob(
  g_pca,     g_volcano,
  g_gsea,    g_km,
  g_roc,     nullGrob(),   # empty cell to keep symmetry
  ncol = 2,
  heights = c(1,1,1)
)

png(
  here("results", "FIGURE_SummaryPanel_Molecular.png"),
  width = 2200,
  height = 2800,
  res = 300
)
grid.draw(panel_no_workflow)
dev.off()
