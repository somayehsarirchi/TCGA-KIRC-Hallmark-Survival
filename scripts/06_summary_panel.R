# =====================================================
# PROFESSIONAL SUMMARY PANEL (Clean Layout)
# =====================================================

library(gridExtra)
library(png)
library(grid)

fig_dir <- file.path("results", "figures")

# Load images
g_pca      <- rasterGrob(readPNG(file.path(fig_dir, "PANEL_pca.png")))
g_volcano  <- rasterGrob(readPNG(file.path(fig_dir, "PANEL_volcano_labeled.png")))
g_gsea     <- rasterGrob(readPNG(file.path(fig_dir, "PANEL_gsea_colored.png")))
g_km       <- rasterGrob(readPNG(file.path(fig_dir, "KM_mgroup.png")))
g_roc      <- rasterGrob(readPNG(file.path(fig_dir, "PANEL_timeROC.png")))

# Arrange in clean 2-column layout
panel_clean <- grid.arrange(
  g_pca,
  g_volcano,  g_gsea,
  g_km,       g_roc,
  ncol = 2,
  heights = c(1,1,1)
)

png(
  file.path("results", "FIGURE_SummaryPanel_Molecular.png"),
  width = 2000,
  height = 2600,
  res = 300
)

grid.draw(panel_clean)
dev.off()
