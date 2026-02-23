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
read_png_grob_safe <- function(path) {
  if (!file.exists(path)) {
    msg <- paste0("Missing file:\n", basename(path))
    return(
      grobTree(
        rectGrob(gp = gpar(fill = "grey95", col = "grey70")),
        textGrob(
          msg,
          x = unit(0.5, "npc"), y = unit(0.5, "npc"),
          just = "center",
          gp = gpar(col = "grey30", fontsize = 12)
        )
      )
    )
  }
  rasterGrob(readPNG(path), interpolate = TRUE)
}

panel_with_label <- function(img_grob, label) {
  grobTree(
    img_grob,
    textGrob(
      label, x = unit(0.01, "npc"), y = unit(0.99, "npc"),
      just = c("left", "top"),
      gp = gpar(fontface = "bold", fontsize = 16, col = "black")
    )
  )
}

caption_grob <- function(text = "KIRC | Hallmark-based molecular risk model") {
  grobTree(
    rectGrob(gp = gpar(fill = "white", col = "grey90")),
    textGrob(
      text,
      x = unit(0.5, "npc"), y = unit(0.5, "npc"),
      just = "center",
      gp = gpar(col = "grey25", fontsize = 12, fontface = "plain")
    )
  )
}

# ------------------------------
# Inputs
# ------------------------------
fig_dir <- DIR$fig

paths <- list(
  A = file.path(fig_dir, "PANEL_pca.png"),
  B = file.path(fig_dir, "PANEL_volcano_labeled.png"),
  C = file.path(fig_dir, "PANEL_gsea_colored.png"),
  D = file.path(fig_dir, "KM_mgroup.png"),
  E = file.path(fig_dir, "PANEL_timeROC.png"),
  F = file.path(fig_dir, "PANEL_calibration.png")  # NEW
)

gA <- panel_with_label(read_png_grob_safe(paths$A), "A")
gB <- panel_with_label(read_png_grob_safe(paths$B), "B")
gC <- panel_with_label(read_png_grob_safe(paths$C), "C")
gD <- panel_with_label(read_png_grob_safe(paths$D), "D")
gE <- panel_with_label(read_png_grob_safe(paths$E), "E")

# Panel F: prefer calibration; if missing, show caption (or placeholder)
if (file.exists(paths$F)) {
  gF <- panel_with_label(read_png_grob_safe(paths$F), "F")
} else {
  # choose fallback behavior:
  # 1) caption (nicer than "missing file"), OR
  # 2) placeholder from read_png_grob_safe(paths$F)
  gF <- caption_grob("Calibration plot not available (run 07_rms_calibration.R)")
}

# ------------------------------
# Layout: 2 columns x 3 rows
# ------------------------------
panel <- arrangeGrob(
  gA, gB,
  gC, gD,
  gE, gF,
  ncol = 2,
  heights = c(1, 1, 1),
  widths  = c(1, 1),
  padding = unit(0.6, "lines")
)

panel_final <- grobTree(
  rectGrob(gp = gpar(fill = "white", col = NA)),
  panel
)

# ------------------------------
# Save
# ------------------------------
out_path <- file.path(DIR$res, "FIGURE_SummaryPanel_Molecular.png")

use_cairo <- capabilities("cairo")
if (use_cairo) {
  png(out_path, width = 2400, height = 3000, res = 300, type = "cairo-png")
} else {
  png(out_path, width = 2400, height = 3000, res = 300)
}

grid.draw(panel_final)
dev.off()

message("[06] DONE: Summary panel saved -> ", out_path)
