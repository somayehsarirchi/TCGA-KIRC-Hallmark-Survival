# =====================================================
# 00_config.R
# Global config (paths + common helpers)
# =====================================================

suppressPackageStartupMessages({
  library(here)
})

set.seed(1)

DIR <- list(
  root   = here::here(),
  cache  = here::here("data_cache"),
  res    = here::here("results"),
  fig    = here::here("results", "figures"),
  tab    = here::here("results", "tables")
)

dir.create(DIR$cache, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR$res,   showWarnings = FALSE, recursive = TRUE)
dir.create(DIR$fig,   showWarnings = FALSE, recursive = TRUE)
dir.create(DIR$tab,   showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, filename, w=7, h=5, dpi=300){
  ggplot2::ggsave(file.path(DIR$fig, filename), plot=p, width=w, height=h, dpi=dpi)
}

message("Project root: ", DIR$root)
