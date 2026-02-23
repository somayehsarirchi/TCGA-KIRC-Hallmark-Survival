# =====================================================
# 00_config.R
# Global config (paths + common helpers)
# =====================================================

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
})

set.seed(1)

# ------------------------------
# Paths
# ------------------------------
DIR <- list(
  root   = here::here(),
  cache  = here::here("data_cache"),
  res    = here::here("results"),
  fig    = here::here("results", "figures"),
  tab    = here::here("results", "tables")
)

assert_dir <- function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
  invisible(x)
}

# Create folders
invisible(lapply(DIR, assert_dir))

# ------------------------------
# Global palette
# ------------------------------
PAL_GROUP <- c("Low"="#2C7BB6", "High"="#D7191C")

# ------------------------------
# Plot defaults (consistent look)
# ------------------------------
theme_set(theme_classic(base_size = 13))
theme_update(
  plot.title = element_text(hjust = 0.5),
  legend.title = element_blank()
)

# ------------------------------
# Helpers
# ------------------------------
save_plot <- function(p, filename, w=7, h=5, dpi=300) {
  # If filename is an absolute/relative path containing '/', keep it.
  # Else save into DIR$fig
  out <- if (grepl("[/\\\\]", filename)) filename else file.path(DIR$fig, filename)
  assert_dir(dirname(out))
  ggplot2::ggsave(out, plot = p, width = w, height = h, dpi = dpi)
  invisible(out)
}

check_required_cols <- function(df, cols, df_name = deparse(substitute(df))) {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("[%s] Missing required columns: %s",
                 df_name, paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}

msg_header <- function(step = NULL) {
  cat("\n==============================\n")
  if (!is.null(step)) cat("STEP: ", step, "\n", sep="")
  cat("Project root: ", DIR$root, "\n", sep="")
  cat("Time: ", format(Sys.time()), "\n", sep="")
  cat("==============================\n\n")
}

message("Project root: ", DIR$root)
