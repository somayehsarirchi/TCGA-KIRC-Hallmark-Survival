suppressPackageStartupMessages({
  library(here)
  library(glue)
})

dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

log_msg <- function(..., level="INFO", log_file=NULL) {
  msg <- glue("[{format(Sys.time(), '%Y-%m-%d %H:%M:%S')}] [{level}] ", paste0(..., collapse=""))
  message(msg)
  if (!is.null(log_file)) cat(msg, "\n", file=log_file, append=TRUE)
}

stopf <- function(..., log_file=NULL) {
  log_msg(..., level="ERROR", log_file=log_file)
  stop(paste0(..., collapse=""), call.=FALSE)
}

save_rds <- function(obj, path, log_file=NULL) {
  dir_create(dirname(path))
  saveRDS(obj, path)
  log_msg("Saved: ", path, log_file=log_file)
}

read_rds <- function(path, log_file=NULL) {
  if (!file.exists(path)) stopf("Missing RDS: ", path, log_file=log_file)
  log_msg("Loaded: ", path, log_file=log_file)
  readRDS(path)
}

assert_nonempty <- function(x, what="object", log_file=NULL) {
  if (is.null(x)) stopf(what, " is NULL", log_file=log_file)
  if (is.data.frame(x) && nrow(x) == 0) stopf(what, " has 0 rows", log_file=log_file)
  if (is.matrix(x) && (nrow(x) == 0 || ncol(x) == 0)) stopf(what, " has 0 dims", log_file=log_file)
  if (is.vector(x) && length(x) == 0) stopf(what, " length 0", log_file=log_file)
  invisible(TRUE)
}

run_step <- function(name, out_files, fn_build, force=FALSE, log_file=NULL) {
  log_msg("---- STEP: ", name, " ----", log_file=log_file)
  ok <- all(file.exists(out_files))
  if (ok && !force) {
    log_msg("Cache hit. Skipping build.", log_file=log_file)
    return(invisible(TRUE))
  }
  log_msg("Cache miss (or force=TRUE). Building...", log_file=log_file)
  fn_build()
  missing <- out_files[!file.exists(out_files)]
  if (length(missing) > 0) stopf("Step '", name, "' did not produce:\n", paste(missing, collapse="\n"), log_file=log_file)
  log_msg("Step completed: ", name, log_file=log_file)
  invisible(TRUE)
}

# --- helpers for robust column detection ---
pick_first_col <- function(df, candidates) {
  candidates <- candidates[candidates %in% colnames(df)]
  if (length(candidates) == 0) return(NA_character_)
  candidates[1]
}

detect_grade_col <- function(meta) {
  # try common TCGA names
  cand <- c("tumor_grade", "grade", "ajcc_tumor_grade", "neoplasm_histologic_grade", "tumor_grade_anatomic_site")
  x <- pick_first_col(meta, cand)
  if (!is.na(x)) return(x)
  # fallback: grep "grade"
  g <- grep("grade", colnames(meta), ignore.case=TRUE, value=TRUE)
  if (length(g) > 0) return(g[1])
  NA_character_
}

make_grade_bin <- function(meta, log_file=NULL) {
  gc <- detect_grade_col(meta)
  if (is.na(gc)) {
    log_msg("WARN: grade column not found. grade_bin will be NA.", level="WARN", log_file=log_file)
    meta$grade_bin <- factor(NA_character_, levels=c("Low","High"))
    return(meta)
  }
  x <- toupper(trimws(as.character(meta[[gc]])))
  meta$tumor_grade <- x
  gb <- dplyr::case_when(
    x %in% c("G1","G2") ~ "Low",
    x %in% c("G3","G4") ~ "High",
    TRUE ~ NA_character_
  )
  meta$grade_bin <- factor(gb, levels=c("Low","High"))
  meta
}