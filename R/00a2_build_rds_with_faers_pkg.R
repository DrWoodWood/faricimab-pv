# ============================================================
# 00a2_build_rds_with_faers_pkg.R
# Download and parse quarterly FAERS data with the faers
# (Bioconductor) package, producing the same per-quarter RDS
# structure as 00a (demo/drug/reac/outc/rpsr/ther/indi).
#
# Difference from 00a:
#   - 00a : reads locally unzipped ASCII folders, no download
#   - 00a2: calls the faers package, which downloads the zips from
#           the FDA website and parses them (robust; the package
#           handles the legacy/modern format differences and column
#           mapping internally)
#
# Notes:
#   1. If you still have the original zip files, put them in
#      ZIP_DIR; faers will reuse them and skip the download
#      (the console shows "already downloaded")
#   2. Quarters are processed, saved and released one at a time;
#      ~90 quarters take about 1-3 hours depending on bandwidth,
#      with several GB of disk usage
#   3. This script does not de-duplicate or standardise;
#      de-duplication is done downstream by load_faers_quarters()
# ============================================================

suppressPackageStartupMessages({
  library(faers)
  library(dplyr)
})

# Set the environment variables FAERS_ZIP_DIR / FAERS_DIR to point
# elsewhere (see README).
ZIP_DIR <- Sys.getenv("FAERS_ZIP_DIR", unset = "data/zip_cache")
OUT_DIR <- Sys.getenv("FAERS_DIR", unset = "data/rds_cache_full")
dir.create(ZIP_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# faers-package fields -> this project's standard schema
# (column names returned by faers_get follow the package docs; the
# list below reflects the vignette's observed names. Missing columns
# are filled with NA without error.)
FIELD_COLS <- list(
  demo = c("primaryid", "caseid", "caseversion", "fda_dt", "event_dt",
           "age", "age_cod", "sex", "occp_cod", "occr_country",
           "reporter_country"),
  drug = c("primaryid", "drug_seq", "role_cod", "drugname", "prod_ai",
           "route"),
  reac = c("primaryid", "pt"),
  outc = c("primaryid", "outc_cod"),
  rpsr = c("primaryid", "rpsr_cod"),
  ther = c("primaryid", "drug_seq", "start_dt", "end_dt"),
  indi = c("primaryid", "drug_seq", "indi_pt")
)

# ------------------------------------------------------------
# Per quarter: download + parse -> extract the 7 tables ->
# harmonise schema -> save RDS
# ------------------------------------------------------------
process_one_quarter <- function(year, qtr) {
  qtag <- paste0(year, qtr)                      # e.g. "2004q1"
  out_file <- file.path(OUT_DIR, paste0(qtag, ".rds"))
  if (file.exists(out_file)) {
    message("Already exists, skipped: ", qtag)
    return(invisible(TRUE))
  }
  message("Processing: ", qtag)

  obj <- faers(year, qtr, dir = ZIP_DIR)         # download + parse one quarter

  tabs <- lapply(names(FIELD_COLS), function(tb) {
    dt <- tryCatch(faers_get(obj, tb), error = function(e) NULL)
    if (is.null(dt)) return(NULL)
    dt <- as_tibble(dt)
    # keep only the standard columns; fill missing ones with NA
    for (cc in setdiff(FIELD_COLS[[tb]], names(dt))) dt[[cc]] <- NA_character_
    dt %>% select(all_of(FIELD_COLS[[tb]])) %>%
      mutate(across(everything(), as.character))
  })
  names(tabs) <- names(FIELD_COLS)

  saveRDS(tabs, out_file)
  rm(obj, tabs); gc(verbose = FALSE)
  invisible(TRUE)
}

# All quarters from 2004Q1 to 2026Q2
grid <- expand.grid(year = 2004:2026, qtr = paste0("q", 1:4),
                    stringsAsFactors = FALSE) %>%
  mutate(qtag = paste0(year, qtr)) %>%
  filter(qtag <= "2026q2") %>%
  arrange(year, qtr)

failed <- character(0)
for (i in seq_len(nrow(grid))) {
  ok <- tryCatch(
    process_one_quarter(grid$year[i], grid$qtr[i]),
    error = function(e) {
      failed <<- c(failed, paste0(grid$year[i], grid$qtr[i]))
      message("  !! Failed: ", grid$year[i], grid$qtr[i], " -- ", e$message)
      FALSE
    }
  )
}

if (length(failed) > 0) {
  message("These quarters failed; simply rerun this script to resume ",
          "(existing quarters are skipped automatically): ",
          paste(failed, collapse = ", "))
} else {
  message("=== All quarters done. Output directory: ", OUT_DIR, " ===")
}
message("Next step: source('R/01_faers_prepare.R')")
