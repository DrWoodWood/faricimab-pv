# ============================================================
# 00a_build_rds_from_ascii.R
# Build the per-quarter RDS cache from FAERS/AERS ASCII files
# (run once)
#
# Input : quarterly folders under ASCII_DIR, e.g. aers_ascii_2004q1 /
#         faers_ascii_2025q1 (each folder may contain a nested
#         ascii/ subfolder)
# Output: one .rds per quarter under OUT_DIR, each a list of the
#         7 standard tables: demo / drug / reac / outc / rpsr / ther / indi
#
# Format issues handled:
#   1. Legacy AERS (2004Q1-2012Q3) has different fields ->
#      renamed: ISR->primaryid, CASE->caseid, GNDR_COD->sex
#   2. Legacy ISR is unique only within a quarter -> a globally
#      unique primaryid is synthesised ("2004Q1_1234567") to avoid
#      ID clashes when quarters are merged
#   3. All columns are read as character; key columns are trimmed
#   4. A quarter that fails to parse does not abort the run;
#      failures are reported at the end
#
# Hardware note: all quarters in one run take ~30-60 min;
# >= 16 GB RAM recommended (quarters are processed and released
# one at a time).
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# Set the environment variable FAERS_ASCII_DIR to the folder that
# contains the quarterly ASCII folders (see README).
ASCII_DIR <- Sys.getenv("FAERS_ASCII_DIR", unset = "data/faers_ascii")
OUT_DIR   <- file.path(ASCII_DIR, "rds_cache_full")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Legacy (2004Q1-2012Q3) -> modern column name mapping
LEGACY_RENAME <- c(isr = "primaryid", `case` = "caseid", gndr_cod = "sex")

TABLE_PATTERNS <- c(demo = "^demo", drug = "^drug", reac = "^reac",
                    outc = "^outc", rpsr = "^rpsr", ther = "^ther",
                    indi = "^indi")

# Columns to trim (drug/PT matching is whitespace-sensitive)
TRIM_COLS <- c("drugname", "prod_ai", "pt", "role_cod", "outc_cod",
               "indi_pt", "sex", "age_cod", "primaryid", "caseid")

# ------------------------------------------------------------
# 0. Official standard schemas (applied by column count for
#    headerless files)
#    Source: FDA AERS/FAERS ASCII file format documentation
# ------------------------------------------------------------
SCHEMAS <- list(
  legacy = list(   # 2004Q1-2012Q3 (AERS)
    demo = list(
      `23` = c("isr","case","i_f_cod","foll_seq","image","event_dt","mfr_dt",
               "fda_dt","rept_cod","mfr_num","mfr_sndr","age","age_cod",
               "gndr_cod","e_sub","wt","wt_cod","rept_dt","occp_cod",
               "death_dt","to_mfr","confid","reporter_country"),
      `24` = c("isr","case","i_f_cod","foll_seq","image","event_dt","mfr_dt",
               "fda_dt","rept_cod","mfr_num","mfr_sndr","age","age_cod",
               "gndr_cod","e_sub","wt","wt_cod","rept_dt","occp_cod",
               "death_dt","to_mfr","confid","reporter_country","occr_country")
    ),
    drug = list(
      `12` = c("isr","drug_seq","role_cod","drugname","val_vbm","route",
               "dose_vbm","dechal","rechal","lot_no","exp_dt","nda_num")
    ),
    reac = list(`3` = c("isr","pt","drug_rec_act")),
    outc = list(`2` = c("isr","outc_cod")),
    rpsr = list(`2` = c("isr","rpsr_cod")),
    ther = list(`6` = c("isr","drug_seq","start_dt","end_dt","dur","dur_cod")),
    indi = list(`3` = c("isr","drug_seq","indi_pt"))
  ),
  modern = list(   # 2012Q4 onwards (FAERS)
    demo = list(
      `24` = c("primaryid","caseid","caseversion","i_f_code","event_dt",
               "mfr_dt","init_fda_dt","fda_dt","rept_cod","auth_num",
               "mfr_num","mfr_sndr","lit_ref","age","age_cod","sex",
               "e_sub","wt","wt_cod","rept_dt","to_mfr","occp_cod",
               "reporter_country","occr_country"),
      `25` = c("primaryid","caseid","caseversion","i_f_code","event_dt",
               "mfr_dt","init_fda_dt","fda_dt","rept_cod","auth_num",
               "mfr_num","mfr_sndr","lit_ref","age","age_cod","age_grp",
               "sex","e_sub","wt","wt_cod","rept_dt","to_mfr","occp_cod",
               "reporter_country","occr_country")
    ),
    drug = list(
      `16` = c("primaryid","drug_seq","role_cod","drugname","prod_ai",
               "val_vbm","route","dose_vbm","cum_dose_chr","cum_dose_unit",
               "dechal","rechal","lot_num","exp_dt","nda_num","dose_freq")
    ),
    reac = list(`3` = c("primaryid","pt","drug_rec_act")),
    outc = list(`2` = c("primaryid","outc_cod")),
    rpsr = list(`2` = c("primaryid","rpsr_cod")),
    ther = list(`6` = c("primaryid","dsg_drug_seq","start_dt","end_dt",
                        "dur","dur_cod")),
    indi = list(`3` = c("primaryid","indi_drug_seq","indi_pt"))
  )
)

# ------------------------------------------------------------
# 1. Read a single ASCII table ($ separated, latin1, all character)
#    Tolerance: some quarterly files lack a header row (e.g.
#    2004Q2); fread would mistake the first data row for a header,
#    so the first line is checked for a valid header first,
#    otherwise the official schema is applied.
# ------------------------------------------------------------
read_ascii_table <- function(path, table_name, legacy, qtag = "") {
  first_line <- readLines(path, n = 1, warn = FALSE)
  first_field <- tolower(trimws(strsplit(first_line, "\\$")[[1]][1]))
  has_header <- grepl("\\$", first_line) &&
    first_field %in% c("isr", "primaryid")

  dt <- tryCatch(
    fread(path, sep = "$", header = has_header, quote = "",
          colClasses = "character", encoding = "Latin-1",
          showProgress = FALSE),
    error = function(e) stop("fread failed on ", basename(path), ": ",
                             e$message, call. = FALSE)
  )

  if (has_header) {
    names(dt) <- tolower(trimws(names(dt)))
  } else {
    sch <- SCHEMAS[[if (legacy) "legacy" else "modern"]][[table_name]][[as.character(ncol(dt))]]
    if (is.null(sch)) {
      stop(qtag, " ", basename(path), ": no header and column count (", ncol(dt),
           ") matches no standard schema; please inspect this file manually",
           call. = FALSE)
    }
    setnames(dt, sch)
    message("  Note: ", qtag, " ", basename(path),
            " has no header row; named with the ", ncol(dt), "-column standard schema")
  }

  # Guard against duplicated column names (tibble forbids them)
  if (anyDuplicated(names(dt))) {
    names(dt) <- make.unique(names(dt), sep = "_dup")
    message("  Note: ", qtag, " ", basename(path),
            " has duplicated column names; made unique")
  }

  # Legacy column name mapping
  for (old in names(LEGACY_RENAME)) {
    if (old %in% names(dt) && !LEGACY_RENAME[[old]] %in% names(dt)) {
      names(dt)[names(dt) == old] <- LEGACY_RENAME[[old]]
    }
  }
  tc <- intersect(TRIM_COLS, names(dt))
  if (length(tc) > 0) dt[, (tc) := lapply(.SD, trimws), .SDcols = tc]
  as_tibble(dt)
}

# ------------------------------------------------------------
# 2. Process one quarterly folder
# ------------------------------------------------------------
process_quarter <- function(qdir, qtag) {
  year <- as.integer(substr(qtag, 1, 4))
  qtr  <- as.integer(substr(qtag, 6, 6))
  legacy <- (year < 2012) || (year == 2012 && qtr <= 3)

  tabs <- lapply(names(TABLE_PATTERNS), function(tb) {
    p <- TABLE_PATTERNS[[tb]]
    f <- list.files(qdir, pattern = paste0(p, ".*\\.txt$"),
                    ignore.case = TRUE, recursive = TRUE,
                    full.names = TRUE)
    if (length(f) == 0) return(NULL)
    if (length(f) > 1) {
      message("  Note: ", qtag, " matched several ", tb, " files; using the first: ",
              basename(f[1]))
      f <- f[1]
    }
    read_ascii_table(f, table_name = tb, legacy = legacy, qtag = qtag)
  })
  names(tabs) <- names(TABLE_PATTERNS)

  # Synthesise a globally unique primaryid for legacy quarters
  if (legacy) {
    tabs <- lapply(tabs, function(dt) {
      if (is.null(dt) || !"primaryid" %in% names(dt)) return(dt)
      dt %>% mutate(primaryid = paste0(toupper(qtag), "_", primaryid))
    })
  } else {
    tabs <- lapply(tabs, function(dt) {
      if (is.null(dt) || !"primaryid" %in% names(dt)) return(dt)
      dt %>% mutate(primaryid = as.character(primaryid))
    })
  }

  saveRDS(tabs, file.path(OUT_DIR, paste0(tolower(qtag), ".rds")))
  tibble(quarter = qtag, legacy = legacy,
         n_demo = if (is.null(tabs$demo)) 0L else nrow(tabs$demo),
         n_drug = if (is.null(tabs$drug)) 0L else nrow(tabs$drug),
         n_reac = if (is.null(tabs$reac)) 0L else nrow(tabs$reac),
         tables = paste(names(tabs)[!vapply(tabs, is.null, logical(1))],
                        collapse = ","))
}

# ------------------------------------------------------------
# 3. Scan quarterly folders and process one by one
# ------------------------------------------------------------
all_dirs <- list.dirs(ASCII_DIR, recursive = FALSE, full.names = TRUE)
qtag <- str_extract(basename(all_dirs), regex("(\\d{4})q([1-4])",
                                              ignore_case = TRUE))
qdirs <- all_dirs[!is.na(qtag)]
qtags <- tolower(qtag[!is.na(qtag)])

if (length(qdirs) == 0) {
  stop("No quarterly folders like aers_ascii_2004q1 found under ", ASCII_DIR)
}

ord <- order(substr(qtags, 1, 4), substr(qtags, 6, 6))
qdirs <- qdirs[ord]; qtags <- qtags[ord]

message("Found ", length(qdirs), " quarterly folders: ",
        qtags[1], " ... ", qtags[length(qtags)])

results <- vector("list", length(qdirs))
failed  <- character(0)

for (i in seq_along(qdirs)) {
  out_f <- file.path(OUT_DIR, paste0(qtags[i], ".rds"))
  if (file.exists(out_f)) {
    message(sprintf("[%d/%d] %s -- already exists, skipped", i, length(qdirs), qtags[i]))
    next   # only missing/deleted quarters are rebuilt; delete the
           # .rds first to force a rebuild
  }
  message(sprintf("[%d/%d] %s", i, length(qdirs), qtags[i]))
  results[[i]] <- tryCatch(
    process_quarter(qdirs[i], qtags[i]),
    error = function(e) {
      failed <<- c(failed, qtags[i])
      message("  !! Parse failed: ", qtags[i], " -- ", e$message)
      NULL
    }
  )
  gc(verbose = FALSE)
}

summary_df <- bind_rows(results)
print(summary_df, n = Inf)

if (length(failed) > 0) {
  message("These quarters failed to parse; check them individually: ",
          paste(failed, collapse = ", "))
}

readr::write_csv(summary_df, file.path(OUT_DIR, "_build_summary.csv"))
message("=== Done. RDS output directory: ", OUT_DIR, " ===")
message("Next step: point FAERS_DIR in scripts 01/06 at this directory")
