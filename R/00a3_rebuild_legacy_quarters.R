# ============================================================
# 00a3_rebuild_legacy_quarters.R
# Rebuild the RDS cache for the legacy AERS quarters
# 2004Q1-2012Q3 (run once, only if needed)
#
# Why this script exists: legacy quarters built by naive parsers are
# often corrupted (the first data row swallowed as the header, or
# records truncated). This script re-parses the original ASCII files
# with the correct 23/24-column AERS schema and overwrites the bad
# RDS files; quarters from 2012Q4 onwards are untouched.
#
# Two pitfalls handled here:
#   1. The legacy DRUG table comes in 12- and 13-column variants
#      (later years added DOSE_FREQ)
#   2. Legacy ASCII files contain embedded newlines inside free-text
#      fields, so one logical record can span several physical
#      lines. fread() alone truncates and misaligns such files.
#      Records are therefore stitched back together first (a new
#      record starts with "digits$"), then parsed.
#
# Legacy AERS demo layout (FDA ASCII format documentation):
#   ISR$CASE$I_F_COD$FOLL_SEQ$IMAGE$EVENT_DT$MFR_DT$FDA_DT$
#   REPT_COD$MFR_NUM$MFR_SNDR$AGE$AGE_COD$GNDR_COD$E_SUB$WT$
#   WT_COD$REPT_DT$OCCP_COD$DEATH_DT$TO_MFR$CONFID$
#   REPORTER_COUNTRY[$OCCR_COUNTRY]
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

ASCII_DIR <- Sys.getenv("FAERS_ASCII_DIR", unset = "data/faers_ascii")
OUT_DIR   <- Sys.getenv("FAERS_DIR",
                        unset = file.path(ASCII_DIR, "rds_cache_full"))

LEGACY_RENAME <- c(isr = "primaryid", `case` = "caseid", gndr_cod = "sex")

TABLE_PATTERNS <- c(demo = "^demo", drug = "^drug", reac = "^reac",
                    outc = "^outc", rpsr = "^rpsr", ther = "^ther",
                    indi = "^indi")

TRIM_COLS <- c("drugname", "prod_ai", "pt", "role_cod", "outc_cod",
               "indi_pt", "sex", "age_cod", "primaryid", "caseid")

# Legacy schemas, indexed by column count
LEGACY_SCHEMAS <- list(
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
             "dose_vbm","dechal","rechal","lot_no","exp_dt","nda_num"),
    `13` = c("isr","drug_seq","role_cod","drugname","val_vbm","route",
             "dose_vbm","dechal","rechal","lot_no","exp_dt","nda_num",
             "dose_freq")
  ),
  reac = list(`3` = c("isr","pt","drug_rec_act")),
  outc = list(`2` = c("isr","outc_cod")),
  rpsr = list(`2` = c("isr","rpsr_cod")),
  ther = list(`6` = c("isr","drug_seq","start_dt","end_dt","dur","dur_cod")),
  indi = list(`3` = c("isr","drug_seq","indi_pt"))
)

# ------------------------------------------------------------
# Stitch-then-parse reader
# The earliest quarters use old Mac-style \r line endings, which
# readLines() does not recognise; read the raw bytes and split on
# \r\n | \r | \n uniformly.
# ------------------------------------------------------------
read_legacy_table <- function(path, table_name, qtag = "") {
  raw   <- readChar(path, file.size(path), useBytes = TRUE)
  lines <- unlist(strsplit(raw, "\r\n|\r|\n", useBytes = TRUE))
  lines <- iconv(lines, from = "latin1", to = "UTF-8", sub = " ")
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) return(NULL)

  # Header detection: the first field of the first line is "ISR"
  # (any case)
  first_field <- tolower(trimws(strsplit(lines[1], "\\$")[[1]][1]))
  has_header <- first_field == "isr"
  if (has_header) {
    header_names <- tolower(trimws(strsplit(lines[1], "\\$")[[1]]))
    lines <- lines[-1]
  }

  # Stitch multi-line records: a new record starts with "digits$";
  # anything else is a continuation of the previous record
  is_start <- grepl("^[0-9]+\\$", lines)
  if (!is_start[1]) {
    stop(qtag, " ", basename(path),
         ": the first data line does not look like a record start; please inspect this file",
         call. = FALSE)
  }
  rec_id   <- cumsum(is_start)
  stitched <- unname(tapply(lines, rec_id, paste, collapse = " "))

  # Field-count consistency check: keep the modal count, drop the
  # rest. Do NOT strip trailing $ here: in early quarters the last
  # fields (death_dt, confid, ...) are often empty, and trailing
  # "$$$$" are legitimate empty fields - stripping them would break
  # the column count.
  nfields <- lengths(strsplit(stitched, "$", fixed = TRUE))
  mode_nf <- as.integer(names(sort(table(nfields), decreasing = TRUE))[1])
  n_bad   <- sum(nfields != mode_nf)
  if (n_bad > 0) {
    message("  Note: ", qtag, " ", basename(path), ": dropped ", n_bad,
            " records with anomalous field counts (",
            round(100 * n_bad / length(nfields), 2), "%)")
    stitched <- stitched[nfields == mode_nf]
  }

  dt <- fread(text = paste(stitched, collapse = "\n"),
              sep = "$", header = FALSE, quote = "",
              colClasses = "character", encoding = "UTF-8",
              fill = TRUE, showProgress = FALSE)

  if (has_header) {
    # Naming priority: the official schema wins, because some
    # quarters ship with header rows that are themselves missing the
    # last column name(s)
    sch <- LEGACY_SCHEMAS[[table_name]][[as.character(ncol(dt))]]
    if (length(header_names) == ncol(dt)) {
      names(dt) <- header_names
      message("  ", qtag, " ", basename(path), ": header row found, used as is")
    } else if (!is.null(sch)) {
      setnames(dt, sch)
      message("  Note: ", qtag, " ", basename(path),
              ": header column count (", length(header_names),
              ") does not match data (", ncol(dt),
              "); named with the ", ncol(dt), "-column official schema")
    } else {
      message("  Note: ", qtag, " ", basename(path),
              ": header column count (", length(header_names),
              ") does not match data (", ncol(dt),
              ") and no schema matches; header padded/truncated")
      hn <- header_names
      if (length(hn) < ncol(dt)) {
        hn <- c(hn, paste0("extra_", seq_len(ncol(dt) - length(hn))))
      }
      names(dt) <- hn[seq_len(ncol(dt))]
    }
  } else {
    sch <- LEGACY_SCHEMAS[[table_name]][[as.character(ncol(dt))]]
    if (is.null(sch)) {
      stop(qtag, " ", basename(path), ": no header and column count (",
           ncol(dt), ") matches no legacy schema; please inspect this file",
           call. = FALSE)
    }
    setnames(dt, sch)
    message("  ", qtag, " ", basename(path),
            ": no header row; named with the ", ncol(dt), "-column legacy schema")
  }

  if (anyDuplicated(names(dt))) {
    names(dt) <- make.unique(names(dt), sep = "_dup")
  }
  for (old in names(LEGACY_RENAME)) {
    if (old %in% names(dt) && !LEGACY_RENAME[[old]] %in% names(dt)) {
      names(dt)[names(dt) == old] <- LEGACY_RENAME[[old]]
    }
  }
  tc <- intersect(TRIM_COLS, names(dt))
  if (length(tc) > 0) dt[, (tc) := lapply(.SD, trimws), .SDcols = tc]
  as_tibble(dt)
}

process_legacy_quarter <- function(qdir, qtag) {
  tabs <- lapply(names(TABLE_PATTERNS), function(tb) {
    f <- list.files(qdir, pattern = paste0(TABLE_PATTERNS[[tb]], ".*\\.txt$"),
                    ignore.case = TRUE, recursive = TRUE, full.names = TRUE)
    if (length(f) == 0) return(NULL)
    if (length(f) > 1) f <- f[1]
    read_legacy_table(f[1], tb, qtag)
  })
  names(tabs) <- names(TABLE_PATTERNS)

  # Legacy ISR is unique only within a quarter -> synthesise a
  # globally unique primaryid
  tabs <- lapply(tabs, function(dt) {
    if (is.null(dt) || !"primaryid" %in% names(dt)) return(dt)
    dt %>% mutate(primaryid = paste0(toupper(qtag), "_", primaryid))
  })

  saveRDS(tabs, file.path(OUT_DIR, paste0(tolower(qtag), ".rds")))
  message("  Done ", qtag, ": demo=",
          if (is.null(tabs$demo)) 0L else nrow(tabs$demo),
          " drug=", if (is.null(tabs$drug)) 0L else nrow(tabs$drug),
          " reac=", if (is.null(tabs$reac)) 0L else nrow(tabs$reac))
}

# ------------------------------------------------------------
# Main loop: 2004Q1-2012Q3 only
# ------------------------------------------------------------
all_dirs <- list.dirs(ASCII_DIR, recursive = FALSE, full.names = TRUE)
qtag_raw <- str_extract(basename(all_dirs), regex("(\\d{4})q([1-4])",
                                                  ignore_case = TRUE))
keep  <- !is.na(qtag_raw)
qdirs <- all_dirs[keep]
qtags <- tolower(qtag_raw[keep])

years <- as.integer(substr(qtags, 1, 4))
qtrs  <- as.integer(substr(qtags, 6, 6))
is_legacy <- years < 2012 | (years == 2012 & qtrs <= 3)

qdirs <- qdirs[is_legacy]
qtags <- qtags[is_legacy]

message("Found ", length(qdirs), " legacy quarterly folders to rebuild")
if (length(qdirs) == 0) {
  stop("No legacy quarterly folders (e.g. aers_ascii_2004q1) found. ",
       "If the original ASCII files have been deleted, re-download the ",
       "2004Q1-2012Q3 quarterly packs from the FDA website first")
}

failed <- character(0)
for (i in seq_along(qdirs)) {
  message("Processing ", qtags[i], " (", i, "/", length(qdirs), ") ...")
  ok <- tryCatch({
    process_legacy_quarter(qdirs[i], qtags[i]); TRUE
  }, error = function(e) {
    message("  FAILED ", qtags[i], ": ", conditionMessage(e)); FALSE
  })
  if (!ok) failed <- c(failed, qtags[i])
}

if (length(failed) > 0) {
  message("\nThese quarters failed to rebuild: ", paste(failed, collapse = ", "))
} else {
  message("\nAll legacy quarters rebuilt successfully")
}

# ------------------------------------------------------------
# Post-rebuild sanity check: inspect 2006q4
# ------------------------------------------------------------
x <- readRDS(file.path(OUT_DIR, "2006q4.rds"))
message("\n=== Sanity check: 2006q4 ===")
message("demo columns: ", paste(names(x$demo), collapse = ", "))
message("fda_dt examples: ", paste(head(x$demo$fda_dt, 5), collapse = ", "))
n_rani <- x$drug %>%
  filter(role_cod == "PS",
         grepl("RANIBIZUMAB|LUCENTIS", toupper(drugname))) %>%
  distinct(primaryid) %>% nrow()
message("Ranibizumab PS reports in 2006q4: ", n_rani,
        "  (tens to hundreds expected; if still 0, please investigate)")

message("\n=== Done. If the sanity check looks right, rerun 05_sensitivity_time_aligned.R ===")
