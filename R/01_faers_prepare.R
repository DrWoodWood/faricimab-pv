# ============================================================
# 01_faers_prepare.R
# FAERS data preparation: drug identification, primary-suspect (PS)
# restriction, deduplication, and construction of the analysis cohorts
# (report level / PT long table / SOC long table / backgrounds / TTO)
#
# Outputs (saved to workdir/):
#   cohort_reports.rds  report-level covariates (for adjusted ROR)
#   pairs_pt.rds        three-drug PT long table (case_id, drug, event)
#   pairs_soc.rds       three-drug SOC long table
#   bg_pt.rds           full-database background PT long table (reports
#                       without any of the three drugs)
#   bg_soc.rds          full-database background SOC long table
#   tto.rds             time-to-onset (case_id, drug, pt, tto_days)
#   antivegf_any.rds    primaryids containing any anti-VEGF (any role)
#   drug_all.rds        all concomitant drugs of cohort reports (any role)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(lubridate); library(purrr)
})
source("R/00_utils.R")

# ------------------------------------------------------------
# 0. Load quarterly FAERS objects and merge.
#    load_faers_quarters() sorts .rds files by name, row-binds
#    demo/drug/reac/outc/ther/indi across quarters, and deduplicates
#    (keeps the newest version of each caseid by max primaryid).
# ------------------------------------------------------------

FAERS_DIR <- Sys.getenv("FAERS_DIR", unset = "data/rds_cache_full")

# The primary analysis only needs post-faricimab-approval quarters;
# filtering at load avoids merging everything back to 2004.
faers_obj <- load_faers_quarters(FAERS_DIR, years = 2022:2026)

# Primary analysis window: 2022Q1-2026Q2 (faricimab approval to latest
# available quarter). To use the full period (e.g. time-aligned
# sensitivity), set WIN_START/WIN_END to NULL.
WIN_START <- "2022-01-01"
WIN_END   <- "2026-06-30"

if (!is.null(WIN_START) && !is.null(WIN_END)) {
  win_ids <- faers_obj$demo %>%
    mutate(.dt = parse_any_date(fda_dt)) %>%
    filter(.dt >= ymd(WIN_START), .dt <= ymd(WIN_END)) %>%
    pull(primaryid)
  faers_obj <- lapply(faers_obj, function(tb) {
    if ("primaryid" %in% names(tb)) filter(tb, primaryid %in% win_ids) else tb
  })
  message("Window filter: ", WIN_START, " to ", WIN_END,
          " -> ", length(win_ids), " reports kept")
}

demo <- as_tibble(faers_obj$demo)
drug <- as_tibble(faers_obj$drug)
reac <- as_tibble(faers_obj$reac)

# Some cached quarters may lack outc / ther / indi (legacy AERS quarters
# only carry demo/drug/reac).
get_table <- function(nm) {
  if (nm %in% names(faers_obj) && !is.null(faers_obj[[nm]])) {
    as_tibble(faers_obj[[nm]])
  } else {
    warning("Table '", nm, "' missing from FAERS object; related analyses will be skipped or degraded",
            call. = FALSE)
    NULL
  }
}
outc <- get_table("outc")   # serious outcomes (death/hospitalisation)
ther <- get_table("ther")   # needed for TTO
indi <- get_table("indi")   # without it, indication = "other/unknown"
has_outc <- !is.null(outc)

# ------------------------------------------------------------
# Defensive column normalisation: modern FAERS calls the sequence
# column dsg_drug_seq (THER) / indi_drug_seq (INDI); legacy AERS uses
# drug_seq. Normalise regardless of upstream naming so downstream TTO
# can join on primaryid + drug_seq. Also coerce drug_seq to character
# to avoid numeric/character join failures.
# ------------------------------------------------------------
norm_seq_col <- function(df, aliases) {
  if (is.null(df)) return(NULL)
  hit <- intersect(aliases, names(df))
  if (!"drug_seq" %in% names(df) && length(hit) > 0) {
    names(df)[names(df) == hit[1]] <- "drug_seq"
  }
  if ("drug_seq" %in% names(df)) {
    df$drug_seq <- as.character(df$drug_seq)
  }
  df
}
ther <- norm_seq_col(ther, c("dsg_drug_seq", "drug_seq"))
indi <- norm_seq_col(indi, c("indi_drug_seq", "drug_seq"))
drug <- norm_seq_col(drug, c("dsg_drug_seq", "drug_seq"))

# MedDRA PT -> primary SOC map, built from the MedDRA 27.1 release zip
# (parsed once and cached as data/meddra_pt_soc_map.csv; the cache is
# reused afterwards). Accepts the zip or an extracted folder.
MEDDRA_ZIP <- Sys.getenv("MEDDRA_ZIP", unset = "data/MedDRA_27_1_English.zip")
pt_soc_map <- build_pt_soc_map(MEDDRA_ZIP)

# ------------------------------------------------------------
# 0b. Exclude non-AE PTs (OFF LABEL USE / NO ADVERSE EVENT /
#     product issues / procedural errors / the indication itself...).
#     The list lives in data/excluded_pts.csv and is editable; it is
#     applied to both target and background so the 2x2 denominators
#     stay consistent. Set to FALSE to disable.
# ------------------------------------------------------------
EXCLUDE_NON_AE_PT <- TRUE
excluded_pts <- character(0)
if (EXCLUDE_NON_AE_PT) {
  if (!file.exists("data/excluded_pts.csv")) {
    stop("data/excluded_pts.csv not found!\n",
         "  Current working directory: ", getwd(), "\n",
         "  Please setwd() to the project root (the folder containing data/ and R/)")
  }
  excluded_pts <- read_csv("data/excluded_pts.csv", show_col_types = FALSE)$pt %>%
    toupper() %>% trimws() %>% unique()
  message("Non-AE PT exclusion enabled: ", length(excluded_pts),
          " PTs (data/excluded_pts.csv)")
}

# Preprocess reac once (uppercase + exclusion); PT/SOC/TTO all use it.
n_before <- nrow(reac)
reac_pt <- reac %>%
  mutate(pt = toupper(pt)) %>%
  filter(!pt %in% excluded_pts)
message("reac: ", n_before, " rows -> ", nrow(reac_pt),
        " after non-AE PT exclusion (removed ", n_before - nrow(reac_pt), ")")

# ------------------------------------------------------------
# 1. Drug identification (brand + generic + biosimilars)
# ------------------------------------------------------------
DRUG_PATTERNS <- list(
  faricimab  = "FARICIMAB|VABYSMO",
  aflibercept= "AFLIBERCEPT|EYLEA",
  ranibizumab= "RANIBIZUMAB|LUCENTIS|BYOOVIZ|CIMERLI|RANIVIZIO"
)
# Full anti-VEGF list used for the masking background (incl. competitors)
ANTIVEGF_ALL <- paste(
  "FARICIMAB|VABYSMO|AFLIBERCEPT|EYLEA|RANIBIZUMAB|LUCENTIS|BYOOVIZ|CIMERLI|",
  "RANIVIZIO|BROLUCIZUMAB|BEOVU|PEGAPTANIB|MACUGEN|ABICIPAR|BEVACIZUMAB|AVASTIN",
  sep = ""
)

classify_drug <- function(nm) {
  nm <- toupper(nm)
  case_when(
    str_detect(nm, DRUG_PATTERNS$faricimab)   ~ "faricimab",
    str_detect(nm, DRUG_PATTERNS$aflibercept) ~ "aflibercept",
    str_detect(nm, DRUG_PATTERNS$ranibizumab) ~ "ranibizumab",
    TRUE ~ NA_character_
  )
}

drug2 <- drug %>%
  mutate(
    drug_std = classify_drug(coalesce(drugname, "")),
    is_antivegf = str_detect(toupper(coalesce(drugname, "")), ANTIVEGF_ALL)
  )

# Reports where a target drug is primary suspect (PS)
ps_index <- drug2 %>%
  filter(role_cod == "PS", !is.na(drug_std)) %>%
  distinct(primaryid, drug_std)

# Drop reports containing >1 target drug (any role) to keep the
# head-to-head comparisons pure.
multi_idx <- drug2 %>%
  filter(!is.na(drug_std)) %>%
  distinct(primaryid, drug_std) %>%
  count(primaryid) %>% filter(n > 1)

cohort <- ps_index %>% anti_join(multi_idx, by = "primaryid")

# Masking set: any report containing any anti-VEGF in any role
antivegf_any <- drug2 %>%
  filter(is_antivegf) %>%
  distinct(primaryid)

# ------------------------------------------------------------
# 2. Report-level covariate table
# ------------------------------------------------------------
parse_fda_dt <- function(x) {
  x <- as.character(x)
  x[nchar(x) != 8] <- NA
  ymd(x, quiet = TRUE)
}

if (!is.null(outc)) {
  serious_flags <- outc %>%
    mutate(code = toupper(outc_cod)) %>%
    summarise(
      death  = any(code == "DE"),
      life_threat = any(code == "LT"),
      hosp = any(code == "HO"),
      disability = any(code == "DS"),
      other_serious = any(code %in% c("OT", "RI", "CA")),
      .by = primaryid
    )
} else {
  message("!! No outc table: death/hospitalisation outcomes set to NA. ",
          "Re-cache quarters with outc retained to enable this analysis.")
  serious_flags <- cohort %>%
    transmute(primaryid, death = NA, life_threat = NA,
              hosp = NA, disability = NA, other_serious = NA)
}

# Indication groups from INDI (PT text matching).
# Performance: semi_join down to the cohort first (full-database INDI
# has ~25M rows vs tens of thousands in the cohort), and vectorise
# str_detect before summarising — never call regex row-by-row inside
# summarise (orders of magnitude slower).
if (!is.null(indi)) {
  indi_grp <- indi %>%
    semi_join(cohort, by = "primaryid") %>%
    mutate(ip = toupper(indi_pt),
           is_dme  = str_detect(ip, "DIABETIC"),
           is_namd = str_detect(ip, "MACULAR DEGENERATION|NEOVASCULAR AGE"),
           is_rvo  = str_detect(ip, "VEIN OCCLUSION|RETINAL VEIN")) %>%
    summarise(
      indication = case_when(
        any(is_dme)  ~ "DME/DR",
        any(is_namd) ~ "nAMD",
        any(is_rvo)  ~ "RVO",
        TRUE ~ "other/unknown"
      ),
      .by = primaryid
    )
} else {
  indi_grp <- cohort %>% transmute(primaryid, indication = "other/unknown")
}

cohort_reports <- cohort %>%
  left_join(demo %>% dplyr::select(primaryid,
                            any_of(c("fda_dt", "age", "age_cod", "sex",
                                     "occp_cod", "occr_country", "reporter_country"))),
            by = "primaryid") %>%
  left_join(serious_flags, by = "primaryid") %>%
  left_join(indi_grp, by = "primaryid")

# Column-completeness fallback: add NA for any demographic column the
# cache lacks so downstream code does not error.
DEMO_EXPECTED <- c("fda_dt", "age", "age_cod", "sex",
                   "occp_cod", "occr_country", "reporter_country")
for (cc in setdiff(DEMO_EXPECTED, names(cohort_reports))) {
  cohort_reports[[cc]] <- NA_character_
}
if (all(is.na(cohort_reports$age)) && all(is.na(cohort_reports$sex))) {
  message("!! Cached demo table has no demographic fields (age/sex/country): ",
          "the corresponding covariates are unavailable for the baseline table ",
          "and adjusted models.")
}

cohort_reports <- cohort_reports %>%
  mutate(
    fda_date = parse_fda_dt(fda_dt),
    report_year = year(fda_date),
    age_years = case_when(
      age_cod == "YR" ~ as.numeric(age),
      age_cod == "DEC" ~ as.numeric(age) * 10,
      age_cod == "MON" ~ as.numeric(age) / 12,
      age_cod == "WK" ~ as.numeric(age) / 52,
      TRUE ~ NA_real_
    ),
    age_grp = case_when(
      is.na(age_years) ~ "unknown",
      age_years < 18 ~ "<18",
      age_years <= 64 ~ "18-64",
      TRUE ~ ">=65"
    ),
    sex = case_when(sex %in% c("F", "M") ~ sex, TRUE ~ "unknown"),
    country = toupper(coalesce(occr_country, reporter_country)),
    region = if_else(country == "US", "US", "non-US"),
    reporter_type = case_when(
      occp_cod %in% c("MD", "PH", "RN", "OT") ~ "health_professional",
      occp_cod == "CN" ~ "consumer",
      TRUE ~ "unknown"
    ),
    # With outc present, NA means the outcome was not reported; without
    # outc, keep NA as a "data missing" marker.
    across(c(death, life_threat, hosp, disability, other_serious),
           ~ if (has_outc) coalesce(.x, FALSE) else .x)
  ) %>%
  rename(drug = drug_std)

# ------------------------------------------------------------
# 3. PT / SOC long tables (target cohort + full-database background)
# ------------------------------------------------------------
reac2 <- reac_pt %>%
  left_join(pt_soc_map, by = "pt") %>%
  distinct(primaryid, pt, soc)

pairs_pt <- reac2 %>%
  inner_join(cohort, by = "primaryid") %>%
  transmute(case_id = primaryid, drug = drug_std, event = pt) %>%
  distinct()

# Verify the exclusion took effect: no excluded PT should remain
n_leak <- sum(pairs_pt$event %in% excluded_pts)
if (n_leak > 0) {
  warning("!! Excluded PTs still present in pairs_pt (", n_leak, " rows)", call. = FALSE)
} else if (length(excluded_pts) > 0) {
  message("OK: no excluded PTs remain in pairs_pt")
}

pairs_soc <- reac2 %>%
  filter(!is.na(soc)) %>%
  inner_join(cohort, by = "primaryid") %>%
  transmute(case_id = primaryid, drug = drug_std, event = soc) %>%
  distinct()

# Full-database background: reports without any of the three drugs in
# any role, so competing class drugs never enter the background.
bg_ids <- demo %>%
  anti_join(antivegf_any, by = "primaryid") %>%
  distinct(primaryid)

bg_pt <- reac2 %>%
  inner_join(bg_ids, by = "primaryid") %>%
  transmute(case_id = primaryid, event = pt) %>% distinct()

bg_soc <- reac2 %>%
  filter(!is.na(soc)) %>%
  inner_join(bg_ids, by = "primaryid") %>%
  transmute(case_id = primaryid, event = soc) %>% distinct()

# ------------------------------------------------------------
# 4. Time-to-onset (TTO): target drug (PS) therapy start
#    (THER.start_dt) -> event date (DEMO.event_dt).
#    Must align on primaryid + drug_seq so the PS drug's own course is
#    used — otherwise a concomitant drug's start date leaks in.
# ------------------------------------------------------------
parse_partial_dt <- function(x) {
  x <- as.character(x)
  x[nchar(x) != 8] <- NA        # keep complete YYYYMMDD only
  ymd(x, quiet = TRUE)
}

# Precise path: ther has drug_seq -> align by case + course
tto_precise <- !is.null(ther) && "event_dt" %in% names(demo) &&
  "drug_seq" %in% names(ther) && "drug_seq" %in% names(drug2)
# Fallback path: drug_seq was trimmed from the cache -> align by case,
# using the earliest drug start date (approximation; reasonable when the
# anti-VEGF is the main therapy; report as sensitivity/approximate).
tto_fallback <- !tto_precise && !is.null(ther) &&
  "event_dt" %in% names(demo) && "start_dt" %in% names(ther)

if (tto_precise) {
  message("TTO: precise mode (primaryid + drug_seq course alignment)")
  ps_seq <- drug2 %>%
    filter(role_cod == "PS", !is.na(drug_std)) %>%
    semi_join(cohort, by = c("primaryid", "drug_std")) %>%
    dplyr::select(primaryid, drug_std, drug_seq)

  tto <- ps_seq %>%
    inner_join(ther %>% dplyr::select(primaryid, drug_seq, start_dt),
               by = c("primaryid", "drug_seq")) %>%
    inner_join(demo %>% dplyr::select(primaryid, any_of("event_dt")), by = "primaryid") %>%
    mutate(
      start = parse_partial_dt(start_dt),
      onset = parse_partial_dt(event_dt),
      tto_days = as.integer(onset - start)
    ) %>%
    filter(!is.na(tto_days), tto_days >= 0, tto_days <= 1825) %>%  # cap at 5 years
    # keep the earliest course per case
    group_by(primaryid) %>%
    slice_min(tto_days, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    inner_join(reac_pt %>% distinct(primaryid, pt),
               by = "primaryid", relationship = "many-to-many") %>%
    transmute(case_id = primaryid, drug_std, pt, tto_days)
} else if (tto_fallback) {
  message("TTO: fallback mode (ther lacks drug_seq; earliest start_dt per primaryid)")
  ther_min <- ther %>%
    mutate(start = parse_partial_dt(start_dt)) %>%
    filter(!is.na(start)) %>%
    summarise(start = min(start), .by = primaryid)

  tto <- cohort %>%
    dplyr::select(primaryid, drug_std) %>%
    inner_join(ther_min, by = "primaryid") %>%
    inner_join(demo %>% dplyr::select(primaryid, any_of("event_dt")), by = "primaryid") %>%
    mutate(
      onset = parse_partial_dt(event_dt),
      tto_days = as.integer(onset - start)
    ) %>%
    filter(!is.na(tto_days), tto_days >= 0, tto_days <= 1825) %>%
    inner_join(reac_pt %>% distinct(primaryid, pt),
               by = "primaryid", relationship = "many-to-many") %>%
    transmute(case_id = primaryid, drug_std, pt, tto_days)
} else {
  message("!! TTO requirements not met, skipping:")
  if (is.null(ther)) message("   - no ther table")
  if (!is.null(ther))
    message("   - ther columns: ", paste(names(ther), collapse = ", "))
  if (!"event_dt" %in% names(demo)) message("   - demo lacks event_dt")
  tto <- NULL
}

# ------------------------------------------------------------
# 5. All concomitant drugs of cohort reports (for adjusted models)
# ------------------------------------------------------------
drug_all <- drug2 %>%
  semi_join(cohort, by = "primaryid") %>%
  dplyr::select(primaryid, any_of(c("drug_seq", "role_cod", "drugname",
                             "drug_std", "is_antivegf")))

# ------------------------------------------------------------
# 6. Save
# ------------------------------------------------------------
dir.create("workdir", showWarnings = FALSE)
# Build fingerprint: 02-06 use this to detect a stale workdir.
saveRDS(list(built_at = as.character(Sys.time()),
             exclude_on = EXCLUDE_NON_AE_PT,
             excluded_pts_n = length(excluded_pts),
             n_cohort = nrow(cohort)),
        "workdir/_build_info.rds")
# Denominator: full cohort size (including cases left with no PTs after
# exclusion) — the standard 2x2 denominator for 02.
saveRDS(cohort %>% count(drug_std, name = "n"), "workdir/cohort_n.rds")
saveRDS(nrow(bg_ids), "workdir/bg_n.rds")
saveRDS(cohort_reports, "workdir/cohort_reports.rds")
saveRDS(pairs_pt,       "workdir/pairs_pt.rds")
saveRDS(pairs_soc,      "workdir/pairs_soc.rds")
saveRDS(bg_pt,          "workdir/bg_pt.rds")
saveRDS(bg_soc,         "workdir/bg_soc.rds")
if (!is.null(tto)) saveRDS(tto, "workdir/tto.rds")
saveRDS(antivegf_any,   "workdir/antivegf_any.rds")
saveRDS(drug_all,       "workdir/drug_all.rds")

message("=== Cohort size ===")
print(cohort %>% count(drug_std))
message("=== 01_faers_prepare done ===")
