# ============================================================
# 07_jader_analysis.R
# JADER (Japanese PMDA adverse-event database) external validation
#
# Data acquisition:
#   PMDA website -> safety info services -> suspected ADR database
#   https://www.pmda.go.jp/safety/info-services/drugs/adr-info/suspected-adr/0003.html
#   Download the yearly CSV packs (demo/drug/reac/hist) of the
#   suspected-ADR database and place the extracted files in data/jader/
#   (all years in one folder).
#
# Notes:
#   - Older yearly files are Shift-JIS (CP932); newer ones are UTF-8.
#     The reader tries UTF-8 first and falls back to CP932.
#   - JADER adverse events are MedDRA/J Japanese PTs; the EN mapping
#     lives in data/meddra_j_en_map.csv (for publication, export a full
#     JP-EN mapping with the MedDRA Desktop Browser).
#   - Faricimab was approved in Japan in March 2022 (Vabysmo);
#     a 2022-04 onwards window is recommended.
# ============================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(purrr); library(tidyr)
})
source("R/00_utils.R")

JADER_DIR <- "data/jader"
dir.create("results/jader", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Read and merge all yearly files (encoding handled automatically)
# ------------------------------------------------------------
read_jader_file <- function(path) {
  df <- tryCatch(
    read_csv(path, locale = locale(encoding = "UTF-8"),
             show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  # fall back to CP932 if UTF-8 fails or column names show mojibake
  if (is.null(df) || any(str_detect(names(df), ""))) {
    df <- read_csv(path, locale = locale(encoding = "CP932"),
                   show_col_types = FALSE, progress = FALSE)
  }
  df
}

read_jader_table <- function(pattern) {
  files <- list.files(JADER_DIR, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("JADER files not found: ", pattern,
                               " - download them into ", JADER_DIR)
  map_dfr(files, read_jader_file)
}

demo_j <- read_jader_table("^demo.*\\.csv$")
drug_j <- read_jader_table("^drug.*\\.csv$")
reac_j <- read_jader_table("^reac.*\\.csv$")

# ------------------------------------------------------------
# 2. Standardise column names (they vary slightly across years;
#    matched by keyword). Japanese column names are kept as-is because
#    they are the actual headers in the source files.
# ------------------------------------------------------------
pick_col <- function(df, patterns, required = TRUE) {
  for (p in patterns) {
    hit <- names(df)[str_detect(names(df), p)]
    if (length(hit) >= 1) return(hit[1])
  }
  if (required) stop("Column not found: ", paste(patterns, collapse = " / "))
  return(NA_character_)   # optional column absent -> NA
}

demo_j <- demo_j %>% rename(
  case_id  = all_of(pick_col(demo_j, c("識別番号"))),
  sex      = all_of(pick_col(demo_j, c("^性別$", "性別"))),
  age      = all_of(pick_col(demo_j, c("^年齢$", "年齢")))
)
# The public-download JADER has no report date column (only
# 報告年度・四半期 = reporting year/quarter), so it is optional.
rd_col <- pick_col(demo_j,
                   c("報告・受付年月日", "受付年月日", "報告年月日",
                     "報告年度・四半期", "報告年度"),
                   required = FALSE)
if (!is.na(rd_col)) {
  demo_j <- demo_j %>% rename(report_dt = all_of(rd_col))
} else {
  demo_j <- demo_j %>% mutate(report_dt = NA_character_)
  message("Note: no report-date column in JADER demo; report_dt set to NA (not used downstream)")
}

drug_j <- drug_j %>% rename(
  case_id   = all_of(pick_col(drug_j, c("識別番号"))),
  generic   = all_of(pick_col(drug_j, c("医薬品（一般名）", "一般名"))),
  brand     = all_of(pick_col(drug_j, c("医薬品（商品名）", "医薬品（販売名）", "商品名", "販売名"))),
  involve   = all_of(pick_col(drug_j, c("医薬品の関与", "関与")))
)

reac_j <- reac_j %>% rename(
  case_id   = all_of(pick_col(reac_j, c("識別番号"))),
  pt_j      = all_of(pick_col(reac_j, c("^有害事象$", "有害事象"))),
  outcome   = all_of(pick_col(reac_j, c("^転帰$", "転帰")))
)

# ------------------------------------------------------------
# 3. Drug identification (suspected drugs only).
#    Japanese strings below are data, not comments — do not translate.
# ------------------------------------------------------------
drug_j <- drug_j %>%
  mutate(
    nm = paste0(coalesce(generic, ""), "|", coalesce(brand, "")),
    drug_std = case_when(
      str_detect(nm, "ファリシマブ|バビスモ")                ~ "faricimab",
      str_detect(nm, "アフリベルセプト|アイリーア")           ~ "aflibercept",
      str_detect(nm, "ラニビズマブ|ルセンティス")             ~ "ranibizumab",
      TRUE ~ NA_character_
    ),
    is_suspect = str_detect(coalesce(involve, ""), "被疑"),
    is_antivegf = str_detect(nm,
      "ファリシマブ|バビスモ|アフリベルセプト|アイリーア|ラニビズマブ|ルセンティス|ブロルシズマブ|ベオビュー|ベバシズマブ|アバスチン")
  )

# Drop reports containing >1 target drug (any role)
multi_idx <- drug_j %>% filter(!is.na(drug_std)) %>%
  distinct(case_id, drug_std) %>% count(case_id) %>% filter(n > 1)

cohort_j <- drug_j %>%
  filter(is_suspect, !is.na(drug_std)) %>%
  distinct(case_id, drug_std) %>%
  anti_join(multi_idx, by = "case_id")

# Background: JADER reports without any anti-VEGF
avegf_ids <- drug_j %>% filter(is_antivegf) %>% distinct(case_id)
bg_ids <- demo_j %>% distinct(case_id) %>% anti_join(avegf_ids, by = "case_id")

# ------------------------------------------------------------
# 4. JP->EN PT mapping + long tables
# ------------------------------------------------------------
pt_map <- read_csv("data/meddra_j_en_map.csv", show_col_types = FALSE) %>%
  transmute(pt_j, pt_en = toupper(pt_en))

reac_joined <- reac_j %>%
  left_join(pt_map, by = "pt_j")
unmapped <- reac_joined %>% filter(is.na(pt_en)) %>% count(pt_j, sort = TRUE)
write_tab(unmapped, "results/jader/unmapped_pt_needing_translation.csv")
message("PT mapping coverage: ", round(mean(!is.na(reac_joined$pt_en)) * 100, 1),
        "% (", sum(!is.na(reac_joined$pt_en)), "/", nrow(reac_joined), " rows)")
reac2 <- reac_joined %>%
  filter(!is.na(pt_en))                       # keep mapped PTs only (validation events)

pairs_j <- reac2 %>%
  inner_join(cohort_j, by = "case_id") %>%
  transmute(case_id, drug = drug_std, event = pt_en) %>% distinct()

bg_j <- reac2 %>%
  inner_join(bg_ids, by = "case_id") %>%
  transmute(case_id, event = pt_en) %>% distinct()

# ------------------------------------------------------------
# 5. Disproportionality: each drug vs JADER background + head-to-head.
#    Denominators use the full cohort / full background (as in 02):
#    event counts a/c only cover mapped PTs, but b/d use all reports,
#    so incomplete mapping does not inflate RORs.
# ------------------------------------------------------------
DRUGS <- c("faricimab", "aflibercept", "ranibizumab")
cohort_n_j <- cohort_j %>% count(drug_std, name = "n")
bg_n_j <- nrow(bg_ids)
message("JADER cohort: ", paste(cohort_n_j$drug_std, cohort_n_j$n, sep = "=", collapse = ", "),
        " | background reports: ", bg_n_j)
jader_vsbg <- map_dfr(DRUGS, function(dg) {
  tgt <- pairs_j %>% filter(drug == dg) %>% distinct(case_id, event)
  if (nrow(tgt) == 0) return(NULL)
  disprop_all(tgt, bg_j, n_min = 3,
              n_target_total = cohort_n_j$n[cohort_n_j$drug_std == dg],
              n_bg_total = bg_n_j) %>% mutate(drug = dg, .before = 1)
}) %>% flag_signal(min_n = 3) %>% add_fisher()

cmp <- list(c("aflibercept", "faricimab"),
            c("ranibizumab", "faricimab"),
            c("ranibizumab", "aflibercept"))
jader_h2h <- map_dfr(cmp, function(cp) {
  if (!all(cp %in% pairs_j$drug)) return(NULL)
  head_to_head(pairs_j, cp[1], cp[2], n_min = 3,
               n_A = cohort_n_j$n[cohort_n_j$drug_std == cp[1]],
               n_B = cohort_n_j$n[cohort_n_j$drug_std == cp[2]])
}) %>%
  group_by(comparison) %>%
  group_modify(~flag_signal(.x, min_n = 3)) %>% ungroup() %>%
  add_fisher() %>%
  # Bidirectional direction column: CI entirely >1 -> higher in the
  # first drug; entirely <1 -> higher in the second drug
  mutate(direction = case_when(
    ROR025 > 1 ~ paste0("higher in ", word(comparison, 1)),
    ROR975 < 1 ~ paste0("higher in ", word(comparison, 3)),
    TRUE ~ "NS"
  ), .after = signal)

write_tab(jader_vsbg, "results/jader/pt_vs_jader_background.csv")
write_tab(jader_h2h,  "results/jader/pt_headtohead.csv")

# ------------------------------------------------------------
# 6. Outcome (death) comparison
# ------------------------------------------------------------
death_j <- reac_j %>%
  filter(str_detect(coalesce(outcome, ""), "死亡")) %>%
  distinct(case_id)

death_tab <- map_dfr(cmp, function(cp) {
  rA <- cohort_j %>% filter(drug_std == cp[1])
  rB <- cohort_j %>% filter(drug_std == cp[2])
  a <- sum(rA$case_id %in% death_j$case_id); b <- nrow(rA) - a
  c <- sum(rB$case_id %in% death_j$case_id); d <- nrow(rB) - c
  ror <- calc_ror(a, b, c, d)
  tibble(comparison = paste(cp[1], "vs", cp[2]), outcome = "death",
         a = a, b = b, c = c, d = d,
         ROR = ror$ROR, ROR025 = ror$ROR025, ROR975 = ror$ROR975)
})
write_tab(death_tab, "results/jader/death_headtohead.csv")

message("=== 07_jader_analysis done ===")
message("Tip: see results/jader/unmapped_pt_needing_translation.csv; adding ")
message("     high-frequency unmapped PTs to data/meddra_j_en_map.csv improves coverage.")
