# ============================================================
# 09_cross_db_validation.R
# Cross-database signal confirmation (the core narrative table)
#
# Evidence grading (suggested wording for Methods):
#   Confirmed   : FAERS-positive AND at least one external database
#                 significant in the same direction
#   Consistent  : FAERS-positive AND external database in the same
#                 direction (ROR>1) but below the significance threshold
#   Unconfirmed : FAERS-positive but not reproduced externally / no
#                 external data
#   Opposite    : external database ROR<1; requires discussion
#
# Prerequisites: 02 / 07 / 08 already run (missing JADER/EV sources are
# skipped automatically)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(purrr); library(ggplot2); library(stringr)
})
source("R/00_utils.R")

dir.create("results/cross_db", showWarnings = FALSE, recursive = TRUE)

safe_read <- function(path) if (file.exists(path)) read_csv(path, show_col_types = FALSE) else NULL

# ------------------------------------------------------------
# 1. Pool head-to-head results across databases
#    (unified schema: source, comparison, level, event, a, ROR, CI;
#     missing sources are skipped)
# ------------------------------------------------------------
norm_src <- function(df, src, lvl) {
  if (is.null(df)) { message("Source missing, skipped: ", src); return(NULL) }
  if (!"event" %in% names(df) && "soc" %in% names(df)) {
    df <- df %>% rename(event = soc)
  }
  df %>% mutate(source = src, level = lvl)
}

faers_h2h_soc <- norm_src(safe_read("results/faers/B1_soc_headtohead.csv"), "FAERS", "SOC")
faers_h2h_pt  <- norm_src(safe_read("results/faers/B2_pt_headtohead.csv"), "FAERS", "PT")
jader_h2h     <- norm_src(safe_read("results/jader/pt_headtohead.csv"), "JADER", "PT")
ev_h2h        <- norm_src(safe_read("results/eudravigilance/soc_headtohead.csv"),
                          "EudraVigilance", "SOC")

src_list <- compact(list(faers_h2h_soc, faers_h2h_pt, jader_h2h, ev_h2h))
if (length(src_list) == 0) stop("No sources read - run 02/07/08 first")

all_res <- bind_rows(lapply(src_list, function(df) {
  df %>% dplyr::select(source, level, comparison, event, a, c, ROR, ROR025, ROR975)
})) %>%
  mutate(event = toupper(event))

write_tab(all_res, "results/cross_db/all_sources_pooled.csv")

# ------------------------------------------------------------
# 2. Evidence grading (bidirectional: ROR>1 = higher in the first drug,
#    ROR<1 = higher in the second drug)
# ------------------------------------------------------------
evidence <- all_res %>%
  mutate(
    # Significant direction: CI excludes 1 and the higher-side count >= 3
    dir_sig = case_when(
      ROR025 > 1 & a >= 3 ~ word(comparison, 1),
      ROR975 < 1 & c >= 3 ~ word(comparison, 3),
      TRUE ~ "NS"),
    # Point-estimate direction (regardless of significance)
    dir_pt = case_when(
      ROR > 1 ~ word(comparison, 1),
      ROR < 1 ~ word(comparison, 3),
      TRUE ~ "NS")
  ) %>%
  group_by(level, comparison, event) %>%
  summarise(
    faers_dir   = dir_sig[source == "FAERS"][1],
    faers_ror   = ROR[source == "FAERS"][1],
    faers_lo    = ROR025[source == "FAERS"][1],
    faers_hi    = ROR975[source == "FAERS"][1],
    ext_sig     = any(source != "FAERS" & dir_sig != "NS" &
                      dir_sig == faers_dir, na.rm = TRUE),
    ext_dir_pos = any(source != "FAERS" & dir_pt != "NS" &
                      dir_pt == dir_pt[source == "FAERS"][1], na.rm = TRUE),
    ext_opp     = any(source != "FAERS" & dir_sig != "NS" &
                      faers_dir != "NS" & dir_sig != faers_dir, na.rm = TRUE),
    n_sources   = n_distinct(source),
    .groups = "drop"
  ) %>%
  mutate(
    faers_pos = !is.na(faers_dir) & faers_dir != "NS",
    evidence_level = case_when(
      !faers_pos                          ~ "FAERS negative",
      faers_pos & ext_sig                 ~ "Confirmed (externally replicated)",
      faers_pos & ext_dir_pos & !ext_sig  ~ "Consistent (same direction)",
      faers_pos & ext_opp                 ~ "Opposite (requires discussion)",
      TRUE                                ~ "Unconfirmed (FAERS only)"
    ))

write_tab(evidence, "results/cross_db/evidence_classification.csv")

# ------------------------------------------------------------
# 3. Multi-database forest plot for key events.
#    Defaults to the systemic SOCs positive in FAERS; switch to PT
#    level if needed.
# ------------------------------------------------------------
focus_soc <- c("BLOOD AND LYMPHATIC SYSTEM DISORDERS", "HEPATOBILIARY DISORDERS",
               "RENAL AND URINARY DISORDERS", "METABOLISM AND NUTRITION DISORDERS",
               "PSYCHIATRIC DISORDERS", "CARDIAC DISORDERS",
               "NERVOUS SYSTEM DISORDERS", "VASCULAR DISORDERS")

plot_df <- all_res %>%
  filter(level == "SOC", event %in% focus_soc,
         comparison %in% c("aflibercept vs faricimab", "ranibizumab vs faricimab"))

for (cp in unique(plot_df$comparison)) {
  d <- plot_df %>% filter(comparison == cp)
  ggsave(paste0("figures/cross_db_", gsub(" ", "_", cp), ".pdf"),
         plot_cross_db(d, event_order = focus_soc),
         width = 12, height = 7)
}

message("=== 09_cross_db_validation done ===")
message("Key output: results/cross_db/evidence_classification.csv")
