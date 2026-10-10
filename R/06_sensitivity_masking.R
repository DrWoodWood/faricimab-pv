# ============================================================
# 06_sensitivity_masking.R
# Masking / competition-bias sensitivity analysis
#   The three drugs act as competing reporting sources and can mask
#   each other's signals. Approach: recompute each drug vs background
#   with all reports containing any anti-VEGF (any role) removed from
#   the background, then compare consistency with the unmasked analysis.
#
# Prerequisite: 01_faers_prepare.R (its bg_pt/bg_soc are already the
# "masked" background; for completeness this script also builds the
# "unmasked" full-database background for comparison)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(purrr); library(ggplot2); library(lubridate)
})
source("R/00_utils.R")

pairs_pt  <- readRDS("workdir/pairs_pt.rds")
pairs_soc <- readRDS("workdir/pairs_soc.rds")
bg_pt_masked  <- readRDS("workdir/bg_pt.rds")     # anti-VEGF reports already removed
bg_soc_masked <- readRDS("workdir/bg_soc.rds")
assert_fresh_workdir()   # errors out on a stale workdir

# Unmasked background: full database (including other anti-VEGF reports)
FAERS_DIR <- Sys.getenv("FAERS_DIR", unset = "data/rds_cache_full")
faers_obj <- load_faers_quarters(FAERS_DIR, years = 2022:2026)

# Same analysis window as 01
WIN_START <- "2022-01-01"
WIN_END   <- "2026-06-30"
win_ids <- faers_obj$demo %>%
  mutate(.dt = parse_any_date(fda_dt)) %>%
  filter(.dt >= ymd(WIN_START), .dt <= ymd(WIN_END)) %>%
  pull(primaryid)
faers_obj <- lapply(faers_obj, function(tb) {
  if ("primaryid" %in% names(tb)) filter(tb, primaryid %in% win_ids) else tb
})

reac <- as_tibble(faers_obj$reac)
demo <- as_tibble(faers_obj$demo)
# Accepts the MedDRA zip or an extracted folder
MEDDRA_ZIP <- Sys.getenv("MEDDRA_ZIP", unset = "data/MedDRA_27_1_English.zip")
pt_soc_map <- build_pt_soc_map(MEDDRA_ZIP)

cohort_ids <- pairs_pt %>% transmute(primaryid = as.character(case_id)) %>% distinct()

reac2 <- reac %>% mutate(pt = toupper(pt)) %>%
  left_join(pt_soc_map, by = "pt") %>% distinct(primaryid, pt, soc)

bg_pt_full <- reac2 %>%
  anti_join(cohort_ids, by = "primaryid") %>%
  transmute(case_id = primaryid, event = pt) %>% distinct()
bg_soc_full <- reac2 %>% filter(!is.na(soc)) %>%
  anti_join(cohort_ids, by = "primaryid") %>%
  transmute(case_id = primaryid, event = soc) %>% distinct()

dir.create("results/masking", showWarnings = FALSE, recursive = TRUE)
DRUGS <- c("faricimab", "aflibercept", "ranibizumab")

compare_backgrounds <- function(pairs, bg_masked, bg_full, level) {
  res_masked <- map_dfr(DRUGS, function(dg) {
    tgt <- pairs %>% filter(drug == dg) %>% distinct(case_id, event)
    disprop_all(tgt, bg_masked, n_min = 3) %>% mutate(drug = dg)
  }) %>% mutate(background = "masked (anti-VEGF reports removed)")

  res_full <- map_dfr(DRUGS, function(dg) {
    tgt <- pairs %>% filter(drug == dg) %>% distinct(case_id, event)
    disprop_all(tgt, bg_full, n_min = 3) %>% mutate(drug = dg)
  }) %>% mutate(background = "unmasked (full database)")

  bind_rows(res_masked, res_full) %>% mutate(level = level, .before = 1)
}

cmp_soc <- compare_backgrounds(pairs_soc, bg_soc_masked, bg_soc_full, "SOC")
cmp_pt  <- compare_backgrounds(pairs_pt,  bg_pt_masked,  bg_pt_full,  "PT")

write_tab(cmp_soc, "results/masking/soc_masked_vs_unmasked.csv")
write_tab(cmp_pt,  "results/masking/pt_masked_vs_unmasked.csv")

# Consistency scatter: log2(ROR_masked) vs log2(ROR_unmasked)
plot_consistency <- function(df, level_lab) {
  d <- df %>%
    dplyr::select(drug, event, background, ROR) %>%
    pivot_wider(names_from = background, values_from = ROR,
                names_prefix = "r_") %>%
    setNames(c("drug", "event", "masked", "unmasked")) %>%
    filter(is.finite(masked), is.finite(unmasked), masked > 0, unmasked > 0)
  rho <- cor(log2(d$masked), log2(d$unmasked), method = "spearman")
  ggplot(d, aes(x = log2(unmasked), y = log2(masked), colour = drug)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_point(alpha = 0.5) +
    labs(title = paste0("Masking sensitivity (", level_lab, "), Spearman rho = ",
                        round(rho, 3)),
         x = "log2 ROR (unmasked background)", y = "log2 ROR (masked background)") +
    theme_bw(base_size = 11)
}

ggsave("figures/masking_consistency_soc.pdf",
       plot_consistency(cmp_soc, "SOC"), width = 7, height = 6)
ggsave("figures/masking_consistency_pt.pdf",
       plot_consistency(cmp_pt, "PT"), width = 7, height = 6)

message("=== 06_sensitivity_masking done ===")
