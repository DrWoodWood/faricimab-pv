# ============================================================
# 05_sensitivity_time_aligned.R
# Post-marketing time-aligned sensitivity analysis (addresses the
# Weber effect / unequal time-on-market critique)
#
# Design: take each drug's first-3-years post-approval reporting window
# and compare against the contemporaneous full-database background:
#   ranibizumab: FDA approval 2006-06-30 -> 2006-07-01 to 2009-06-30
#   aflibercept: FDA approval 2011-11-18 -> 2011-11-18 to 2014-11-17
#   faricimab  : FDA approval 2022-01-28 -> 2022-01-28 to 2025-01-27
#
# Prerequisite: historical FAERS quarters (from 2006 onwards) must be
# present in FAERS_DIR (download via the faers package or parse the
# legacy ASCII packs from the FDA website into the same .rds structure;
# see 00a3_rebuild_legacy_quarters.R).
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(lubridate); library(purrr)
})
source("R/00_utils.R")
mem.maxVSize(48 * 1024^3)   # full-history merge needs headroom

# Historical data: merge all quarters in the same quarterly cache
# directory (must include the 2006-2014 historical quarters).
FAERS_DIR <- Sys.getenv("FAERS_DIR", unset = "data/rds_cache_full")
faers_hist <- load_faers_quarters(FAERS_DIR)

demo_h <- as_tibble(faers_hist$demo)
drug_h <- as_tibble(faers_hist$drug)
reac_h <- as_tibble(faers_hist$reac)

# Accepts the MedDRA zip or an extracted folder
MEDDRA_ZIP <- Sys.getenv("MEDDRA_ZIP", unset = "data/MedDRA_27_1_English.zip")
pt_soc_map <- build_pt_soc_map(MEDDRA_ZIP)

dir.create("results/time_aligned", showWarnings = FALSE, recursive = TRUE)

DRUG_PATTERNS <- list(
  faricimab   = "FARICIMAB|VABYSMO",
  aflibercept = "AFLIBERCEPT|EYLEA",
  ranibizumab = "RANIBIZUMAB|LUCENTIS|BYOOVIZ|CIMERLI|RANIVIZIO"
)
ANTIVEGF_ALL <- paste(
  "FARICIMAB|VABYSMO|AFLIBERCEPT|EYLEA|RANIBIZUMAB|LUCENTIS|BYOOVIZ|CIMERLI|",
  "RANIVIZIO|BROLUCIZUMAB|BEOVU|PEGAPTANIB|MACUGEN|ABICIPAR|BEVACIZUMAB|AVASTIN", sep = "")

WINDOWS <- tribble(
  ~drug,        ~win_start,     ~win_end,
  "ranibizumab", "2006-07-01",  "2009-06-30",
  "aflibercept", "2011-11-18",  "2014-11-17",
  "faricimab",   "2022-01-28",  "2025-01-27"
) %>% mutate(across(c(win_start, win_end), ymd))

parse_fda_dt <- function(x) {
  x <- as.character(x); x[nchar(x) != 8] <- NA; ymd(x, quiet = TRUE)
}

# ------------------------------------------------------------
# Per-window analysis: target drug (PS) vs contemporaneous full
# database (all anti-VEGF reports removed from the background)
# ------------------------------------------------------------
aligned_results <- pmap_dfr(WINDOWS, function(drug, win_start, win_end) {
  d_demo <- demo_h %>%
    mutate(dt = parse_fda_dt(fda_dt)) %>%
    filter(dt >= win_start, dt <= win_end)

  d_drug <- drug_h %>%
    semi_join(d_demo, by = "primaryid") %>%
    mutate(nm = toupper(coalesce(drugname, "")),
           is_target = str_detect(nm, DRUG_PATTERNS[[drug]]),
           is_antivegf = str_detect(nm, ANTIVEGF_ALL))

  tgt_ids <- d_drug %>% filter(is_target, role_cod == "PS") %>% distinct(primaryid)
  avegf_ids <- d_drug %>% filter(is_antivegf) %>% distinct(primaryid)
  bg_ids <- d_demo %>% anti_join(avegf_ids, by = "primaryid") %>% distinct(primaryid)

  d_reac <- reac_h %>%
    mutate(pt = toupper(pt)) %>%
    left_join(pt_soc_map, by = "pt") %>%
    filter(!is.na(soc)) %>%
    distinct(primaryid, soc)

  tgt <- d_reac %>% inner_join(tgt_ids, by = "primaryid") %>%
    transmute(case_id = primaryid, event = soc) %>% distinct()
  bg  <- d_reac %>% inner_join(bg_ids, by = "primaryid") %>%
    transmute(case_id = primaryid, event = soc) %>% distinct()

  disprop_all(tgt, bg, n_min = 3) %>%
    flag_signal(min_n = 3) %>%
    mutate(drug = drug, era = "first_3_years_post_approval",
           window = paste(win_start, "to", win_end), .before = 1)
})

write_tab(aligned_results, "results/time_aligned/soc_first3years.csv")

# ------------------------------------------------------------
# Side-by-side with the primary analysis (faricimab 2022-2025,
# full-window background)
# ------------------------------------------------------------
primary_faricimab <- read_csv("results/faers/A1_soc_vs_fulldb.csv",
                              show_col_types = FALSE) %>%
  filter(drug == "faricimab") %>%
  mutate(era = "faricimab_current_full_window")

comp_plot_df <- aligned_results %>%
  dplyr::select(drug, event, ROR, ROR025, ROR975, era)

# Time-aligned comparison across the three drugs for key systemic SOCs
systemic_soc <- c("CARDIAC DISORDERS", "VASCULAR DISORDERS",
                  "NERVOUS SYSTEM DISORDERS", "BLOOD AND LYMPHATIC SYSTEM DISORDERS",
                  "RENAL AND URINARY DISORDERS", "HEPATOBILIARY DISORDERS",
                  "METABOLISM AND NUTRITION DISORDERS", "PSYCHIATRIC DISORDERS")

p <- comp_plot_df %>%
  filter(event %in% systemic_soc, is.finite(ROR), ROR025 > 0) %>%
  ggplot(aes(x = ROR, y = drug, colour = drug)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = ROR025, xmax = ROR975), height = 0.2) +
  geom_point(size = 2.4) +
  scale_x_log10() +
  facet_wrap(~event) +
  labs(x = "ROR vs contemporaneous background (log scale)", y = NULL,
       title = "Post-marketing first-3-years aligned comparison") +
  theme_bw(base_size = 10) + theme(legend.position = "none")

ggsave("figures/time_aligned_soc.pdf", p, width = 11, height = 8)

message("=== 05_sensitivity_time_aligned done ===")
