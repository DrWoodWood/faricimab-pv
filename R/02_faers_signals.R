# ============================================================
# 02_faers_signals.R
# FAERS primary analysis:
#   A. Each drug vs full-database background (SOC / PT levels, FDR)
#   B. Head-to-head comparisons (aflibercept vs faricimab,
#      ranibizumab vs faricimab, ranibizumab vs aflibercept)
#   C. Serious-outcome (death etc.) disproportionality
#   D. DME (EMA designated medical events) screen
#   E. Figures (forest plots / volcano plots)
#
# Prerequisite: run 01_faers_prepare.R first
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(ggrepel)
})
source("R/00_utils.R")

cohort_reports <- readRDS("workdir/cohort_reports.rds")
pairs_pt  <- readRDS("workdir/pairs_pt.rds")
pairs_soc <- readRDS("workdir/pairs_soc.rds")
bg_pt     <- readRDS("workdir/bg_pt.rds")
bg_soc    <- readRDS("workdir/bg_soc.rds")
assert_fresh_workdir()   # errors out on a stale workdir / wrong folder

# Denominator: full cohort size; reports left with no events after PT
# exclusion stay in the denominator (standard 2x2 definition).
cohort_n <- readRDS("workdir/cohort_n.rds")
bg_n     <- readRDS("workdir/bg_n.rds")
n_of <- function(dg) cohort_n$n[cohort_n$drug_std == dg]

dir.create("results/faers", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE)

DRUGS <- c("faricimab", "aflibercept", "ranibizumab")

# ------------------------------------------------------------
# A. Each drug vs full-database background
# ------------------------------------------------------------
run_vs_bg <- function(pairs, bg, level = c("PT", "SOC")) {
  level <- match.arg(level)
  map_dfr(DRUGS, function(dg) {
    tgt <- pairs %>% filter(drug == dg) %>% distinct(case_id, event)
    disprop_all(tgt, bg, n_min = 3,
                n_target_total = n_of(dg), n_bg_total = bg_n) %>%
      mutate(drug = dg, .before = 1)
  }) %>%
    flag_signal(min_n = 3) %>%     # BH FDR within drug x level family
    add_fisher()
}

sig_soc_vsbg <- run_vs_bg(pairs_soc, bg_soc, "SOC")
sig_pt_vsbg  <- run_vs_bg(pairs_pt,  bg_pt,  "PT")

message("Signals vs background: PT ", nrow(sig_pt_vsbg), " rows / SOC ",
        nrow(sig_soc_vsbg), " rows")

write_tab(sig_soc_vsbg, "results/faers/A1_soc_vs_fulldb.csv")
write_tab(sig_pt_vsbg,  "results/faers/A2_pt_vs_fulldb.csv")

# ------------------------------------------------------------
# B. Head-to-head (each drug serves as the other's background)
# ------------------------------------------------------------
cmp <- list(c("aflibercept", "faricimab"),
            c("ranibizumab", "faricimab"),
            c("ranibizumab", "aflibercept"))

h2h_soc <- map_dfr(cmp, ~head_to_head(pairs_soc, .x[1], .x[2],
                                      n_A = n_of(.x[1]), n_B = n_of(.x[2]))) %>%
  group_by(comparison) %>%
  group_modify(~flag_signal(.x, min_n = 3)) %>% ungroup() %>%
  add_fisher()

h2h_pt <- map_dfr(cmp, ~head_to_head(pairs_pt, .x[1], .x[2],
                                     n_A = n_of(.x[1]), n_B = n_of(.x[2]))) %>%
  group_by(comparison) %>%
  group_modify(~flag_signal(.x, min_n = 3)) %>% ungroup() %>%
  add_fisher()

write_tab(h2h_soc, "results/faers/B1_soc_headtohead.csv")
write_tab(h2h_pt,  "results/faers/B2_pt_headtohead.csv")

# ------------------------------------------------------------
# C. Serious-outcome disproportionality (death / life-threatening /
#    hospitalisation / disability). Requires outc; skipped if the cache
#    lacks it (01 sets the columns to NA).
# ------------------------------------------------------------
if (all(is.na(cohort_reports$death))) {
  message("!! Serious-outcome data unavailable (no outc table); skipping section C")
} else {
serious_h2h <- map_dfr(cmp, function(cp) {
  rA <- cohort_reports %>% filter(drug == cp[1])
  rB <- cohort_reports %>% filter(drug == cp[2])
  map_dfr(c("death", "life_threat", "hosp", "disability"), function(oc) {
    a <- sum(rA[[oc]]); b <- nrow(rA) - a
    c <- sum(rB[[oc]]); d <- nrow(rB) - c
    ror <- calc_ror(a, b, c, d); prr <- calc_prr(a, b, c, d); ic <- calc_ic(a, b, c, d)
    tibble(comparison = paste(cp[1], "vs", cp[2]),
           outcome = oc, a = a, b = b, c = c, d = d,
           ROR = ror$ROR, ROR025 = ror$ROR025, ROR975 = ror$ROR975,
           PRR = prr$PRR, CHISQ = prr$CHISQ, P = prr$PCHISQ,
           IC = ic$IC, IC025 = ic$IC025)
  })
}) %>%
  group_by(comparison) %>%
  mutate(P_FDR = p.adjust(P, method = "BH")) %>%
  ungroup()

write_tab(serious_h2h, "results/faers/C1_serious_outcomes_headtohead.csv")
}

# ------------------------------------------------------------
# D. DME (EMA designated medical events) screen.
#    data/dme_list.csv is a subset example; replace with the full EMA
#    DME list for a formal analysis.
# ------------------------------------------------------------
if (file.exists("data/dme_list.csv")) {
  dme <- read_csv("data/dme_list.csv", show_col_types = FALSE)$pt %>% toupper()
  dme_signals <- sig_pt_vsbg %>% filter(event %in% dme, signal)
  write_tab(dme_signals, "results/faers/D1_dme_signals.csv")
}

# ------------------------------------------------------------
# E. Figures
# ------------------------------------------------------------
# Volcano: PT-level (vs full database) for each drug
for (dg in DRUGS) {
  d <- sig_pt_vsbg %>% filter(drug == dg)
  ggsave(paste0("figures/volcano_", dg, ".pdf"),
         plot_volcano(d, title = paste0("FAERS PT-level signals: ", dg)),
         width = 9, height = 7)
}

# Forest: systemic SOC head-to-head (ocular/procedural SOCs excluded)
systemic_soc <- c("CARDIAC DISORDERS", "VASCULAR DISORDERS",
                  "NERVOUS SYSTEM DISORDERS", "BLOOD AND LYMPHATIC SYSTEM DISORDERS",
                  "RENAL AND URINARY DISORDERS", "HEPATOBILIARY DISORDERS",
                  "METABOLISM AND NUTRITION DISORDERS", "PSYCHIATRIC DISORDERS",
                  "INFECTIONS AND INFESTATIONS", "NEOPLASMS BENIGN, MALIGNANT AND UNSPECIFIED (INCL CYSTS AND POLYPS)",
                  "RESPIRATORY, THORACIC AND MEDIASTINAL DISORDERS",
                  "GASTROINTESTINAL DISORDERS", "IMMUNE SYSTEM DISORDERS")

for (cp in cmp) {
  cname <- paste(cp[1], "vs", cp[2])
  d <- h2h_soc %>%
    filter(comparison == cname, event %in% systemic_soc)
  ggsave(paste0("figures/forest_soc_", gsub(" ", "_", cname), ".pdf"),
         plot_forest(d, "event", title = paste0("SOC-level head-to-head: ", cname)),
         width = 8.5, height = 6)
}

# EBGM (optional; runs only if openEBGM is installed). Restricted to the
# three-drug cohort for a within-class comparison and to bound compute.
eb <- tryCatch({
  pairs3 <- pairs_pt %>%
    transmute(id = case_id, drug = drug, event = event)
  run_ebgm(pairs3)
}, error = function(e) { message("EBGM skipped: ", e$message); NULL })

if (!is.null(eb)) write_tab(eb, "results/faers/A3_ebgm_within_class.csv")

message("=== 02_faers_signals done ===")
