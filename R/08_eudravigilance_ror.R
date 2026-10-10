# ============================================================
# 08_eudravigilance_ror.R
# EudraVigilance aggregate-data validation (SOC level)
#
# Data acquisition (public, no registration needed):
#   1. Open https://www.adrreports.eu/en/search_subst.html
#   2. Search each active substance: FARICIMAB / AFLIBERCEPT / RANIBIZUMAB
#   3. On the results page, "Number of Individual Cases by Reaction
#      Groups" tab: transcribe/export the case count per SOC plus the
#      drug's total case count
#   4. Fill data/eudravigilance_soc_counts.csv (template in data/)
#
# Limitation: the public interface only exposes SOC-level aggregate
# counts; PT-level detail requires an EVDAS line-listing request to the
# EMA. This source therefore serves as a third, SOC-level validator.
# ============================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(purrr); library(tidyr)
})
source("R/00_utils.R")

ev <- read_csv("data/eudravigilance_soc_counts.csv", show_col_types = FALSE)
# Expected columns: drug, soc, n_cases_soc, total_cases_drug
#   drug: faricimab / aflibercept / ranibizumab
#   soc : SOC name in English (uppercase MedDRA SOC recommended)
#   n_cases_soc: cases of that SOC for the drug (cases may span SOCs)
#   total_cases_drug: total cases for the drug

dir.create("results/eudravigilance", showWarnings = FALSE, recursive = TRUE)

ev <- ev %>% mutate(soc = toupper(soc))

cmp <- list(c("aflibercept", "faricimab"),
            c("ranibizumab", "faricimab"),
            c("ranibizumab", "aflibercept"))

ev_h2h <- map_dfr(cmp, function(cp) {
  A <- ev %>% filter(drug == cp[1])
  B <- ev %>% filter(drug == cp[2])
  totA <- unique(A$total_cases_drug); totB <- unique(B$total_cases_drug)
  full_join(A %>% dplyr::select(soc, a = n_cases_soc),
            B %>% dplyr::select(soc, c = n_cases_soc), by = "soc") %>%
    mutate(across(c(a, c), ~coalesce(.x, 0L)),
           b = totA - a, d = totB - c) %>%
    rowwise() %>%
    mutate(
      ROR = calc_ror(a, b, c, d)$ROR,
      ROR025 = calc_ror(a, b, c, d)$ROR025,
      ROR975 = calc_ror(a, b, c, d)$ROR975
    ) %>%
    ungroup() %>%
    mutate(comparison = paste(cp[1], "vs", cp[2]), .before = 1) %>%
    filter(a >= 1 | c >= 1)
})

write_tab(ev_h2h, "results/eudravigilance/soc_headtohead.csv")
message("=== 08_eudravigilance_ror done ===")
