# ============================================================
# 03_adjusted_ror.R
# Multivariable logistic-regression adjusted ROR (aROR)
#   - Unit of analysis: report; outcome = whether the report contains
#     the event of interest (SOC or priority PT)
#   - Exposure: drug (faricimab as reference)
#   - Covariates: age group, sex, report year, region (US/non-US),
#     reporter type, indication group, concomitant-medication classes
#     (anticoagulant / antiplatelet / steroid / statin / immunosuppressant)
#   - Outputs aROR with 95% CI alongside the unadjusted ROR
#
# Prerequisite: 01_faers_prepare.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(broom); library(purrr); library(ggplot2)
})
source("R/00_utils.R")

cohort_reports <- readRDS("workdir/cohort_reports.rds")
pairs_soc <- readRDS("workdir/pairs_soc.rds")
pairs_pt  <- readRDS("workdir/pairs_pt.rds")
drug_all  <- readRDS("workdir/drug_all.rds")
assert_fresh_workdir()   # errors out on a stale workdir / wrong folder

dir.create("results/adjusted", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Concomitant-medication flags (other drugs in the same report with
#    role C/SS/I; extend the lists as needed)
# ------------------------------------------------------------
CONMED <- list(
  con_anticoag = paste("WARFARIN|COUMADIN|APIXABAN|ELIQUIS|RIVAROXABAN|XARELTO|",
                       "DABIGATRAN|PRADAXA|EDOXABAN|HEPARIN|ENOXAPARIN|LOVENOX|FONDAPARINUX"),
  con_antiplatelet = paste("ASPIRIN|CLOPIDOGREL|PLAVIX|PRASUGREL|TICAGRELOR|BRILINTA|",
                           "DIPYRIDAMOLE|CILOSTAZOL"),
  con_steroid  = "PREDNISONE|PREDNISOLONE|METHYLPREDNISOLONE|DEXAMETHASONE|TRIAMCINOLONE",
  con_statin   = "ATORVASTATIN|ROSUVASTATIN|SIMVASTATIN|PRAVASTATIN|LOVASTATIN",
  con_immunosup= "METHOTREXATE|AZATHIOPRINE|MYCOPHENOLATE|CYCLOSPORINE|TACROLIMUS|ADALIMUMAB|INFLIXIMAB"
)

conmed_flags <- drug_all %>%
  filter(role_cod %in% c("C", "I", "SS")) %>%
  # Vectorise the matching first, then summarise — never call str_detect
  # row-by-row inside summarise.
  mutate(nm = toupper(coalesce(drugname, "")),
         m_anticoag     = str_detect(nm, CONMED$con_anticoag),
         m_antiplatelet = str_detect(nm, CONMED$con_antiplatelet),
         m_steroid      = str_detect(nm, CONMED$con_steroid),
         m_statin       = str_detect(nm, CONMED$con_statin),
         m_immunosup    = str_detect(nm, CONMED$con_immunosup)) %>%
  summarise(
    con_anticoag     = any(m_anticoag),
    con_antiplatelet = any(m_antiplatelet),
    con_steroid      = any(m_steroid),
    con_statin       = any(m_statin),
    con_immunosup    = any(m_immunosup),
    .by = primaryid
  )

analysis_df <- cohort_reports %>%
  left_join(conmed_flags, by = "primaryid") %>%
  mutate(across(starts_with("con_"), ~coalesce(.x, FALSE)),
         drug = factor(drug, levels = c("faricimab", "aflibercept", "ranibizumab")),
         age_grp = factor(age_grp), sex = factor(sex),
         region = factor(region), reporter_type = factor(reporter_type),
         indication = factor(coalesce(indication, "other/unknown")))

# ------------------------------------------------------------
# 2. Dynamic formula: keep only usable covariates (a covariate that is
#    entirely missing, has a single level, or any level with <5 cases
#    is dropped with a message)
# ------------------------------------------------------------
ALL_COVARS <- c("age_grp", "sex", "region", "reporter_type",
                "report_year", "indication",
                "con_anticoag", "con_antiplatelet", "con_steroid",
                "con_statin", "con_immunosup")

build_formula <- function(dat) {
  usable <- ALL_COVARS[vapply(ALL_COVARS, function(v) {
    x <- dat[[v]]
    if (is.null(x)) return(FALSE)
    x <- x[!is.na(x)]
    if (length(unique(x)) < 2) return(FALSE)
    min(table(x)) >= 5
  }, logical(1))]
  dropped <- setdiff(ALL_COVARS, usable)
  if (length(dropped) > 0) {
    message("  Covariates dropped (insufficient data): ", paste(dropped, collapse = ", "))
  }
  rhs <- paste(c("drug", usable), collapse = " + ")
  as.formula(paste("y ~", rhs))
}

# ------------------------------------------------------------
# 3. Adjusted model for a single event
# ------------------------------------------------------------
fit_one <- function(term_name, pairs, term_col = "event") {
  ev_ids <- pairs %>%
    filter(.data[[term_col]] == term_name) %>%
    distinct(case_id) %>% pull(case_id)

  dat <- analysis_df %>%
    mutate(y = as.integer(primaryid %in% ev_ids))

  # Skip events with too few cases
  if (sum(dat$y) < 10) return(NULL)

  fml <- build_formula(dat)

  fit <- tryCatch(
    glm(fml, data = dat, family = binomial()),
    error = function(e) { message("  glm failed (", term_name, "): ", e$message); NULL }
  )
  if (is.null(fit)) return(NULL)

  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(str_detect(term, "^drug")) %>%
    transmute(
      event = term_name,
      contrast = str_remove(term, "^drug"),
      aROR = estimate, aROR025 = conf.low, aROR975 = conf.high,
      p_adj_model = p.value,
      n_event = sum(dat$y)
    )
}

# Unadjusted ROR comparator (head-to-head, reference = faricimab)
fit_unadj <- function(term_name, pairs, drugA, term_col = "event") {
  idsA <- analysis_df %>% filter(drug == drugA) %>% pull(primaryid)
  idsF <- analysis_df %>% filter(drug == "faricimab") %>% pull(primaryid)
  ev   <- pairs %>% filter(.data[[term_col]] == term_name) %>%
    distinct(case_id) %>% pull(case_id)
  a <- sum(idsA %in% ev); b <- length(idsA) - a
  c <- sum(idsF %in% ev); d <- length(idsF) - c
  r <- calc_ror(a, b, c, d)
  tibble(event = term_name, contrast = drugA,
         ROR = r$ROR, ROR025 = r$ROR025, ROR975 = r$ROR975,
         n_event = a + c)
}

# ------------------------------------------------------------
# 3. SOC-level adjustment (systemic SOCs)
# ------------------------------------------------------------
systemic_soc <- c("CARDIAC DISORDERS", "VASCULAR DISORDERS",
                  "NERVOUS SYSTEM DISORDERS", "BLOOD AND LYMPHATIC SYSTEM DISORDERS",
                  "RENAL AND URINARY DISORDERS", "HEPATOBILIARY DISORDERS",
                  "METABOLISM AND NUTRITION DISORDERS", "PSYCHIATRIC DISORDERS",
                  "INFECTIONS AND INFESTATIONS", "IMMUNE SYSTEM DISORDERS",
                  "RESPIRATORY, THORACIC AND MEDIASTINAL DISORDERS",
                  "GASTROINTESTINAL DISORDERS", "INVESTIGATIONS",
                  "NEOPLASMS BENIGN, MALIGNANT AND UNSPECIFIED (INCL CYSTS AND POLYPS)")

adj_soc <- map_dfr(systemic_soc, ~fit_one(.x, pairs_soc))

unadj_soc <- map_dfr(systemic_soc, function(s) {
  map_dfr(c("aflibercept", "ranibizumab"),
          ~fit_unadj(s, pairs_soc, .x))
})

cmp_soc <- adj_soc %>%
  left_join(unadj_soc %>% dplyr::select(-n_event), by = c("event", "contrast")) %>%
  mutate(level = "SOC")

write_tab(cmp_soc, "results/adjusted/soc_adjusted_vs_unadjusted.csv")

# ------------------------------------------------------------
# 4. Priority-PT adjustment (clinically relevant systemic events,
#    extendable; key ocular events included as comparators)
# ------------------------------------------------------------
priority_pts <- c(
  "MYOCARDIAL INFARCTION", "ACUTE MYOCARDIAL INFARCTION",
  "CEREBROVASCULAR ACCIDENT", "CEREBRAL INFARCTION", "ISCHAEMIC STROKE",
  "TRANSIENT ISCHAEMIC ATTACK",
  "DEATH", "SUDDEN DEATH", "SUDDEN CARDIAC DEATH",
  "ACUTE KIDNEY INJURY", "RENAL FAILURE", "NEPHROTIC SYNDROME", "PROTEINURIA",
  "HEPATIC FUNCTION ABNORMAL", "ACUTE HEPATIC FAILURE", "DRUG-INDUCED LIVER INJURY",
  "THROMBOCYTOPENIA", "ANAEMIA", "AGRANULOCYTOSIS", "PANCYTOPENIA",
  "HYPERTENSION", "ARTERIAL THROMBOSIS", "DEEP VEIN THROMBOSIS", "PULMONARY EMBOLISM",
  "ARTERIAL THROMBOEMBOLISM", "EMBOLISM",
  "DEPRESSION", "ANXIETY", "SUICIDAL IDEATION",
  "PERIPHERAL NEUROPATHY", "POLYNEUROPATHY",
  # key ocular events (comparators)
  "ENDOPHTHALMITIS", "RETINAL VASCULITIS", "UVEITIS", "IRIDOCYCLITIS",
  "VITREOUS HAEMORRHAGE", "RETINAL HAEMORRHAGE", "INTRAOCULAR PRESSURE INCREASED"
)

adj_pt <- map_dfr(priority_pts, ~fit_one(.x, pairs_pt))
unadj_pt <- map_dfr(priority_pts, function(p) {
  map_dfr(c("aflibercept", "ranibizumab"), ~fit_unadj(p, pairs_pt, .x))
})

cmp_pt <- adj_pt %>%
  left_join(unadj_pt %>% dplyr::select(-n_event), by = c("event", "contrast")) %>%
  mutate(level = "PT")

write_tab(cmp_pt, "results/adjusted/pt_adjusted_vs_unadjusted.csv")

# ------------------------------------------------------------
# 5. Adjusted vs unadjusted comparison plots
# ------------------------------------------------------------
plot_adj_compare <- function(df, title) {
  d <- df %>%
    dplyr::select(event, contrast, ROR, ROR025, ROR975, aROR, aROR025, aROR975) %>%
    pivot_longer(cols = c(ROR, aROR), names_to = "model", values_to = "est") %>%
    mutate(
      lo = if_else(model == "ROR", ROR025, aROR025),
      hi = if_else(model == "ROR", ROR975, aROR975),
      model = if_else(model == "ROR", "unadjusted", "adjusted")
    ) %>%
    filter(is.finite(est), lo > 0)
  ggplot(d, aes(x = est, y = event, colour = model, shape = model)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15,
                   position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    scale_x_log10() +
    facet_wrap(~contrast) +
    labs(x = "ROR (95% CI, log scale)", y = NULL, title = title) +
    theme_bw(base_size = 10)
}

ggsave("figures/adjusted_vs_unadjusted_soc.pdf",
       plot_adj_compare(cmp_soc, "SOC-level: adjusted vs unadjusted ROR"),
       width = 11, height = 7)
ggsave("figures/adjusted_vs_unadjusted_pt.pdf",
       plot_adj_compare(cmp_pt, "PT-level: adjusted vs unadjusted ROR"),
       width = 11, height = 9)

message("=== 03_adjusted_ror done ===")
