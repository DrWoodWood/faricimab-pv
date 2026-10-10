# ============================================================
# 04_tto_weibull.R
# Time-to-onset (TTO) analysis + Weibull shape-parameter testing
#   - Outcomes: all AEs per drug, plus key event groups
#     (IOI / haemorrhage / arterial thromboembolism / death)
#   - Weibull shape parameter beta:
#       upper 95% CI < 1  -> early failure
#       95% CI covers 1   -> random failure
#       lower 95% CI > 1  -> wear-out failure
#   - Between-drug TTO comparison: Kruskal-Wallis + ECDF plots
#
# Prerequisite: 01_faers_prepare.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(purrr)
  library(fitdistrplus); library(ggplot2)
})
source("R/00_utils.R")

if (!file.exists("workdir/tto.rds")) {
  stop("workdir/tto.rds not found - 01_faers_prepare.R did not produce TTO data ",
       "(cache lacks the ther table or demo lacks event_dt). ",
       "Fix the cache and re-run 01 before running this script.", call. = FALSE)
}
tto <- readRDS("workdir/tto.rds")
assert_fresh_workdir()   # errors out on a stale workdir
dir.create("results/tto", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 0. Quick sanity check of the TTO data
# ------------------------------------------------------------
message("--- TTO data summary ---")
message("rows: ", nrow(tto))
print(table(tto$drug_std, useNA = "ifany"))
print(summary(tto$tto_days))
message("tto_days == 0 (same-day onset): ",
        round(mean(tto$tto_days == 0, na.rm = TRUE) * 100, 1), "%")

# ------------------------------------------------------------
# 1. Event-group definitions (MedDRA PTs, extendable)
# ------------------------------------------------------------
EVENT_GROUPS <- list(
  all_AE = NULL,   # NULL = all events
  IOI = c("ENDOPHTHALMITIS", "UVEITIS", "IRITIS", "IRIDOCYCLITIS", "VITRITIS",
          "CHORIORETINITIS", "ANTERIOR CHAMBER INFLAMMATION", "HYPOPYON",
          "PANOPHTHALMITIS", "PSEUDOENDOPHTHALMITIS", "INTRAOCULAR INFLAMMATION",
          "RETINAL VASCULITIS"),
  ocular_haemorrhage = c("VITREOUS HAEMORRHAGE", "RETINAL HAEMORRHAGE",
                         "CHOROIDAL HAEMORRHAGE", "SUBRETINAL HAEMATOMA",
                         "HYPHEMA", "CONJUNCTIVAL HAEMORRHAGE"),
  arterial_thromboembolism = c("MYOCARDIAL INFARCTION", "ACUTE MYOCARDIAL INFARCTION",
                               "CEREBROVASCULAR ACCIDENT", "CEREBRAL INFARCTION",
                               "ISCHAEMIC STROKE", "TRANSIENT ISCHAEMIC ATTACK",
                               "ARTERIAL THROMBOSIS", "ARTERIAL THROMBOEMBOLISM"),
  death = "DEATH"
)

# ------------------------------------------------------------
# 2. Weibull fit for one group.
#    TTO = 0 (same-day onset) is set to 0.5 day per FAERS TTO
#    convention. Always returns one row: on failure the shape columns
#    are NA and fit_note explains why.
# ------------------------------------------------------------
fit_weibull <- function(days) {
  days <- days[is.finite(days) & days >= 0]
  n_total <- length(days)
  n_zero  <- sum(days == 0)

  base <- list(
    n = n_total, n_zero = n_zero,
    median_tto = if (n_total > 0) median(days) else NA_real_,
    q1 = if (n_total > 0) unname(quantile(days, 0.25)) else NA_real_,
    q3 = if (n_total > 0) unname(quantile(days, 0.75)) else NA_real_,
    shape = NA_real_, shape_lo = NA_real_, shape_hi = NA_real_,
    scale = NA_real_, fail_type = NA_character_, fit_note = NA_character_
  )
  if (n_total == 0)   { base$fit_note <- "no valid records";        return(base) }
  if (n_total < 30)   { base$fit_note <- "n<30, not fitted";        return(base) }

  days_fit <- ifelse(days == 0, 0.5, days)   # same-day onset as 0.5 day
  fit <- tryCatch(fitdist(days_fit, "weibull"), error = function(e) e)
  if (inherits(fit, "error")) {
    base$fit_note <- paste0("fitdist failed: ", conditionMessage(fit))
    return(base)
  }
  ci <- tryCatch(confint(fit), error = function(e) NULL)
  if (is.null(ci)) {
    base$fit_note <- "confint failed"
    return(base)
  }
  base$shape    <- unname(fit$estimate["shape"])
  base$shape_lo <- unname(ci["shape", 1])
  base$shape_hi <- unname(ci["shape", 2])
  base$scale    <- unname(fit$estimate["scale"])
  base$fail_type <- case_when(
    base$shape_hi < 1 ~ "early failure",
    base$shape_lo > 1 ~ "wear-out failure",
    TRUE              ~ "random failure"
  )
  base$fit_note <- "OK"
  base
}

# ------------------------------------------------------------
# 3. Loop over drug x event group
# ------------------------------------------------------------
results <- map_dfr(names(EVENT_GROUPS), function(g) {
  pts <- EVENT_GROUPS[[g]]
  map_dfr(unique(tto$drug_std), function(dg) {
    d <- tto %>% filter(drug_std == dg)
    if (!is.null(pts)) d <- d %>% filter(pt %in% pts)
    # keep the earliest onset per case per group
    d <- d %>% group_by(case_id) %>% slice_min(tto_days, n = 1) %>% ungroup()
    r <- fit_weibull(d$tto_days)
    message(sprintf("  [%s x %s] n=%d (same-day %d) -> %s",
                    dg, g, r$n, r$n_zero, r$fit_note))
    tibble(drug = dg, event_group = g, !!!r)
  })
})

write_tab(results, "results/tto/weibull_by_drug_group.csv")

# ------------------------------------------------------------
# 4. Between-drug TTO comparison (Kruskal-Wallis)
# ------------------------------------------------------------
kw <- map_dfr(names(EVENT_GROUPS), function(g) {
  pts <- EVENT_GROUPS[[g]]
  d <- tto
  if (!is.null(pts)) d <- d %>% filter(pt %in% pts)
  d <- d %>% group_by(case_id) %>% slice_min(tto_days, n = 1) %>% ungroup()
  if (n_distinct(d$drug_std) < 2 || nrow(d) < 50) return(NULL)
  kt <- kruskal.test(tto_days ~ drug_std, data = d)
  tibble(event_group = g, statistic = unname(kt$statistic),
         p_value = kt$p.value, n = nrow(d))
})
write_tab(kw, "results/tto/kruskal_wallis.csv")

# ------------------------------------------------------------
# 5. Empirical cumulative distribution (ECDF) plots
# ------------------------------------------------------------
for (g in c("all_AE", "IOI", "arterial_thromboembolism", "death")) {
  pts <- EVENT_GROUPS[[g]]
  d <- tto
  if (!is.null(pts)) d <- d %>% filter(pt %in% pts)
  d <- d %>% group_by(case_id) %>% slice_min(tto_days, n = 1) %>% ungroup()
  if (nrow(d) < 30) next
  p <- ggplot(d, aes(x = tto_days, colour = drug_std)) +
    stat_ecdf(linewidth = 0.9) +
    coord_cartesian(xlim = c(0, 365)) +
    labs(x = "Time to onset (days)", y = "Cumulative proportion",
         colour = "Drug", title = paste0("TTO cumulative curve: ", g)) +
    theme_bw(base_size = 12)
  ggsave(paste0("figures/tto_ecdf_", g, ".pdf"), p, width = 7, height = 5)
}

message("=== 04_tto_weibull done ===")
