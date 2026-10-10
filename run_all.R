# ============================================================
# run_all.R
# Runs the full analysis pipeline (scripts 01-10) in order.
# Run from the repository root: source("run_all.R")
#
# One-time caching scripts (00a / 00a2 / 00a3) are NOT included here;
# run them once beforehand to build the quarterly .rds cache.
# ============================================================

steps <- c(
  "R/01_faers_prepare.R",
  "R/02_faers_signals.R",
  "R/03_adjusted_ror.R",
  "R/04_tto_weibull.R",
  "R/05_sensitivity_time_aligned.R",   # needs historical quarters in the cache
  "R/06_sensitivity_masking.R",
  "R/07_jader_analysis.R",             # needs data/jader/
  "R/08_eudravigilance_ror.R",         # needs data/eudravigilance_soc_counts.csv
  "R/09_cross_db_validation.R",
  "R/10_figure_quarterly_trend.R"
)

t0 <- Sys.time()
for (s in steps) {
  message("\n==============================")
  message("Running ", s)
  message("==============================")
  tryCatch(
    source(s, local = new.env()),
    error = function(e) {
      message("!! ", s, " failed: ", e$message)
      message("   (If this is a sensitivity/validation script, check its ",
              "input prerequisites in the README and re-run.)")
    }
  )
}
message("\n=== Pipeline finished in ",
        round(difftime(Sys.time(), t0, units = "mins"), 1), " minutes ===")
