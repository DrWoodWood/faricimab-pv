# ============================================================
# install_packages.R
# One-time setup: install all R packages required by the pipeline.
# Run: source("install_packages.R")
# ============================================================

core_pkgs <- c(
  "dplyr", "tidyr", "readr", "stringr", "lubridate", "purrr",
  "ggplot2", "ggrepel", "broom", "fitdistrplus", "data.table"
)

optional_pkgs <- c(
  "openEBGM",   # empirical Bayes (EBGM) disproportionality; 02 runs without it
  "faers"       # alternative FAERS download/caching route (00a2); install from
                # GitHub if CRAN is unavailable: remotes::install_github(...)
)

install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    message("Installing: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org")
  } else {
    message("All packages already installed.")
  }
}

install_if_missing(core_pkgs)

# Optional packages: attempt, but never fail the setup over them
for (p in optional_pkgs) {
  ok <- tryCatch({
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org")
    }
    TRUE
  }, error = function(e) {
    message("Optional package '", p, "' could not be installed (skipping): ", e$message)
    FALSE
  })
}

message("=== Package setup complete ===")
