# ============================================================
# 10_figure_quarterly_trend.R
# Figure 1: quarterly reporting trends for the three drugs in FAERS
#   Data : workdir/cohort_reports.rds (from 01; has fda_date/drug)
#   Style: faricimab blue solid circle / aflibercept orange solid square /
#          ranibizumab green solid triangle
# Prerequisite: 01_faers_prepare.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(lubridate); library(readr)
})
source("R/00_utils.R")

if (!file.exists("workdir/cohort_reports.rds")) {
  stop("workdir/cohort_reports.rds not found - run 01_faers_prepare.R first", call. = FALSE)
}
cr <- readRDS("workdir/cohort_reports.rds")
assert_fresh_workdir()
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Aggregate by quarter
# ------------------------------------------------------------
qtr <- cr %>%
  filter(!is.na(fda_date)) %>%
  mutate(quarter = as.Date(paste0(year(fda_date), "-",
                                  sprintf("%02d", (quarter(fda_date) - 1) * 3 + 1), "-01"))) %>%
  count(drug, quarter, name = "n")

# Fill quarters with zero reports so lines drop to 0 instead of breaking
full_grid <- expand.grid(
  drug = unique(qtr$drug),
  quarter  = seq(min(qtr$quarter), max(qtr$quarter), by = "3 months")
)
qtr <- full_grid %>%
  left_join(qtr, by = c("drug", "quarter")) %>%
  mutate(n = coalesce(n, 0L),
         drug = factor(drug,
                           levels = c("faricimab", "aflibercept", "ranibizumab")))

message("Quarter range: ", min(qtr$quarter), " to ", max(qtr$quarter))
print(qtr %>% group_by(drug) %>% summarise(total = sum(n), .groups = "drop"))

# ------------------------------------------------------------
# 2. Plot (Figure 1 colours and point shapes)
# ------------------------------------------------------------
cols  <- c(faricimab = "#2166AC", aflibercept = "#E08214", ranibizumab = "#1B9E77")
shps  <- c(faricimab = 16, aflibercept = 15, ranibizumab = 17)   # circle/square/triangle

p <- ggplot(qtr, aes(x = quarter, y = n, colour = drug,
                     shape = drug, group = drug)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_shape_manual(values = shps, name = NULL) +
  scale_x_date(date_breaks = "6 months",
               labels = function(x) paste0(lubridate::year(x), " Q", lubridate::quarter(x)),
               expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Number of reports",
       title = "Quarterly trends of adverse event reports, FAERS") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank())

ggsave("figures/quarterly_trend.pdf", p, width = 9, height = 5)
message("=== done: figures/quarterly_trend.pdf ===")
