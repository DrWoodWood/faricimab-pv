# Faricimab pharmacovigilance: triangulated disproportionality analysis

Reproducible R pipeline for a pharmacovigilance study of **faricimab (Vabysmo)**
versus **aflibercept (Eylea)** and **ranibizumab (Lucentis)** across three
spontaneous-reporting systems:

- **FAERS** (FDA Adverse Event Reporting System) — primary analysis
- **JADER** (Japanese PMDA database) — external validation
- **EudraVigilance** (EMA, public aggregate interface) — SOC-level validation

## What the pipeline does

1. **Disproportionality analysis** — ROR, PRR and IC (with Haldane 0.5
   correction and Norén shrinkage), Benjamini–Hochberg FDR control, at both
   MedDRA SOC and PT levels, for each drug vs the full database and in
   pairwise head-to-head comparisons.
2. **Adjusted ROR (aROR)** — per-event multivariable logistic regression
   (report level), adjusting for age group, sex, report year, region,
   reporter type, indication group and concomitant-medication classes.
3. **Time-to-onset (TTO)** — Weibull shape-parameter classification
   (early / random / wear-out failure) and Kruskal–Wallis between-drug
   comparisons.
4. **Sensitivity analyses** — post-marketing first-3-years era-aligned
   comparison (Weber effect) and masking/competition-bias assessment
   (masked vs unmasked background, Spearman consistency).
5. **Cross-database triangulation** — evidence grading of every FAERS
   head-to-head signal against JADER and EudraVigilance
   (Confirmed / Consistent / Unconfirmed / Opposite).

## Repository structure

```
├── R/
│   ├── 00_utils.R                      # shared helpers: ROR/PRR/IC, signal flags, loaders, plots
│   ├── 00a_build_rds_from_ascii.R      # one-time: parse FAERS ASCII quarterly packs -> .rds cache
│   ├── 00a2_build_rds_with_faers_pkg.R # alternative one-time caching via the faers package
│   ├── 00a3_rebuild_legacy_quarters.R  # one-time: re-parse legacy AERS quarters (2004Q1-2012Q3)
│   ├── 01_faers_prepare.R              # cohort building, drug identification, backgrounds, TTO
│   ├── 02_faers_signals.R              # FAERS disproportionality: vs background + head-to-head
│   ├── 03_adjusted_ror.R               # adjusted ROR (logistic regression)
│   ├── 04_tto_weibull.R                # TTO + Weibull shape analysis
│   ├── 05_sensitivity_time_aligned.R   # era-aligned (first 3 years post-approval) sensitivity
│   ├── 06_sensitivity_masking.R        # masking / competition-bias sensitivity
│   ├── 07_jader_analysis.R             # JADER validation
│   ├── 08_eudravigilance_ror.R         # EudraVigilance SOC-level validation
│   ├── 09_cross_db_validation.R        # cross-database evidence grading
│   └── 10_figure_quarterly_trend.R     # Figure 1: quarterly reporting trends
├── data/                               # input templates (see Data below)
├── install_packages.R                  # installs all required R packages
├── run_all.R                           # runs the full pipeline (01-10) in order
├── .gitignore
└── LICENSE
```

## Requirements

- R >= 4.2 (developed and tested on macOS; works on Linux/Windows)
- Packages: run `source("install_packages.R")` once. Core dependencies are
  `dplyr tidyr readr stringr lubridate purrr ggplot2 ggrepel broom
  fitdistrplus data.table`; `openEBGM` and the `faers` package are optional.

## Data (not included)

Raw data are **not** redistributed in this repository:

| Data | Source | Licence note |
|---|---|---|
| FAERS quarterly ASCII packs | [FDA website](https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html) (2004Q1 onwards) | public |
| MedDRA 27.1 (English) | [MedDRA MSSO](https://www.meddra.org/) | **licence required** — PT-to-SOC maps are built locally and must not be committed |
| JADER yearly CSVs | [PMDA](https://www.pmda.go.jp/safety/info-services/drugs/adr-info/suspected-adr/0003.html) | public; PTs are MedDRA/J (Japanese) |
| EudraVigilance SOC counts | [adrreports.eu](https://www.adrreports.eu/en/search_subst.html) | public aggregate counts, transcribed manually |

The `data/` folder ships **templates** (`*.template.csv`). Copy each template
to the filename without `.template` and fill it in (the real filenames are
git-ignored):

- `data/excluded_pts.csv` — non-AE MedDRA PTs excluded from all analyses
  (one PT per line, column `pt`)
- `data/dme_list.csv` — EMA designated medical events (column `pt`)
- `data/meddra_j_en_map.csv` — JADER Japanese-PT to English-PT mapping
  (columns `pt_j, pt_en`)
- `data/eudravigilance_soc_counts.csv` — columns
  `drug, soc, n_cases_soc, total_cases_drug`

## Configuration

Paths are read from environment variables with project-relative defaults,
so no code changes are needed if you follow the layout above. To override
(e.g. point to an existing FAERS cache), add to `~/.Renviron`:

```r
FAERS_DIR=/path/to/rds_cache_full        # quarterly .rds cache
FAERS_ASCII_DIR=/path/to/faers_ascii     # FDA ASCII packs (00a)
FAERS_ZIP_DIR=/path/to/zip_cache         # quarterly zips (00a2)
MEDDRA_ZIP=/path/to/MedDRA_27_1_English.zip
```

## Quick start

```r
# 0. One-time: build the quarterly .rds cache from FDA ASCII packs
source("R/00a_build_rds_from_ascii.R")
#    (legacy quarters 2004Q1-2012Q3 need the dedicated parser:)
source("R/00a3_rebuild_legacy_quarters.R")

# 1. Full analysis pipeline
source("run_all.R")
```

`run_all.R` executes scripts 01-10 in order. Scripts 05-09 are skipped
automatically with a message if their inputs (historical quarters, JADER
files, EudraVigilance counts) are unavailable.

## Outputs

- `results/` — CSV tables (disproportionality, adjusted ROR, TTO/Weibull,
  sensitivity analyses, cross-database grading)
- `figures/` — PDF figures (volcano, forest, ECDF, masking consistency,
  era-aligned, quarterly trend, cross-database forest)
- `workdir/` — intermediate `.rds` objects (git-ignored)

## Methodological notes

- Signal definition: lower 95% CI of ROR > 1 with event count >= 3,
  BH-FDR controlled within each drug x level family.
- Zero cells handled with the Haldane-Anscombe 0.5 correction.
- Head-to-head direction is reported bidirectionally: a significant ROR > 1
  favours the first-named drug, ROR < 1 favours the comparator.
- TTO = 0 (same-day onset) is analysed as 0.5 day, per FAERS TTO convention.
- All scripts source `R/00_utils.R`; run everything from the repository root.

## License

Code is released under the MIT License (see `LICENSE`). Data remain subject
to their own terms — in particular, MedDRA content requires a valid licence
and must not be redistributed.
