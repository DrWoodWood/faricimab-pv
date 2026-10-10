# ============================================================
# 00_utils.R
# Shared utilities: disproportionality metrics / FDR adjustment / plotting
# Faricimab vs aflibercept vs ranibizumab pharmacovigilance project
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(lubridate)
})

# ------------------------------------------------------------
# 1. Basic 2x2 table metrics
#    a = target drug & target event
#    b = target drug & other events
#    c = background   & target event
#    d = background   & other events
# ------------------------------------------------------------

# Haldane-Anscombe 0.5 correction (applied when any cell is 0).
# Inputs must be coerced to double first: against a full-database
# background d ~ 1e7, and integer multiplication a*d overflows
# (R integer limit ~2.1e9), leaving downstream CHISQ/P/IC all NA.
haldane <- function(a, b, c, d) {
  a <- as.numeric(a); b <- as.numeric(b)
  c <- as.numeric(c); d <- as.numeric(d)
  zero <- (a == 0 | b == 0 | c == 0 | d == 0)
  a <- ifelse(zero, a + 0.5, a)
  b <- ifelse(zero, b + 0.5, b)
  c <- ifelse(zero, c + 0.5, c)
  d <- ifelse(zero, d + 0.5, d)
  list(a = a, b = b, c = c, d = d)
}

# ROR (Rothman 2004)
calc_ror <- function(a, b, c, d) {
  z <- haldane(a, b, c, d)
  ror <- (z$a * z$d) / (z$b * z$c)
  se  <- sqrt(1/z$a + 1/z$b + 1/z$c + 1/z$d)
  list(ROR = ror,
       ROR025 = exp(log(ror) - 1.96 * se),
       ROR975 = exp(log(ror) + 1.96 * se))
}

# PRR + chi-squared (Evans 2001; Yates continuity correction)
calc_prr <- function(a, b, c, d) {
  # chi-squared uses raw counts but must be double-precision,
  # otherwise a*d overflows to NA
  a <- as.numeric(a); b <- as.numeric(b)
  c <- as.numeric(c); d <- as.numeric(d)
  z <- haldane(a, b, c, d)
  prr <- (z$a / (z$a + z$b)) / (z$c / (z$c + z$d))
  se  <- sqrt(1/z$a - 1/(z$a + z$b) + 1/z$c - 1/(z$c + z$d))
  N <- a + b + c + d
  chisq <- (N * (abs(a * d - b * c) - N / 2)^2) /
    ((a + b) * (c + d) * (a + c) * (b + d))
  chisq <- ifelse(is.nan(chisq) | is.infinite(chisq), NA_real_, chisq)
  list(PRR = prr,
       PRR025 = exp(log(prr) - 1.96 * se),
       PRR975 = exp(log(prr) + 1.96 * se),
       CHISQ = chisq,
       PCHISQ = pchisq(chisq, df = 1, lower.tail = FALSE),
       # P on the log10 scale: at large N the plain P value underflows
       # to 0; the log-scale version is kept as a fallback for plotting
       LOG10P = pchisq(chisq, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10))
}

# BCPNN information component (Bate 1998; shrinkage O/E form, Norén 2006)
# Priors: alpha1 = alpha2 = 1 -> alpha = 2, alpha11 = 0.5
calc_ic <- function(a, b, c, d) {
  # double-precision again: (a+b)*(a+c) can exceed 2^31 against a
  # full-database background
  a <- as.numeric(a); b <- as.numeric(b)
  c <- as.numeric(c); d <- as.numeric(d)
  N <- a + b + c + d
  E <- (a + b) * (a + c) / N          # expected count
  alpha <- 2; alpha11 <- 0.5
  IC <- log2((a + alpha11) / (E + alpha))
  # Delta-method variance of log(O/E) under the Norén shrinkage
  # framework; always positive. Do not use the form containing
  # (E + alpha - a - alpha11): for strong signals E << a makes that
  # term negative, sqrt() returns NaN, and IC025 is lost for exactly
  # the strongest signals.
  var_ic <- (1 / log(2))^2 * (1 / (a + alpha11) + 1 / (E + alpha))
  list(IC = IC,
       IC025 = IC - 1.96 * sqrt(var_ic),
       IC975 = IC + 1.96 * sqrt(var_ic))
}

# ------------------------------------------------------------
# 2. Core engine: compute all metrics per event for one target
#    drug vs a background
#    target     : data.frame(case_id, event), already distinct
#                 (one row per case-event pair)
#    background : data.frame(case_id, event), disjoint from target
#    n_min      : minimum report count threshold (usually 3)
#    n_target_total / n_bg_total: denominators (total case counts).
#      By default inferred from the pairs; after excluding non-AE
#      PTs, pass the full cohort counts (cohort_n / bg_n saved by 01)
#      so that reports with no remaining event still stay in the
#      denominator (standard 2x2 convention).
# ------------------------------------------------------------
disprop_all <- function(target, background, n_min = 3,
                        n_target_total = NULL, n_bg_total = NULL) {
  n_t <- if (!is.null(n_target_total)) n_target_total else dplyr::n_distinct(target$case_id)
  n_b <- if (!is.null(n_bg_total)) n_bg_total else dplyr::n_distinct(background$case_id)

  tb <- target %>% count(event, name = "a")
  bb <- background %>% count(event, name = "c")

  res <- full_join(tb, bb, by = "event") %>%
    mutate(across(c(a, c), ~coalesce(.x, 0L))) %>%
    filter(a >= n_min) %>%
    mutate(
      b = n_t - a,
      d = n_b - c,
      N = a + b + c + d,
      n_target_total = n_t,
      n_bg_total = n_b
    )

  ror <- calc_ror(res$a, res$b, res$c, res$d)
  prr <- calc_prr(res$a, res$b, res$c, res$d)
  ic  <- calc_ic(res$a, res$b, res$c, res$d)

  res %>%
    mutate(
      ROR = ror$ROR, ROR025 = ror$ROR025, ROR975 = ror$ROR975,
      PRR = prr$PRR, PRR025 = prr$PRR025, PRR975 = prr$PRR975,
      CHISQ = prr$CHISQ, P = prr$PCHISQ, LOG10P = prr$LOG10P,
      IC = ic$IC, IC025 = ic$IC025, IC975 = ic$IC975
    ) %>%
    arrange(desc(ROR))
}

# ------------------------------------------------------------
# 3. Signal detection (multi-criteria, configurable)
#    Default: a >= 3, ROR025 > 1, PRR >= 2 with CHISQ >= 4, IC025 > 0
#    FDR: Benjamini-Hochberg within each analysis family
# ------------------------------------------------------------
flag_signal <- function(df,
                        min_n = 3,
                        use_prr = TRUE, use_ic = TRUE) {
  df <- df %>%
    mutate(
      sig_ror = (a >= min_n) & (ROR025 > 1),
      sig_prr = (a >= min_n) & (PRR >= 2) & (CHISQ >= 4),
      sig_ic  = (a >= min_n) & (IC025 > 0),
      signal  = sig_ror &
        (if (use_prr) sig_prr else TRUE) &
        (if (use_ic)  sig_ic  else TRUE)
    )
  # Within-family BH adjustment (caller guarantees df belongs to one
  # testing family, e.g. one database at one level)
  df <- df %>% mutate(P_FDR = p.adjust(P, method = "BH"))
  # BH-FDR on the log10 scale: BH is rank-based, so the log-scale
  # computation is equivalent and avoids P-value underflow to 0 at
  # large N (used for plotting; exact values can reach the hundreds)
  if ("LOG10P" %in% names(df)) {
    lp <- df$LOG10P                       # log10(P); smaller = more significant
    lp_fdr <- rep(NA_real_, length(lp))
    ok <- !is.na(lp)
    if (any(ok)) {
      lpv <- lp[ok]
      o  <- order(lpv)                    # most significant first
      m  <- length(lpv)
      adj <- log10(m) - log10(seq_len(m)) + lpv[o]
      adj <- rev(cummin(rev(adj)))        # BH step-up cumulative minimum
      lp_fdr[ok] <- adj[order(o)]
    }
    df$LOG10P_FDR <- lp_fdr               # log10(FDR P); negative
  }
  df
}

# Fisher exact test for small cells (more robust when any cell < 5)
add_fisher <- function(df, cell_min = 5) {
  idx <- which(df$a < cell_min | df$b < cell_min |
                 df$c < cell_min | df$d < cell_min)
  df$P_FISHER <- NA_real_
  if (length(idx) > 0) {
    df$P_FISHER[idx] <- vapply(idx, function(i) {
      m <- matrix(c(df$a[i], df$b[i], df$c[i], df$d[i]), nrow = 2)
      fisher.test(m)$p.value
    }, numeric(1))
  }
  df
}

# ------------------------------------------------------------
# 4. Head-to-head comparison: drug A vs drug B
#    (each drug serves as the other's background)
# ------------------------------------------------------------
head_to_head <- function(pairs, drugA, drugB, drug_col = "drug",
                         n_min = 3, n_A = NULL, n_B = NULL) {
  ta <- pairs %>% filter(.data[[drug_col]] == drugA) %>%
    distinct(case_id, event)
  tb <- pairs %>% filter(.data[[drug_col]] == drugB) %>%
    distinct(case_id, event)
  disprop_all(ta, tb, n_min = n_min,
              n_target_total = n_A, n_bg_total = n_B) %>%
    mutate(comparison = paste(drugA, "vs", drugB, sep = " "), .before = 1)
}

# ------------------------------------------------------------
# 5. EBGM/MGPS (optional; requires openEBGM)
#    pairs: report-level long data.frame(id, drug, event)
# ------------------------------------------------------------
run_ebgm <- function(pairs) {
  if (!requireNamespace("openEBGM", quietly = TRUE)) {
    message("openEBGM not installed; skipping EBGM. Install with install.packages('openEBGM')")
    return(NULL)
  }
  colnames(pairs) <- c("id", "var1", "var2")
  dat <- openEBGM::processRaw(pairs, stratify = FALSE, zeroes = FALSE)
  hyp <- openEBGM::autoHyper(dat, squashed = TRUE, zeroes = FALSE)
  theta <- openEBGM::Qn(hyp)
  ebgm <- openEBGM::ebgm(theta, dat)
  eb05 <- openEBGM::quantBisect(5, theta = theta, N = dat$N, E = dat$E)
  tibble(drug = dat$var1, event = dat$var2, N = dat$N, E = dat$E,
         EBGM = ebgm, EB05 = eb05) %>%
    filter(N >= 3)
}

# ------------------------------------------------------------
# 6. Plotting
# ------------------------------------------------------------

# Forest plot (log scale)
plot_forest <- function(df, label_col, title = "",
                        ror_col = "ROR", lo_col = "ROR025", hi_col = "ROR975") {
  df <- df %>%
    mutate(lab = .data[[label_col]],
           x = .data[[ror_col]],
           lo = .data[[lo_col]],
           hi = .data[[hi_col]]) %>%
    filter(is.finite(x), is.finite(lo), is.finite(hi), lo > 0)
  ggplot(df, aes(x = x, y = reorder(lab, x))) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, colour = "grey40") +
    geom_point(size = 2.2, colour = "#B2182B") +
    scale_x_log10() +
    labs(x = "ROR (95% CI, log scale)", y = NULL, title = title) +
    theme_bw(base_size = 12)
}

# Publication-style volcano plot: log2(ROR) vs -log10(FDR P)
#   point size ~ report count a; red = meets all signal criteria;
#   labels are drawn only for the top signals to avoid overcrowding
plot_volcano <- function(df, label_col = "event", top_n = 12, title = "",
                         y_cap = 300) {
  empty_plot <- function(msg) {
    message("!! plot_volcano: ", title, " - ", msg)
    ggplot() +
      annotate("text", x = 0, y = 0, size = 5, label = paste0(title, "\n", msg)) +
      theme_bw(base_size = 12)
  }
  if (is.null(df) || nrow(df) == 0)
    return(empty_plot("NO DATA - check upstream"))

  df <- df %>% filter(is.finite(ROR), ROR > 0)
  # Prefer the log10-scale FDR (exact, no underflow); fall back to
  # plain P_FDR with a lower-bound floor
  if ("LOG10P_FDR" %in% names(df) && any(!is.na(df$LOG10P_FDR))) {
    df <- df %>% filter(!is.na(LOG10P_FDR)) %>% mutate(nlp = -LOG10P_FDR)
  } else {
    df <- df %>%
      filter(!is.na(P_FDR)) %>%
      mutate(P_FDR_plot = ifelse(P_FDR <= 0,
                                 min(P_FDR[P_FDR > 0], na.rm = TRUE) / 10,
                                 P_FDR),
             nlp = -log10(P_FDR_plot))
  }
  df <- df %>% mutate(l2r = log2(ROR), sig = signal %in% TRUE)
  if (nrow(df) == 0 || all(is.na(df$nlp)))
    return(empty_plot("ALL P_FDR MISSING"))

  # Clip the y axis: extremely significant events (nlp can reach the
  # thousands) are capped at y_cap and marked with triangles so a
  # single point does not flatten the whole plot
  df <- df %>%
    mutate(clipped = nlp > y_cap,
           nlp_plot = pmin(nlp, y_cap * 1.02))
  n_clip <- sum(df$clipped)

  # Label priority: top_n signals by -log10(P_FDR); if fewer than
  # top_n signals exist, fill with the overall top events so the
  # plot always carries some labels
  lab <- df %>% filter(sig) %>% arrange(desc(nlp)) %>% slice_head(n = top_n)
  if (nrow(lab) < top_n) {
    lab <- df %>% arrange(desc(nlp)) %>% slice_head(n = top_n)
  }

  y05 <- -log10(0.05)

  p <- ggplot(df, aes(x = l2r, y = nlp_plot)) +
    geom_point(data = df %>% filter(!sig),
               colour = "grey82", alpha = 0.5, size = 1.1) +
    geom_point(data = df %>% filter(sig, !clipped),
               aes(size = log10(a)), colour = "#B2182B", alpha = 0.85) +
    geom_point(data = df %>% filter(sig, clipped),
               aes(size = log10(a)), colour = "#B2182B", alpha = 0.9,
               shape = 17) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40",
               linewidth = 0.4) +
    geom_hline(yintercept = y05, linetype = 3, colour = "grey60",
               linewidth = 0.35) +
    annotate("text", x = -Inf, y = y05, label = "FDR = 0.05",
             hjust = -0.05, vjust = -0.4, size = 3, colour = "grey50") +
    ggrepel::geom_text_repel(
      data = lab,
      aes(label = stringr::str_wrap(stringr::str_to_title(.data[[label_col]]), 22)),
      size = 3.1, fontface = "bold", colour = "grey15",
      box.padding = 0.5, point.padding = 0.3,
      segment.colour = "grey55", segment.size = 0.3,
      force = 2, force_pull = 0.5, max.overlaps = Inf,
      min.segment.length = 0, seed = 42
    ) +
    scale_size_continuous(range = c(1.6, 4.5), guide = "none") +
    coord_cartesian(ylim = c(NA, y_cap * 1.08)) +
    labs(x = expression(log[2](ROR)),
         y = expression(-log[10]("FDR-adjusted"~italic(P))),
         title = title,
         caption = if (n_clip > 0)
           paste0("Triangles: ", n_clip, " event(s) with -log10(FDR P) > ",
                  y_cap, ", truncated for display") else NULL) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.caption = element_text(colour = "grey50", size = 8, hjust = 0),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
      axis.title = element_text(size = 12)
    )
  p
}

# Cross-database forest plot (one event, several sources)
plot_cross_db <- function(df, event_order = NULL) {
  # df must contain: event, source, ROR, ROR025, ROR975
  d <- df %>% filter(is.finite(ROR), ROR025 > 0)
  if (!is.null(event_order)) d$event <- factor(d$event, levels = event_order)
  ggplot(d, aes(x = ROR, y = source, colour = source)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = ROR025, xmax = ROR975), height = 0.2) +
    geom_point(size = 2.2) +
    scale_x_log10() +
    facet_wrap(~event, scales = "free_y") +
    labs(x = "ROR (95% CI, log scale)", y = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
}

# ------------------------------------------------------------
# 7. Stale-workdir guard: called at the start of 02/03/04/06/10.
#    If the exclusion list exists but pairs_pt still contains those
#    PTs, the workdir in the current directory was built by an old
#    version of 01 (or the wrong project folder is open).
# ------------------------------------------------------------
assert_fresh_workdir <- function() {
  meta_file <- "workdir/_build_info.rds"
  # Hard requirement: workdir must carry the build fingerprint written
  # by 01. Workdirs from older pipelines lack this file -> stop here.
  if (!file.exists(meta_file)) {
    stop("workdir lacks the build fingerprint _build_info.rds - this workdir was generated by an old version of 01",
         "\n  Current working directory: ", getwd(),
         "\n  Please check: 1) rerun 01_faers_prepare.R in this directory with the current code;",
         "\n         2) whether multiple copies of the project folder exist and the wrong one is open",
         call. = FALSE)
  }
  info <- readRDS(meta_file)
  if (!isTRUE(info$exclude_on) || (info$excluded_pts_n %||% 0) == 0) {
    stop("workdir was built by 01 with non-AE-PT exclusion disabled; please rerun 01_faers_prepare.R",
         call. = FALSE)
  }
  # Double check: when the exclusion list is available, verify no
  # excluded PT leaks through
  excl_file <- "data/excluded_pts.csv"
  pf <- "workdir/pairs_pt.rds"
  if (file.exists(excl_file) && file.exists(pf)) {
    excl <- toupper(trimws(readr::read_csv(excl_file,
                                           show_col_types = FALSE)$pt))
    leak <- sum(readRDS(pf)$event %in% excl)
    if (leak > 0) {
      stop("workdir was built by an old version of 01 (pairs_pt still contains ", leak,
           " rows of excluded non-AE PTs)\n  Current working directory: ", getwd(),
           "\n  Please rerun 01_faers_prepare.R in this directory",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

# ------------------------------------------------------------
# 8. Output helper
# ------------------------------------------------------------
write_tab <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(df, path)
  message("Written: ", path, " (", nrow(df), " rows)")
  invisible(path)
}

# ------------------------------------------------------------
# 9. Load and merge per-quarter RDS files
#    For objects saved one quarter at a time (each .rds is a list of
#    tables: demo/drug/reac/...).
#    Memory strategy: trim each quarter to the required columns
#    immediately after reading (rpsr is never loaded), then
#    incrementally combine with data.table::rbindlist so full-width
#    tables never pile up in memory.
#    After merging, de-duplicate across quarters following FDA
#    recommendations: for the same caseid keep the version with the
#    most recent fda_dt and highest caseversion.
# ------------------------------------------------------------

# Columns required by the analyses (00a stores all columns; only
# these are kept here)
NEEDED_COLS <- list(
  demo = c("primaryid", "caseid", "caseversion", "fda_dt", "event_dt",
           "age", "age_cod", "sex", "occp_cod", "occr_country",
           "reporter_country"),
  drug = c("primaryid", "drug_seq", "role_cod", "drugname", "prod_ai"),
  reac = c("primaryid", "pt"),
  outc = c("primaryid", "outc_cod"),
  ther = c("primaryid", "drug_seq", "dsg_drug_seq", "start_dt"),
  indi = c("primaryid", "drug_seq", "indi_drug_seq", "indi_pt")
)

# In the current FAERS format the sequence columns of THER/INDI have
# different names from DRUG; harmonise them
SEQ_RENAME <- c(dsg_drug_seq = "drug_seq", indi_drug_seq = "drug_seq")

load_faers_quarters <- function(dir,
                                pattern = "\\.rds$",
                                tables = names(NEEDED_COLS),
                                keep_cols = NEEDED_COLS,
                                years = NULL) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE,
                      ignore.case = TRUE)
  if (length(files) == 0) stop("No .rds files found in: ", dir)
  files <- sort(files)

  # Pre-filter by the year parsed from the file name
  # (e.g. "2022q1.rds" -> 2022). Files whose names yield no year are
  # kept, to avoid dropping anything by mistake.
  if (!is.null(years)) {
    yr <- suppressWarnings(
      as.integer(stringr::str_extract(basename(files), "(19|20)\\d{2}"))
    )
    keep <- is.na(yr) | yr %in% years
    skipped <- files[!keep & !is.na(yr)]
    if (length(skipped) > 0) {
      message("Year filter: skipping ", length(skipped), " files: ",
              paste(head(basename(skipped), 3), collapse = ", "),
              if (length(skipped) > 3) " ..." else "")
    }
    files <- files[keep]
    if (length(files) == 0) stop("No files left after year filtering; check the years argument or file names")
  }

  message("Found ", length(files), " quarterly files:")
  message(paste("  -", basename(files), collapse = "\n"))

  # Read quarter by quarter -> trim columns immediately -> collect
  # incrementally (full-width tables never all reside in memory)
  buckets <- setNames(vector("list", length(tables)), tables)
  for (tb in tables) buckets[[tb]] <- vector("list", length(files))

  for (i in seq_along(files)) {
    o <- readRDS(files[i])
    if (!is.list(o) || !"demo" %in% names(o)) {
      stop(files[i], " is not a 'list of demo/drug/reac/... tables'; ",
           "inspect it with names(readRDS(...))")
    }
    for (tb in tables) {
      dt <- o[[tb]]
      if (is.null(dt)) next
      want <- intersect(keep_cols[[tb]], names(dt))
      sub <- dplyr::select(as_tibble(dt), all_of(want))
      # Harmonise sequence column names:
      # dsg_drug_seq / indi_drug_seq -> drug_seq
      for (old_nm in intersect(names(SEQ_RENAME), names(sub))) {
        if (!SEQ_RENAME[[old_nm]] %in% names(sub)) {
          names(sub)[names(sub) == old_nm] <- SEQ_RENAME[[old_nm]]
        }
      }
      buckets[[tb]][[i]] <- sub
    }
    rm(o)
    if (i %% 4 == 0) gc(verbose = FALSE)
  }

  merged <- lapply(buckets, function(lst) {
    lst <- lst[!vapply(lst, is.null, logical(1))]
    if (length(lst) == 0) return(NULL)
    as_tibble(data.table::rbindlist(lst, use.names = TRUE, fill = TRUE))
  })
  merged <- merged[!vapply(merged, is.null, logical(1))]

  # Cross-quarter de-duplication (driven by demo; other tables are
  # filtered to the retained primaryids).
  # Rule (FDA recommendation): for the same caseid keep the record
  # with the most recent fda_dt, breaking ties by highest
  # caseversion. primaryid may be a synthetic character ID (legacy
  # AERS uses the form "2004Q1_1234567"), so no numeric comparison.
  if (all(c("caseid", "primaryid") %in% names(merged$demo))) {
    d <- merged$demo
    if (!"caseversion" %in% names(d)) d$caseversion <- NA
    if (!"fda_dt" %in% names(d))      d$fda_dt <- NA
    keep_ids <- d %>%
      mutate(.fdt = parse_any_date(fda_dt),
             .cv  = suppressWarnings(as.numeric(caseversion))) %>%
      arrange(suppressWarnings(as.numeric(caseid)),
              desc(.fdt), desc(.cv), desc(primaryid)) %>%
      distinct(caseid, .keep_all = TRUE) %>%
      pull(primaryid)
    merged <- lapply(merged, function(tb) {
      if ("primaryid" %in% names(tb)) {
        tb %>% filter(primaryid %in% keep_ids)
      } else tb
    })
    message("Cross-quarter de-duplication done; latest-version cases kept: ", length(keep_ids))
  } else {
    message("Note: demo has no caseid/primaryid columns; skipping cross-quarter de-duplication")
  }

  merged
}

# Robust FAERS date parsing (handles yyyymmdd as character or
# numeric, and Date)
parse_any_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  x <- as.character(x)
  x[nchar(x) != 8] <- NA
  ymd(x, quiet = TRUE)
}

# ------------------------------------------------------------
# 10. Build the PT -> primary SOC map (two channels)
#     Channel A (preferred): pass an unzipped MedDRA folder with the
#         faers package installed
#         -> faers::meddra(path, primary_soc = TRUE) + meddra_hierarchy()
#     Channel B (fallback): pass a zip or folder; parse mdhier.asc +
#         soc.asc manually (mdhier field 11, pt_soc_code, is the
#         primary SOC code of each PT)
#     The result is cached as a CSV on first run and read from the
#     cache afterwards.
# ------------------------------------------------------------

# Robust unzip: on macOS/Linux prefer the system unzip (handles
# zip64/large files); fall back to R's internal unzip on failure
unzip_robust <- function(zip, exdir) {
  if (.Platform$OS.type == "unix") {
    rc <- tryCatch(
      system2("unzip",
              args = c("-o", "-q", shQuote(zip), "-d", shQuote(exdir))),
      warning = function(w) 1L, error = function(e) 1L
    )
    if (identical(as.integer(rc), 0L)) return(invisible(TRUE))
    message("System unzip returned a non-zero status; falling back to R's internal unzip ...")
  }
  utils::unzip(zip, exdir = exdir)
  invisible(TRUE)
}

build_pt_soc_map <- function(meddra_path,
                             cache_csv = "data/meddra_pt_soc_map.csv") {
  if (file.exists(cache_csv)) {
    message("Reading cached MedDRA map: ", cache_csv)
    return(readr::read_csv(cache_csv, show_col_types = FALSE) %>%
             mutate(pt = toupper(pt), soc = toupper(soc)))
  }
  if (!file.exists(meddra_path)) stop("MedDRA path does not exist: ", meddra_path)

  # ---- Channel A: faers package + unzipped folder ----
  if (dir.exists(meddra_path) && requireNamespace("faers", quietly = TRUE)) {
    message("Reading MedDRA via the faers package (primary_soc = TRUE): ", meddra_path)
    m <- faers::meddra(meddra_path, primary_soc = TRUE)
    h <- as_tibble(faers::meddra_hierarchy(m))
    pt_col  <- intersect(c("pt_name", "pt"), names(h))[1]
    soc_col <- intersect(c("soc_name", "soc"), names(h))[1]
    if (is.na(pt_col) || is.na(soc_col)) {
      stop("Could not identify PT/SOC columns in meddra_hierarchy output")
    }
    map <- h %>%
      transmute(pt = .data[[pt_col]], soc = .data[[soc_col]]) %>%
      distinct()

  } else {
    # ---- Channel B: manual parsing of the .asc files ----
    if (dir.exists(meddra_path)) {
      message("MedDRA path is a folder; looking for .asc files inside it")
      td <- meddra_path
    } else {
      message("Building PT->SOC map from MedDRA zip (one-off): ", meddra_path)
      td <- tempfile("meddra_"); dir.create(td)
      on.exit(unlink(td, recursive = TRUE), add = TRUE)
      unzip_robust(meddra_path, exdir = td)
    }

    f_hier <- list.files(td, pattern = "^mdhier\\.asc$", recursive = TRUE,
                         full.names = TRUE, ignore.case = TRUE)
    f_soc  <- list.files(td, pattern = "^soc\\.asc$", recursive = TRUE,
                         full.names = TRUE, ignore.case = TRUE)
    if (length(f_hier) == 0 || length(f_soc) == 0) {
      stop("mdhier.asc / soc.asc not found. ",
           "Tip: unzip the MedDRA zip manually and point MEDDRA_ZIP to the unzipped folder",
           call. = FALSE)
    }
    f_hier <- f_hier[1]; f_soc <- f_soc[1]

    read_asc <- function(f) {
      dt <- data.table::fread(f, sep = "$", header = FALSE,
                              colClasses = "character", encoding = "UTF-8",
                              showProgress = FALSE)
      # A trailing $ produces an all-empty last column; drop it
      dt[, which(vapply(dt, function(x) all(is.na(x) | x == ""), logical(1))) := NULL]
      dt
    }

    hier <- read_asc(f_hier)
    soc  <- read_asc(f_soc)
    message("mdhier.asc: ", ncol(hier), " cols, ", nrow(hier), " rows; ",
            "soc.asc: ", ncol(soc), " cols, ", nrow(soc), " rows")

    # Standard 11-field mdhier layout (MedDRA ASCII specification):
    # 1 pt_code, 2 hlt_code, 3 hlgt_code, 4 soc_code, 5 pt_name,
    # 6 hlt_name, 7 hlgt_name, 8 soc_name, 9 soc_abbrev,
    # 10 pt_soc_code (primary SOC code), 11 primary_soc_fg (Y/N)
    names(hier) <- c("pt_code", "hlt_code", "hlgt_code", "soc_code",
                     "pt_name", "hlt_name", "hlgt_name", "soc_name",
                     "soc_abbrev", "pt_soc_code", "primary_soc_fg")[seq_len(ncol(hier))]
    names(soc) <- c("soc_code", "soc_name", "soc_abbrev",
                    "null_field")[seq_len(ncol(soc))]

    # Field sanity check: pt_soc_code should be all numeric codes
    # (if it has been misaligned with primary_soc_fg it will be all
    # Y/N); fix automatically on mismatch
    if (all(c("pt_soc_code", "primary_soc_fg") %in% names(hier))) {
      samp <- trimws(na.omit(hier$pt_soc_code))
      if (length(samp) > 0 &&
          mean(grepl("^[0-9]+$", samp)) < 0.5 &&
          mean(grepl("^[0-9]+$", trimws(na.omit(hier$soc_code)))) > 0.9) {
        message("!! pt_soc_code field misalignment detected; using soc_code instead")
        hier$pt_soc_code <- hier$soc_code
      }
    }

    # Fallback: if a given mdhier version lacks field 11
    # (pt_soc_code), use field 4 soc_code instead (may contain
    # multi-axial mappings; take the first row per PT and warn)
    use_primary <- "pt_soc_code" %in% names(hier) &&
      any(!is.na(hier$pt_soc_code) & trimws(hier$pt_soc_code) != "")
    if (!use_primary) {
      message("!! Warning: mdhier has no valid pt_soc_code field; ",
              "falling back to the hierarchy SOC field (may contain multi-axial mappings)")
    }

    map <- hier %>%
      transmute(
        pt = trimws(pt_name),
        soc_code = trimws(if (use_primary) pt_soc_code else soc_code)
      ) %>%
      filter(!is.na(soc_code), soc_code != "") %>%
      distinct() %>%
      # Trim the join key on both sides (guard against CR/space/BOM)
      left_join(
        soc %>% transmute(soc_code = trimws(soc_code), soc = soc_name),
        by = "soc_code"
      ) %>%
      filter(!is.na(soc)) %>%
      group_by(pt) %>%
      slice_head(n = 1) %>%     # one SOC per PT
      ungroup() %>%
      dplyr::select(pt, soc)

    # Guard against an empty cache: 0 rows means parsing failed
    if (nrow(map) == 0) {
      stop("Map has 0 rows; the join failed", call. = FALSE)
    }
    if (nrow(map) < 20000) {
      warning("!! Unusually few PTs in the map: ", nrow(map),
              " (MedDRA 27.x should have ~26,000); check the .asc parsing", call. = FALSE)
    }
  }

  message("Map built: ", nrow(map), " PTs")
  dir.create(dirname(cache_csv), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(map, cache_csv)
  message("Cached to: ", cache_csv, " (read from cache on later runs)")

  map %>% mutate(pt = toupper(pt), soc = toupper(soc))
}
