################################################################################
## Investigate which simulation cells fail to produce a correlation estimate.
## We look at `cor_est` (NA = correlation could not be calculated) and look for
## systematic patterns by n, correlation, datatype, model and method.
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

results_dir <- "simresults/2026_04_23_arrayresults/results"
agg_file    <- "simresults/aggregated.rds"

if(file.exists(agg_file)){
  cat("Loading existing aggregated.rds\n")
  aggregated <- readRDS(agg_file)
} else {
  cat("Aggregated file not found - rebuilding from per-task RDS files...\n")
  files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
  cat("Found", length(files), "rds files\n")

  read_one <- function(path){
    out <- tryCatch(readRDS(path), error = function(e) NULL,
                    warning = function(w) tryCatch(readRDS(path),
                                                   error = function(e) NULL))
    if(is.null(out) || !is.data.frame(out) || nrow(out) == 0) return(NULL)
    # encode model/datatype from filename so we don't lose them
    bn <- basename(path)
    m  <- regmatches(bn, regexpr("model_\\d+", bn))
    d  <- regmatches(bn, regexpr("datatype_\\d+", bn))
    out$model_id    <- as.integer(sub("model_", "", m))
    out$datatype_id <- as.integer(sub("datatype_", "", d))
    out$source_file <- bn
    out
  }

  read_list <- vector("list", length(files))
  for(i in seq_along(files)){
    read_list[[i]] <- read_one(files[i])
    if(i %% 200 == 0) cat("  ...", i, "/", length(files), "\n")
  }
  read_list <- read_list[!vapply(read_list, is.null, logical(1))]
  aggregated <- bind_rows(read_list)
  cat("Aggregated rows:", nrow(aggregated), "\n")
  saveRDS(aggregated, agg_file)
}

cat("\n=== Dimensions ===\n")
cat(nrow(aggregated), "rows x", ncol(aggregated), "cols\n")
cat("\n=== Columns ===\n")
print(names(aggregated))

# Make sure model/datatype columns exist (older builds might not have them)
if(!"model_id" %in% names(aggregated)){
  aggregated$model_id <- as.integer(sub(".*model_(\\d+)_.*", "\\1",
                                        aggregated$source_file))
}
if(!"datatype_id" %in% names(aggregated)){
  aggregated$datatype_id <- as.integer(sub(".*datatype_(\\d+)\\.rds$", "\\1",
                                           aggregated$source_file))
}

# alpha defines a per-rep replicate; collapse to one row per (rep, method)
# for the NA-pattern analysis -- cor_est does not depend on alpha.
df <- aggregated %>% filter(alpha == 0.05)

cat("\n=== Overall NA rate of cor_est by method ===\n")
overall <- df %>%
  group_by(method) %>%
  summarize(n_total = n(),
            n_na    = sum(is.na(cor_est)),
            pct_na  = round(100 * mean(is.na(cor_est)), 3),
            .groups = "drop")
print(overall)

cat("\n=== NA rate by method x sample size ===\n")
by_n <- df %>%
  group_by(method, n) %>%
  summarize(n_total = n(),
            n_na    = sum(is.na(cor_est)),
            pct_na  = round(100 * mean(is.na(cor_est)), 2),
            .groups = "drop") %>%
  arrange(method, n)
print(by_n, n = Inf)

cat("\n=== NA rate by method x datatype ===\n")
by_dt <- df %>%
  group_by(method, datatype) %>%
  summarize(n_total = n(),
            n_na    = sum(is.na(cor_est)),
            pct_na  = round(100 * mean(is.na(cor_est)), 2),
            .groups = "drop") %>%
  arrange(method, datatype)
print(by_dt, n = Inf)

cat("\n=== NA rate by method x correlation (true value) ===\n")
by_corr <- df %>%
  group_by(method, correlation) %>%
  summarize(n_total = n(),
            n_na    = sum(is.na(cor_est)),
            pct_na  = round(100 * mean(is.na(cor_est)), 2),
            .groups = "drop") %>%
  arrange(method, correlation)
print(by_corr, n = Inf)

cat("\n=== NA rate by method x model_id ===\n")
by_model <- df %>%
  group_by(method, model_id) %>%
  summarize(n_total = n(),
            n_na    = sum(is.na(cor_est)),
            pct_na  = round(100 * mean(is.na(cor_est)), 2),
            .groups = "drop") %>%
  arrange(method, model_id)
print(by_model, n = Inf)

cat("\n=== Cells with HIGH NA rate (>=10%): method x n x correlation x datatype ===\n")
hot <- df %>%
  group_by(method, n, correlation, datatype) %>%
  summarize(n_total = n(),
            n_na    = sum(is.na(cor_est)),
            pct_na  = round(100 * mean(is.na(cor_est)), 2),
            .groups = "drop") %>%
  filter(pct_na >= 10) %>%
  arrange(desc(pct_na))
print(hot, n = Inf)

# ============================================================================
# Look at the "conphi" (CFA) method specifically -- those rows store error/
# warning text from the lavaan fits.
# ============================================================================
cphi <- df %>% filter(method == "conphi")

cat("\n=== conphi: counts of unique error_unconstrained messages ===\n")
err_tbl <- cphi %>%
  filter(!is.na(error_unconstrained) | is.na(cor_est)) %>%
  count(error_unconstrained, sort = TRUE)
print(err_tbl, n = 30)

cat("\n=== conphi: counts of unique warning_unconstrained messages (when NA) ===\n")
warn_tbl <- cphi %>%
  filter(is.na(cor_est)) %>%
  count(warning_unconstrained, sort = TRUE)
print(warn_tbl, n = 30)

cat("\n=== conphi NA breakdown: error vs warning vs neither ===\n")
cphi_na <- cphi %>% filter(is.na(cor_est))
cat("Total conphi NA rows:", nrow(cphi_na), "\n")
cat("  with error message:  ",
    sum(!is.na(cphi_na$error_unconstrained)), "\n")
cat("  with warning only:   ",
    sum(is.na(cphi_na$error_unconstrained) &
        !is.na(cphi_na$warning_unconstrained)), "\n")
cat("  neither (silent NA): ",
    sum(is.na(cphi_na$error_unconstrained) &
        is.na(cphi_na$warning_unconstrained)), "\n")

# ============================================================================
# Compare HTMT (delta) NA pattern vs conphi NA pattern -- do they fail in the
# same reps?
# ============================================================================
wide <- df %>%
  select(simruns, batch, n, correlation, datatype, model_id, method, cor_est) %>%
  pivot_wider(names_from = method, values_from = cor_est)

cat("\n=== Co-occurrence of NA across methods (per rep) ===\n")
co <- wide %>%
  summarize(n_reps          = n(),
            na_delta        = sum(is.na(delta)),
            na_conphi       = sum(is.na(conphi)),
            na_boot         = sum(is.na(boot)),
            na_bcaboot      = sum(is.na(bcaboot)),
            na_bcboot       = sum(is.na(bcboot)),
            both_delta_conphi = sum(is.na(delta) & is.na(conphi)),
            delta_only      = sum(is.na(delta) & !is.na(conphi)),
            conphi_only     = sum(!is.na(delta) & is.na(conphi)))
print(co)

cat("\n=== HTMT (delta) NA: cells with high rate ===\n")
htmt_hot <- df %>%
  filter(method == "delta") %>%
  group_by(n, correlation, datatype) %>%
  summarize(n_total = n(),
            n_na = sum(is.na(cor_est)),
            pct_na = round(100*mean(is.na(cor_est)), 2),
            .groups = "drop") %>%
  filter(n_na > 0) %>%
  arrange(desc(pct_na))
print(htmt_hot, n = Inf)

cat("\nDone.\n")
