########################## Overview table ####################################
## Builds a per-condition summary of basic statistics from the aggregated
## simulation results. One row per (method, n, correlation, model, datatype,
## alpha) cell.
###############################################################################

library(dplyr)
library(tidyr)

agg_path <- "simresults/2026_04_25_arrayresults/aggregated.rds"
out_rds  <- "simresults/2026_04_25_arrayresults/overview_table.rds"
out_csv  <- "simresults/2026_04_25_arrayresults/overview_table.csv"

dfall <- readRDS(agg_path)

# model id is not a column on its own — pull it out of the source filename
# (e.g. "model_3_sample_size_200_batch_5_datatype_1.rds").
dfall$model <- as.integer(sub(".*model_(\\d+)_sample_size_.*", "\\1",
                              dfall$source_file))

# CI width and bound-side coverage (matches plotting.R conventions).
dfall$ci_width    <- dfall$upperbound - dfall$lowerbound
dfall$upperwithin <- dfall$correlation < dfall$upperbound
dfall$lowerwithin <- dfall$correlation > dfall$lowerbound

ag1 <- dfall %>%
  group_by(method, n, correlation, model, datatype, alpha) %>%
  summarize(
    n_reps          = n(),
    # point estimate
    mean_cor_est    = mean(cor_est,    na.rm = TRUE),
    sd_cor_est      = sd(cor_est,      na.rm = TRUE),
    mean_se         = mean(se_cor_est, na.rm = TRUE),
    # CI bounds and width
    mean_lower      = mean(lowerbound, na.rm = TRUE),
    mean_upper      = mean(upperbound, na.rm = TRUE),
    mean_width      = mean(ci_width,   na.rm = TRUE),
    # coverage in % (matches plotting.R scale)
    cov_corr_pct    = mean(coveragecorr, na.rm = TRUE) * 100,
    cov_one_pct     = mean(coverageone,  na.rm = TRUE) * 100,
    upper_within_pct = mean(upperwithin, na.rm = TRUE) * 100,
    lower_within_pct = mean(lowerwithin, na.rm = TRUE) * 100,
    # runtime
    mean_time       = mean(time, na.rm = TRUE),
    # diagnostics
    mean_missing    = mean(missing, na.rm = TRUE),
    na_cor_est_pct  = mean(is.na(cor_est))    * 100,
    na_lower_pct    = mean(is.na(lowerbound)) * 100,
    na_upper_pct    = mean(is.na(upperbound)) * 100,
    warning_pct     = mean(!is.na(warning))   * 100,
    error_pct       = mean(!is.na(error))     * 100,
    .groups = "drop"
  ) %>%
  arrange(method, datatype, alpha, model, n)

cat("ag1 dimensions:", nrow(ag1), "rows x", ncol(ag1), "cols\n")
cat("cells per method:\n"); print(table(ag1$method))

saveRDS(ag1, out_rds)
write.csv2(ag1, out_csv, row.names = FALSE)
cat("Wrote:", out_rds, "\n")
cat("Wrote:", out_csv, "\n")
