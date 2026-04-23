########################## Build failures.rds #################################
## Reads the aggregation failure log and builds a conditions-style data.frame
## with one row per failed SLURM task, matching conditions.rds in columns and
## types so it can be swapped in as a drop-in rerun grid.
################################################################################

fail_csv <- "simresults/2026_04_22_arrayresults/aggregation_failures.csv"
out_file <- "arraysimulation/failures.rds"

fails <- read.csv2(fail_csv, stringsAsFactors = FALSE)

pat <- "model_([0-9]+)_sample_size_([0-9]+)_batch_([0-9]+)_datatype_([0-9]+)\\.rds"
m <- regmatches(fails$file, regexec(pat, fails$file))
parsed <- do.call(rbind, lapply(m, function(x){
  if(length(x) == 5) as.integer(x[2:5]) else rep(NA_integer_, 4)
}))
colnames(parsed) <- c("model", "n", "batch", "datatype")

# Match conditions.rds column order and types (rep_batch int, model int,
# n num, datatype num)
failures <- data.frame(
  rep_batch = as.integer(parsed[, "batch"]),
  model     = as.integer(parsed[, "model"]),
  n         = as.numeric(parsed[, "n"]),
  datatype  = as.numeric(parsed[, "datatype"])
)

# Sort to mirror conditions.rds ordering (rep_batch varies fastest, then model,
# then n, then datatype — i.e. expand.grid's natural order).
failures <- failures[order(failures$datatype, failures$n,
                           failures$model, failures$rep_batch), ]
rownames(failures) <- NULL

cat("Parsed failures:", nrow(failures), "\n")
cat("\ndistribution:\n")
print(table(failures$datatype, failures$n, dnn = c("datatype", "n")))

saveRDS(failures, out_file)
cat("\nWrote:", out_file, "\n")
