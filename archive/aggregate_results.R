########################## Aggregate array-job RDS outputs ####################
## Reads every per-task .rds file written by sim.R into one data.frame.
## Files that fail to read (corrupted, wrong object, etc.) are skipped and
## logged to a failure CSV rather than killing the aggregation.
################################################################################

library(dplyr)

results_dir  <- "simresults/2026_04_25_arrayresults/results"
out_dir      <- "simresults/2026_04_25_arrayresults/"
out_file     <- file.path(out_dir, "aggregated.rds")
fail_file    <- file.path(out_dir, "aggregation_failures.csv")

stopifnot(dir.exists(results_dir))

files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
cat("Found", length(files), "rds files in", results_dir, "\n")

# Read one file defensively: returns either a data.frame or a single-row
# data.frame describing the failure (status = "error"). This keeps the two
# outcomes in the same list so we can separate them at the end.
read_one <- function(path){
  tryCatch({
    obj <- readRDS(path)
    if(!is.data.frame(obj)){
      return(data.frame(file = basename(path),
                        status = "not_a_data_frame",
                        message = paste("class:", paste(class(obj), collapse = "/")),
                        stringsAsFactors = FALSE))
    }
    if(nrow(obj) == 0){
      return(data.frame(file = basename(path),
                        status = "empty",
                        message = "0 rows",
                        stringsAsFactors = FALSE))
    }
    obj$source_file <- basename(path)
    obj
  },
  error = function(e){
    data.frame(file = basename(path),
               status = "read_error",
               message = conditionMessage(e),
               stringsAsFactors = FALSE)
  },
  warning = function(w){
    # Treat warnings as non-fatal: try to read once more without the handler,
    # but surface the warning alongside success.
    obj <- try(readRDS(path), silent = TRUE)
    if(inherits(obj, "try-error") || !is.data.frame(obj)){
      return(data.frame(file = basename(path),
                        status = "read_warning_then_error",
                        message = conditionMessage(w),
                        stringsAsFactors = FALSE))
    }
    obj$source_file <- basename(path)
    attr(obj, "read_warning") <- conditionMessage(w)
    obj
  })
}

read_list <- vector("list", length(files))
for(i in seq_along(files)){
  read_list[[i]] <- read_one(files[i])
  if(i %% 100 == 0) cat("  ...", i, "/", length(files), "\n")
}

is_failure <- vapply(read_list, function(x){
  is.data.frame(x) && all(c("file", "status", "message") %in% names(x)) &&
    nrow(x) == 1 && !"correlation" %in% names(x)
}, logical(1))

failures <- read_list[is_failure]
successes <- read_list[!is_failure]

cat("\nSuccessful reads: ", length(successes), "\n")
cat("Failed reads:     ", length(failures),  "\n")

if(length(failures) > 0){
  fail_df <- bind_rows(failures)
  write.csv2(fail_df, fail_file, row.names = FALSE)
  cat("Failure log written to:", fail_file, "\n")
  print(table(fail_df$status))
}

if(length(successes) == 0){
  stop("No files could be read. Aborting before writing an empty aggregate.")
}

aggregated <- bind_rows(successes)
cat("\nAggregated data.frame:", nrow(aggregated), "rows x",
    ncol(aggregated), "cols\n")

# time is a difftime in the per-task outputs; convert to numeric seconds so
# the column can be bound/saved cleanly downstream.
if("time" %in% names(aggregated) && inherits(aggregated$time, "difftime")){
  aggregated$time <- as.numeric(aggregated$time, units = "secs")
}

saveRDS(aggregated, out_file)
cat("Aggregated RDS written to:", out_file, "\n")


per_file <- aggregated %>%
  group_by(source_file) %>%
  summarize(n_realizations = n_distinct(simruns), .groups = "drop")

cat("\nFiles with fewer than 100 realizations:\n")
print(per_file[per_file$n_realizations < 100, ])


