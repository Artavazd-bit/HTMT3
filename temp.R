library(dplyr)
library(tidyverse)

files <- list.files(path = "simresults/2026_04_01_arrayresults/results", pattern = ".rds", full.names = TRUE)
safe_read <- function(f) tryCatch(readRDS(f), error = function(e) { message("Skipping corrupt: ", basename(f)); NULL })
df <- files %>% map(safe_read) %>% compact() %>% list_rbind()

# Check whether each condition (correlation x n x datatype) has unique seeds
# Each seed appears multiple times due to method x alpha combinations, so check per simrun
seed_check <- df %>%
  distinct(correlation, n, datatype, seed, simruns) %>%
  group_by(correlation, n, datatype) %>%
  summarise(
    total_rows = n(),
    unique_seeds = n_distinct(seed),
    has_duplicates = n() != n_distinct(seed),
    .groups = "drop"
  )

cat("Conditions with duplicate seeds (after collapsing method/alpha):\n")
dupes <- seed_check %>% filter(has_duplicates)
if (nrow(dupes) == 0) {
  cat("None - all seeds are unique within each condition.\n")
} else {
  print(dupes)
}

cat("\nSeeds per condition:\n")
print(seed_check)
