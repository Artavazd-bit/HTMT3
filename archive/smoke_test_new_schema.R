################################################################################
## Smoke test the new sim.R schema on a single small condition (the n=25 +
## severenonnormal cell, where ~5% of reps had NA cor_est in the previous run).
## We just want to confirm:
##   (a) the script no longer references the four old columns
##   (b) the new `error` and `warning` columns exist and are populated when the
##       monoblock-negative warning fires.
################################################################################

suppressPackageStartupMessages({
  library(lavaan)
  library(dplyr)
  library(foreach)
  library(covsim)
})

# point the helpers at the local copy
source("arraysimulation/setup.R")

# pick the cell that previously produced the most NAs: model 1 (rho=0.7),
# n=25, severenonnormal (data_id = 3).
model_id <- 1
data_id  <- 3
n        <- 25
current_model <- simModels$model[model_id]
correlation   <- simModels$correlation[model_id]
latent1 <- "xi_1"; latent2 <- "xi_2"

fit     <- sem(current_model, do.fit = FALSE)
popcov  <- lavInspect(fit, "implied")$cov
popskew <- rep(3, 6); popkurt <- rep(21, 6)

flat_cond <- function(x){
  if(length(x) == 0) return(NA_character_)
  paste(vapply(x, conditionMessage, character(1)), collapse = " | ")
}
combine_two <- function(a, b){
  parts <- c(if(!is.na(a)) paste0("con: ", a),
             if(!is.na(b)) paste0("uncon: ", b))
  if(length(parts) == 0) NA_character_ else paste(parts, collapse = " || ")
}

alpha_vec <- c(0.05, 0.1, 0.01)

# 30 reps is enough to see at least one monoblock-negative warning
n_reps <- 30
partial <- foreach(r = 1:n_reps, .combine = "rbind") %do% {
  set.seed(1000 + r)
  data <- as.data.frame(covsim::rPLSIM(n, sigma.target = popcov,
                                       skewness = popskew,
                                       excesskurtosis = popkurt,
                                       numsegments = 8)[[1]][[1]])
  colnames(data) <- colnames(popcov)

  na_simuresults <- list(
    delta   = list(HTMT = NA, se = NA, lowerbound = NA, upperbound = NA,
                   time = NA, missing = NA),
    boot    = list(se = NA, lowerbound = NA, upperbound = NA,
                   missing = NA, time = NA),
    bcaboot = list(se = NA, lowerbound = NA, upperbound = NA,
                   time = NA, missing = NA),
    bcboot  = list(se = NA, lowerbound = NA, upperbound = NA,
                   missing = NA, time = NA)
  )
  htmt_run <- safe_run(
    fun = wrapper, data = data, model = model_est,
    latent1 = latent1, latent2 = latent2, alpha = alpha_vec,
    scale = FALSE, htmt2 = FALSE, nboot = 50)  # nboot small for speed

  simuresults    <- if(is.na(htmt_run$error)) htmt_run$res else na_simuresults
  htmt_warn_flat <- flat_cond(htmt_run$warnings)
  htmt_err_flat  <- if(is.na(htmt_run$error)) NA_character_
                    else as.character(htmt_run$error)

  conphi <- c_phi(model_constrained = model_constrained,
                  model_unconstrained = model_est, data = data)
  uncon_valid <- inherits(conphi$unconfit, "lavaan")
  cor_est_cfa <- if(uncon_valid){
    pe <- parameterestimates(conphi$unconfit)
    pe$est[pe$lhs == latent1 & pe$rhs == latent2]
  } else NA_real_

  con_warn_flat   <- flat_cond(conphi$con_warning)
  uncon_warn_flat <- flat_cond(conphi$uncon_warning)
  con_err_flat    <- if(is.na(conphi$con_error)) NA_character_ else as.character(conphi$con_error)
  uncon_err_flat  <- if(is.na(conphi$uncon_error)) NA_character_ else as.character(conphi$uncon_error)
  conphi_warn <- combine_two(con_warn_flat, uncon_warn_flat)
  conphi_err  <- combine_two(con_err_flat,  uncon_err_flat)

  data.frame(
    rep     = r,
    method  = c("conphi", "delta", "boot", "bcaboot", "bcboot"),
    cor_est = c(cor_est_cfa,
                simuresults$delta$HTMT,
                simuresults$delta$HTMT,
                simuresults$delta$HTMT,
                simuresults$delta$HTMT),
    error   = c(conphi_err,  htmt_err_flat, htmt_err_flat, htmt_err_flat, htmt_err_flat),
    warning = c(conphi_warn, htmt_warn_flat, htmt_warn_flat, htmt_warn_flat, htmt_warn_flat),
    stringsAsFactors = FALSE
  )
}

cat("\n=== Schema check ===\n")
cat("Columns in partial:\n")
print(names(partial))
stopifnot("error"   %in% names(partial))
stopifnot("warning" %in% names(partial))
stopifnot(!any(c("error_constrained","error_unconstrained",
                 "warning_constrained","warning_unconstrained") %in% names(partial)))
cat("OK: only single error+warning columns.\n")

cat("\n=== HTMT rows: warning content where cor_est is NA ===\n")
htmt_na <- partial %>% filter(method == "delta" & is.na(cor_est))
cat("HTMT NA reps:", nrow(htmt_na), "\n")
print(htmt_na[, c("rep", "warning", "error")])

cat("\n=== HTMT rows: warning content where cor_est is NOT NA ===\n")
htmt_ok <- partial %>% filter(method == "delta" & !is.na(cor_est))
cat("HTMT OK reps:", nrow(htmt_ok),
    " | with warning text:", sum(!is.na(htmt_ok$warning)), "\n")
print(head(htmt_ok[!is.na(htmt_ok$warning), c("rep","warning")], 5))

cat("\n=== Distinct warning strings on HTMT side ===\n")
print(table(partial$warning[partial$method == "delta"], useNA = "ifany"))

cat("\nDone.\n")
