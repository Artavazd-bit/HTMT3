.libPaths(c("/home/jab49wd/R-projects/R/library", .libPaths()))

library(lavaan)
library(dplyr)
library(foreach)
library(covsim)

source("code/setup.R")
conditions <- readRDS("code/conditions.rds")
n_batches <- max(conditions$rep_batch)

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
cond <- conditions[task_id,]

total_reps <- 1000
number_bootruns <- 1000
alpha_vec <- c(0.05, 0.1, 0.01)

rep_batch <- cond$rep_batch
reps_per_batch <- total_reps / n_batches
start_rep <- (cond$rep_batch - 1) * reps_per_batch + 1  
end_rep <- (cond$rep_batch) * reps_per_batch


n <- cond$n
n_id <- match(n, c(25, 50, 100, 200, 400, 800, 1600, 3200, 6400))
model_id <- cond$model
data_id <- cond$datatype
current_model <- simModels$model[model_id]
correlation <- simModels$correlation[model_id]
latent1 <- "xi_1"
latent2 <- "xi_2"


if(data_id == 1){
  datatype <- "normal"
}else if(data_id == 2){
  datatype <- "nonnormal"
  
  fit <- sem(current_model, do.fit = FALSE)
  popcov <- lavInspect(fit, "implied")$cov
  popskew <- c(0.7, 0.9, 1.2, 1.5, 0.8, 1.3)
  popkurt <- c(3.5, 3.6, 3.7, 3.5, 3.6, 3.7)
}

flat_cond <- function(x){
  if(length(x) == 0) return(NA_character_)
  paste(vapply(x, conditionMessage, character(1)), collapse = " | ")
}

cat(sprintf("Task %d | model=%d n=%d batch=%d datatype=%s | reps %d-%d | corr=%.2f\n",
            task_id, model_id, n, rep_batch, datatype, start_rep, end_rep, correlation))

partial <- foreach(r = start_rep : end_rep, .combine = "rbind")%do%
      {
        seed <- r + (model_id) * 10000 + (data_id) * 100000 + (n_id) * 1000000
        seed <- as.integer(seed)
        if(data_id == 1){
        data <- lavaan::simulateData(model = current_model,
                               sample.nobs = n,
                               seed = seed,
                               skewness = NULL, 
                               kurtosis = NULL, 
                               empirical = FALSE,
                               return.type = "data.frame"
        )
        }else if(data_id == 2){
          set.seed(seed)
          data <- rIG(N = n, sigma.target = popcov, skewness = popskew, 
                      excesskurtosis = popkurt, typeA = "symm")
          data <- as.data.frame(data)
          colnames(data) <- colnames(popcov)
        }
        simuresults <- withCallingHandlers(
          tryCatch({
            wrapper(data = data,
                    model = model_est,
                    latent1 = latent1,
                    latent2 = latent2,
                    alpha = alpha_vec,
                    scale = FALSE,
                    htmt2 = FALSE,
                    nboot = number_bootruns)
          }, error = function(e){
            cat("An error occured:", e$message, "\n")
            list(delta = list(HTMT = NA,
                              se = NA,
                              lowerbound = NA,
                              upperbound = NA,
                              time = NA,
                              missing = NA),
                 boot = list(se = NA,
                             lowerbound = NA,
                             upperbound = NA,
                             missing = NA,
                             time = NA),
                 bcaboot = list(se = NA,
                                lowerbound = NA,
                                upperbound = NA,
                                time = NA,
                                missing = NA),
                 bcboot = list(se = NA,
                               lowerbound = NA,
                               upperbound = NA,
                               missing = NA,
                               time = NA))
          }),
          warning = function(w){
            cat("Warning:", w$message, "\n")
            invokeRestart("muffleWarning")
          }
        )
        
        conphi <- c_phi(model_constrained = model_constrained,
                        model_unconstrained = model_est,
                        data = data
                        )

        # FIX (review pt 3): flatten list-valued warnings/errors to character scalars
        # so every column in conphirow stays atomic and rbind(res, conphirow) below
        # doesn't hit a type mismatch against the atomic NA columns in `row`.
        
        con_warn_flat   <- flat_cond(conphi$con_warning)
        uncon_warn_flat <- flat_cond(conphi$uncon_warning)
        con_err_flat    <- if(is.na(conphi$con_error))   NA_character_ else as.character(conphi$con_error)
        uncon_err_flat  <- if(is.na(conphi$uncon_error)) NA_character_ else as.character(conphi$uncon_error)
        
        res <- foreach(type = c("delta", "boot", "bcaboot", "bcboot"), .combine = "rbind")%:%
                foreach(indexgamma = 1:length(alpha_vec), .combine = "rbind") %do%
                {
                  HTMT <- simuresults$delta$HTMT
                  seHTMT <- simuresults[[type]]$se
                  lowerbound <- unname(simuresults[[type]]$lowerbound[indexgamma])
                  upperbound <- unname(simuresults[[type]]$upperbound[indexgamma])
                  missing <- unname(simuresults[[type]]$missing)
                  time <- simuresults[[type]]$time
                  row <- data.frame(correlation = correlation
                                    , n = n
                                    , HTMT = HTMT
                                    , seHTMT = seHTMT
                                    , simruns = r
                                    , batch = rep_batch
                                    , datatype = datatype
                                    , method = type
                                    , alpha = alpha_vec[indexgamma]
                                    , lowerbound = lowerbound
                                    , upperbound =  upperbound
                                    , coveragecorr = lowerbound <= correlation & correlation <=  upperbound
                                    , coverageone =  lowerbound <= 1 & 1 <= upperbound
                                    , time = time
                                    , seed = seed
                                    , missing = missing
                                    , teststat_constrained = NA_real_
                                    , teststat_unconstrained = NA_real_
                                    , error_constrained = NA_character_
                                    , error_unconstrained = NA_character_
                                    , warning_constrained = NA_character_
                                    , warning_unconstrained = NA_character_
                                    , df_constrained = NA_real_
                                    , df_unconstrained = NA_real_
                  )
                  
                  row
                }
        conphirow <- data.frame(correlation = correlation
                                , n = n
                                , HTMT = NA
                                , seHTMT = NA
                                , simruns = r
                                , batch = rep_batch
                                , datatype = datatype
                                , method = "conphi"
                                , alpha = NA
                                , lowerbound = NA
                                , upperbound =  NA
                                , coveragecorr = NA
                                , coverageone =  NA
                                , time = NA
                                , seed = seed
                                , missing = NA 
                                , teststat_constrained = conphi$con_teststat
                                , teststat_unconstrained = conphi$uncon_teststat
                                , error_constrained = con_err_flat
                                , error_unconstrained = uncon_err_flat
                                , warning_constrained = con_warn_flat
                                , warning_unconstrained = uncon_warn_flat
                                , df_constrained = conphi$con_df
                                , df_unconstrained = conphi$uncon_df
        )
          res <- rbind(res, conphirow)
          res
      }

outfile <- paste0("~/simulationhtmt/results/model_", cond$model, "_sample_size_",
                  cond$n, "_batch_", cond$rep_batch, "_datatype_", data_id, ".rds")
tmpfile <- paste0(outfile, ".tmp")

cat("About to saveRDS to tmp file:", tmpfile, "\n")
success <- FALSE
for(i in 1:3){
  tryCatch(saveRDS(partial, tmpfile, version = 2),
           error = function(e) cat("saveRDS attempt", i, "threw:", conditionMessage(e), "\n"))
  if(file.exists(tmpfile) && file.info(tmpfile)$size > 0){
    success <- TRUE
    break
  }
  cat("attempt", i, "failed, tmpfile missing or empty -- retrying in 5s\n")
  Sys.sleep(5)
}
if(!success) stop("saveRDS failed after 3 attempts: ", tmpfile)

Sys.sleep(3)   # NFS flush buffer before rename
cat("renaming to:", outfile, "\n")
if(!file.rename(tmpfile, outfile)){
  stop("file.rename failed even though tmpfile exists: ", tmpfile)
}


        
        
