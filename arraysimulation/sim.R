.libPaths(c("/home/jab49wd/R-projects/R/library", .libPaths()))

library(lavaan)
library(dplyr)
library(foreach)
library(covsim)

source("code/setup.R")
conditions <- readRDS("code/conditions.rds")
n_batches <- 10 

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
  datatype <- "moderatenonnormal"
  
  fit <- sem(current_model, do.fit = FALSE)
  popcov <- lavInspect(fit, "implied")$cov
  popskew <- c(rep(2, 6))
  popkurt <- c(rep(7,6))
}else if(data_id == 3){
  datatype <- "severenonnormal"
  
  fit <- sem(current_model, do.fit = FALSE)
  popcov <- lavInspect(fit, "implied")$cov
  popskew <- c(rep(3,6))
  popkurt<-  c(rep(21,6))
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
        }else if(data_id != 1){
          set.seed(seed)
          data <- covsim::rPLSIM(n, sigma.target = popcov, 
                                 skewness=popskew, excesskurtosis=popkurt,
                                 numsegments = 8)
          data <- as.data.frame(data[[1]][[1]])
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
        
        testdiffstart <- Sys.time()
        conphi <- c_phi(model_constrained = model_constrained,
                        model_unconstrained = model_est,
                        data = data
                        )
        conphipvalue <- 1 - pchisq(conphi$con_teststat - conphi$uncon_teststat,
                                   df = conphi$con_df - conphi$uncon_df)

        uncon_valid <- inherits(conphi$unconfit, "lavaan")
        if (uncon_valid) {
          parest_cfa  <- parameterestimates(conphi$unconfit)
          cor_est_cfa <- parest_cfa$est[parest_cfa$lhs == latent1 & parest_cfa$rhs == latent2]
          cor_est_se  <- parest_cfa$se[parest_cfa$lhs == latent1 & parest_cfa$rhs == latent2]
        } else {
          cor_est_cfa <- NA_real_
          cor_est_se  <- NA_real_
        }
        testdiffend <- Sys.time()
        testdiff <- as.numeric(testdiffend - testdiffstart, units = "secs")
        
        con_warn_flat   <- flat_cond(conphi$con_warning)
        uncon_warn_flat <- flat_cond(conphi$uncon_warning)
        con_err_flat    <- if(is.na(conphi$con_error))   NA_character_ else as.character(conphi$con_error)
        uncon_err_flat  <- if(is.na(conphi$uncon_error)) NA_character_ else as.character(conphi$uncon_error)
        
        res <- foreach(indexgamma = 1:length(alpha_vec), .combine = "rbind") %do% 
                {
                  alpha <- alpha_vec[indexgamma]
                  if (uncon_valid) {
                    parest     <- parameterestimates(conphi$unconfit, level = 1 - alpha)
                    lowerbound <- parest$ci.lower[parest$lhs == latent1 & parest$rhs == latent2]
                    upperbound <- parest$ci.upper[parest$lhs == latent1 & parest$rhs == latent2]
                  } else {
                    lowerbound <- NA_real_
                    upperbound <- NA_real_
                  }
                  
                  conphirow <- data.frame(correlation = correlation
                                          , n = n
                                          , cor_est = cor_est_cfa
                                          , se_cor_est = cor_est_se
                                          , simruns = r
                                          , batch = rep_batch
                                          , datatype = datatype
                                          , method = "conphi"
                                          , alpha = alpha
                                          , lowerbound = lowerbound
                                          , upperbound =  upperbound
                                          , coveragecorr = lowerbound <= correlation & correlation <=  upperbound
                                          , coverageone =  lowerbound <= 1 & 1 <= upperbound
                                          , time = testdiff
                                          , seed = seed
                                          , missing = NA 
                                          , teststat_constrained = conphi$con_teststat
                                          , teststat_unconstrained = conphi$uncon_teststat
                                          , p_value = conphipvalue
                                          , error_constrained = con_err_flat
                                          , error_unconstrained = uncon_err_flat
                                          , warning_constrained = con_warn_flat
                                          , warning_unconstrained = uncon_warn_flat
                                          , df_constrained = conphi$con_df
                                          , df_unconstrained = conphi$uncon_df
                  )
                  
                  htmt <- foreach(type = c("delta", "boot", "bcaboot", "bcboot"), .combine = "rbind")%do%
                    {
                      HTMT <- simuresults$delta$HTMT
                      seHTMT <- simuresults[[type]]$se
                      lowerbound <- unname(simuresults[[type]]$lowerbound[indexgamma])
                      upperbound <- unname(simuresults[[type]]$upperbound[indexgamma])
                      missing <- unname(simuresults[[type]]$missing)
                      time <- as.numeric(simuresults[[type]]$time, units = "secs")
                      row <- data.frame(correlation = correlation
                                        , n = n
                                        , cor_est = HTMT
                                        , se_cor_est = seHTMT
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
                                        , p_value = NA_real_
                                        , error_constrained = NA_character_
                                        , error_unconstrained = NA_character_
                                        , warning_constrained = NA_character_
                                        , warning_unconstrained = NA_character_
                                        , df_constrained = NA_real_
                                        , df_unconstrained = NA_real_
                      )
                      row
                    }
                  rbind(conphirow, htmt)
                }
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

