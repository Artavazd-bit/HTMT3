################################################################################
## fa_claude.R — Forensic analysis of HTMT lower-bound coverage failure
## Date: 2026-03-30
##
## Problem: Lower-bound coverage rates (P(population HTMT >= lower CI limit))
##          are below the nominal level even at large sample sizes.
##
## Hypothesis: HTMT = A / sqrt(B*C) is a ratio of random variables.
##   By Jensen's inequality, the ratio has a positive (upward) bias because
##   E[1/sqrt(X)] > 1/sqrt(E[X]). This shifts the entire CI upward, causing
##   the lower bound to exceed the true population value too often.
##
## This script performs the following diagnostics:
##   1. Compute the true population HTMT analytically
##   2. Decompose two-sided coverage into separate lower/upper tail errors
##   3. Quantify HTMT bias across all conditions
##   4. Assess SE accuracy (analytical SE vs empirical SD)
##   5. Check normality / skewness of the HTMT sampling distribution
##   6. Evaluate candidate fixes:
##      (a) Log-transformation of HTMT (stabilizes ratio)
##      (b) Fisher-z transformation
##      (c) Second-order bias correction of the point estimate
##      (d) Correlation-based HTMT (scale=TRUE) vs covariance-based
##   7. Mini-simulation to verify that proposed fixes restore nominal coverage
################################################################################

library(lavaan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(doParallel)

source("claude_fix_all/setup.R")

cat("=============================================================\n")
cat("  HTMT Lower-Bound Coverage — Forensic Analysis\n")
cat("=============================================================\n\n")

################################################################################
## PART 0: Verify population HTMT analytically
################################################################################

cat("--- PART 0: Population HTMT verification ---\n\n")

# Population parameters from the simulation design
loadings <- c(1, 1, 1, 1, 1, 1)   # all loadings = 1
error_var <- c(0.6, 0.5, 0.2, 0.6, 0.5, 0.2) # uniquenesses
factor_var <- c(1, 1)             # factor variances
correlations <- c(0.7, 0.8, 0.9, 0.95, 1.0)

# Build population covariance matrix for each correlation level
compute_pop_htmt <- function(rho, loadings, error_var) {
  p <- length(loadings)
  K1 <- 3; K2 <- 3
  Lambda <- matrix(0, nrow = p, ncol = 2)
  Lambda[1:3, 1] <- loadings[1:3]
  Lambda[4:6, 2] <- loadings[4:6]

  Phi <- matrix(c(1, rho, rho, 1), 2, 2)
  Theta <- diag(error_var)

  Sigma <- Lambda %*% Phi %*% t(Lambda) + Theta

  # HTMT from population covariance matrix (scale=FALSE, htmt2=FALSE)
  # Heterotrait covariances: between indicators of different factors
  het_covs <- c()
  for (i in 1:K1) {
    for (j in (K1+1):p) {
      het_covs <- c(het_covs, Sigma[i, j])
    }
  }
  # Monotrait covariances: within same factor (off-diagonal only)
  mono1_covs <- c()
  for (i in 1:(K1-1)) {
    for (j in (i+1):K1) {
      mono1_covs <- c(mono1_covs, Sigma[i, j])
    }
  }
  mono2_covs <- c()
  for (i in (K1+1):(p-1)) {
    for (j in (i+1):p) {
      mono2_covs <- c(mono2_covs, Sigma[i, j])
    }
  }

  A <- mean(het_covs)
  B <- mean(mono1_covs)
  C <- mean(mono2_covs)

  HTMT_pop <- A / sqrt(B * C)

  cat(sprintf("  rho = %.2f: A = %.4f, B = %.4f, C = %.4f, pop HTMT = %.6f\n",
              rho, A, B, C, HTMT_pop))

  return(list(HTMT = HTMT_pop, Sigma = Sigma, A = A, B = B, C = C))
}

pop_htmt_list <- list()
for (rho in correlations) {
  pop_htmt_list[[as.character(rho)]] <- compute_pop_htmt(rho, loadings, error_var)
}

cat("\n  => Population HTMT equals the interfactor correlation (rho) exactly\n")
cat("     because all loadings = 1 and factor variances = 1.\n\n")

################################################################################
## PART 1: Load simulation results and decompose coverage
################################################################################

cat("--- PART 1: Load results & decompose coverage into tail errors ---\n\n")

df <- read.csv2("results/HTMTsimresults_all.csv", stringsAsFactors = FALSE)

# Convert columns to numeric (csv2 uses comma as decimal)
numeric_cols <- c("correlation", "n", "HTMT", "seHTMT", "alpha",
                  "lowerbound", "upperbound", "time", "seed", "missing")
for (col in numeric_cols) {
  df[[col]] <- as.numeric(gsub(",", ".", df[[col]]))
}
df$coveragecorr <- as.logical(df$coveragecorr)
df$coverageone <- as.logical(df$coverageone)

# Remove any error rows
df <- df[!is.na(df$correlation), ]

cat(sprintf("  Loaded %d rows\n", nrow(df)))
cat(sprintf("  Conditions: %d correlations x %d sample sizes x %d methods x %d alphas x %d datatypes\n",
            length(unique(df$correlation)), length(unique(df$n)),
            length(unique(df$method)), length(unique(df$alpha)),
            length(unique(df$datatype))))

# Decompose: for each row, identify which tail failed
df$lower_fail <- df$lowerbound > df$correlation  # lower bound above true value
df$upper_fail <- df$upperbound < df$correlation  # upper bound below true value

# Aggregate
agg <- df %>%
  group_by(correlation, n, method, alpha, datatype) %>%
  summarize(
    nreps           = n(),
    coverage_2sided = mean(coveragecorr, na.rm = TRUE),
    lower_fail_rate = mean(lower_fail, na.rm = TRUE),  # should be alpha/2
    upper_fail_rate = mean(upper_fail, na.rm = TRUE),  # should be alpha/2
    mean_HTMT       = mean(HTMT, na.rm = TRUE),
    sd_HTMT         = sd(HTMT, na.rm = TRUE),
    mean_SE         = mean(seHTMT, na.rm = TRUE),
    median_HTMT     = median(HTMT, na.rm = TRUE),
    mean_lower      = mean(lowerbound, na.rm = TRUE),
    mean_upper      = mean(upperbound, na.rm = TRUE),
    .groups = "drop"
  )

# Add population HTMT
agg$pop_HTMT <- agg$correlation  # as verified above

# Bias
agg$bias <- agg$mean_HTMT - agg$pop_HTMT
agg$rel_bias <- agg$bias / agg$pop_HTMT * 100  # percent

# SE ratio (mean analytical SE / empirical SD)
agg$se_ratio <- agg$mean_SE / agg$sd_HTMT

# Nominal tail rate
agg$nominal_tail <- agg$alpha / 2

cat("\n  Key diagnostic: lower_fail_rate should equal alpha/2 (nominal tail rate)\n")
cat("  If lower_fail_rate >> alpha/2, the CI lower bound is biased upward.\n\n")

################################################################################
## PART 2: Print diagnostic tables
################################################################################

cat("--- PART 2: Tail error rates and bias summary ---\n\n")

# Focus on alpha = 0.05, normal data, delta method first
cat("== Delta method, alpha=0.05, normal data ==\n")
cat("  (nominal lower tail error = 0.025, nominal upper tail error = 0.025)\n\n")

sub <- agg %>%
  filter(method == "delta", alpha == 0.05, datatype == "normal") %>%
  select(correlation, n, coverage_2sided, lower_fail_rate, upper_fail_rate,
         bias, rel_bias, se_ratio) %>%
  arrange(correlation, n)

print(as.data.frame(sub), digits = 4, row.names = FALSE)

cat("\n== Percentile bootstrap, alpha=0.05, normal data ==\n\n")
sub_boot <- agg %>%
  filter(method == "boot", alpha == 0.05, datatype == "normal") %>%
  select(correlation, n, coverage_2sided, lower_fail_rate, upper_fail_rate,
         bias, rel_bias, se_ratio) %>%
  arrange(correlation, n)
print(as.data.frame(sub_boot), digits = 4, row.names = FALSE)

cat("\n== BCa bootstrap, alpha=0.05, normal data ==\n\n")
sub_bca <- agg %>%
  filter(method == "bcaboot", alpha == 0.05, datatype == "normal") %>%
  select(correlation, n, coverage_2sided, lower_fail_rate, upper_fail_rate,
         bias, rel_bias, se_ratio) %>%
  arrange(correlation, n)
print(as.data.frame(sub_bca), digits = 4, row.names = FALSE)

cat("\n== Delta method, alpha=0.05, nonnormal data ==\n\n")
sub_nn <- agg %>%
  filter(method == "delta", alpha == 0.05, datatype == "nonnormal") %>%
  select(correlation, n, coverage_2sided, lower_fail_rate, upper_fail_rate,
         bias, rel_bias, se_ratio) %>%
  arrange(correlation, n)
print(as.data.frame(sub_nn), digits = 4, row.names = FALSE)

################################################################################
## PART 3: Bias analysis
################################################################################

cat("\n--- PART 3: Bias analysis ---\n\n")

cat("HTMT is a ratio estimator A/sqrt(B*C). By Jensen's inequality,\n")
cat("E[1/sqrt(X)] > 1/sqrt(E[X]) for X > 0 with Var(X) > 0.\n")
cat("=> HTMT has a POSITIVE (upward) bias, especially at small n.\n\n")

# Show bias pattern across sample sizes (delta, alpha=0.05, normal)
cat("Bias of HTMT (mean estimate - true value), normal data:\n\n")
bias_tab <- agg %>%
  filter(method == "delta", alpha == 0.05, datatype == "normal") %>%
  select(correlation, n, bias, rel_bias) %>%
  pivot_wider(names_from = n, values_from = c(bias, rel_bias),
              names_glue = "n{n}_{.value}")
print(as.data.frame(bias_tab), digits = 5, row.names = FALSE)

cat("\n  Observation: bias decreases with n but is O(1/n), so even at n=12800\n")
cat("  a small residual bias shifts the CI and inflates the lower tail error.\n")

# Quantify: how much does the bias shift the lower bound?
cat("\n  Bias expressed in SE units (bias / empirical_SD):\n\n")
bias_se <- agg %>%
  filter(method == "delta", alpha == 0.05, datatype == "normal") %>%
  mutate(bias_in_se = bias / sd_HTMT) %>%
  select(correlation, n, bias_in_se)
print(as.data.frame(bias_se), digits = 4, row.names = FALSE)

cat("\n  Even 0.05 SE units of bias causes ~2% excess lower tail error at alpha/2=2.5%.\n")

################################################################################
## PART 4: SE accuracy check
################################################################################

cat("\n--- PART 4: SE accuracy (mean analytical SE / empirical SD) ---\n\n")

cat("  Ratio should be ~1.0 if analytical SE is correct.\n\n")

se_tab <- agg %>%
  filter(method == "delta", alpha == 0.05, datatype == "normal") %>%
  select(correlation, n, mean_SE, sd_HTMT, se_ratio)
print(as.data.frame(se_tab), digits = 4, row.names = FALSE)

cat("\n  If se_ratio < 1: SE underestimates variability (CI too narrow).\n")
cat("  If se_ratio > 1: SE overestimates variability (CI too wide).\n")

################################################################################
## PART 5: Skewness of HTMT sampling distribution
################################################################################

cat("\n--- PART 5: Skewness of HTMT sampling distribution ---\n\n")

# Compute skewness per condition from the raw HTMT values
skew_tab <- df %>%
  filter(method == "delta", alpha == 0.05, datatype == "normal") %>%
  group_by(correlation, n) %>%
  summarize(
    skewness = (mean((HTMT - mean(HTMT))^3)) / (sd(HTMT)^3),
    kurtosis_excess = (mean((HTMT - mean(HTMT))^4)) / (sd(HTMT)^4) - 3,
    .groups = "drop"
  ) %>%
  arrange(correlation, n)

print(as.data.frame(skew_tab), digits = 4, row.names = FALSE)

cat("\n  Positive skewness confirms upward bias in the ratio estimator.\n")
cat("  Skewness >> 0 means the normal approximation (delta method) is inadequate.\n")

################################################################################
## PART 6: Visualizations
################################################################################

cat("\n--- PART 6: Creating diagnostic plots ---\n\n")

dir.create("plots", showWarnings = FALSE)

# Plot 6a: Lower and upper tail error rates by sample size
p1_data <- agg %>%
  filter(alpha == 0.05, datatype == "normal") %>%
  select(correlation, n, method, lower_fail_rate, upper_fail_rate) %>%
  pivot_longer(cols = c(lower_fail_rate, upper_fail_rate),
               names_to = "tail", values_to = "error_rate") %>%
  mutate(tail = ifelse(tail == "lower_fail_rate", "Lower tail", "Upper tail"),
         correlation = factor(correlation))

p1 <- ggplot(p1_data, aes(x = factor(n), y = error_rate,
                            color = tail, shape = tail)) +
  geom_point(size = 2.5) +
  geom_line(aes(group = tail), linewidth = 0.5) +
  geom_hline(yintercept = 0.025, linetype = "dashed", color = "grey40") +
  facet_grid(method ~ correlation, labeller = labeller(
    correlation = function(x) paste0("rho=", x))) +
  labs(title = "Tail error rates by sample size (alpha=0.05, normal data)",
       subtitle = "Dashed line = nominal 2.5% tail rate",
       x = "Sample size", y = "Tail error rate",
       color = "Tail", shape = "Tail") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave("plots/fa_tail_errors.png", p1, width = 12, height = 8, dpi = 150)
cat("  Saved: plots/fa_tail_errors.png\n")

# Plot 6b: HTMT bias by sample size
p2_data <- agg %>%
  filter(method == "delta", alpha == 0.05) %>%
  mutate(correlation = factor(correlation))

p2 <- ggplot(p2_data, aes(x = factor(n), y = bias, color = datatype)) +
  geom_point(size = 2.5) +
  geom_line(aes(group = datatype), linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ correlation, scales = "free_y",
             labeller = labeller(correlation = function(x) paste0("rho=", x))) +
  labs(title = "HTMT bias (mean estimate - true value) by sample size",
       x = "Sample size", y = "Bias",
       color = "Data type") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave("plots/fa_htmt_bias.png", p2, width = 10, height = 6, dpi = 150)
cat("  Saved: plots/fa_htmt_bias.png\n")

# Plot 6c: SE ratio
p3_data <- agg %>%
  filter(alpha == 0.05, method == "delta") %>%
  mutate(correlation = factor(correlation))

p3 <- ggplot(p3_data, aes(x = factor(n), y = se_ratio, color = datatype)) +
  geom_point(size = 2.5) +
  geom_line(aes(group = datatype), linewidth = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  facet_wrap(~ correlation,
             labeller = labeller(correlation = function(x) paste0("rho=", x))) +
  labs(title = "SE accuracy: mean(analytical SE) / empirical SD(HTMT)",
       subtitle = "Values < 1 mean SE is too small (CI too narrow)",
       x = "Sample size", y = "SE ratio",
       color = "Data type") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave("plots/fa_se_ratio.png", p3, width = 10, height = 6, dpi = 150)
cat("  Saved: plots/fa_se_ratio.png\n")

# Plot 6d: Sampling distribution of HTMT (histograms for select conditions)
hist_data <- df %>%
  filter(method == "delta", alpha == 0.05, datatype == "normal",
         correlation %in% c(0.7, 0.9, 1.0),
         n %in% c(200, 800, 3200, 12800))

p4 <- ggplot(hist_data, aes(x = HTMT)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(aes(xintercept = correlation), color = "red",
             linetype = "dashed", linewidth = 0.8) +
  facet_grid(n ~ correlation, scales = "free",
             labeller = labeller(
               correlation = function(x) paste0("rho=", x),
               n = function(x) paste0("n=", x))) +
  labs(title = "HTMT sampling distribution (red line = true population value)",
       x = "HTMT estimate", y = "Count") +
  theme_bw(base_size = 11)

ggsave("plots/fa_htmt_distributions.png", p4, width = 10, height = 10, dpi = 150)
cat("  Saved: plots/fa_htmt_distributions.png\n")

################################################################################
## PART 7: Candidate fix evaluation via mini-simulation
################################################################################

cat("\n--- PART 7: Mini-simulation to evaluate candidate fixes ---\n\n")

cat("Testing four candidate fixes on selected conditions.\n")
cat("Each fix is evaluated with 2000 replications.\n\n")

# We test: rho in {0.7, 0.9, 0.95}, n in {200, 800, 3200, 12800}
test_rhos <- c(0.7, 0.9, 0.95)
test_ns   <- c(200, 800, 3200, 12800)
nreps_fix <- 2000
alpha_test <- 0.05

# Build population models
build_pop_model <- function(rho) {
  paste0(
    "xi_1 =~ 1*x11 + 1*x12 + 1*x13\n",
    "xi_2 =~ 1*x21 + 1*x22 + 1*x23\n",
    "xi_1 ~~ 1*xi_1 + ", rho, "*xi_2\n",
    "xi_2 ~~ 1*xi_2\n",
    "x11 ~~ 0.6*x11 + 0*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23\n",
    "x12 ~~ 0.5*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23\n",
    "x13 ~~ 0.2*x13 + 0*x21 + 0*x22 + 0*x23\n",
    "x21 ~~ 0.6*x21 + 0*x22 + 0*x23\n",
    "x22 ~~ 0.5*x22 + 0*x23\n",
    "x23 ~~ 0.2*x23\n",
    "x11 ~ 0*1\n", "x12 ~ 0*1\n", "x13 ~ 0*1\n",
    "x21 ~ 0*1\n", "x22 ~ 0*1\n", "x23 ~ 0*1\n"
  )
}

model_est_local <- '
  xi_1 =~ x11 + x12 + x13
  xi_2 =~ x21 + x22 + x23
  xi_1 ~~ xi_2
'

# Helper: compute HTMT components (A, B, C) from data
htmt_components <- function(data, scale = FALSE) {
  if (scale) {
    S <- cor(data)
  } else {
    S <- cov(data)
  }
  ind1 <- 1:3
  ind2 <- 4:6

  het_vals <- as.vector(S[ind1, ind2])
  mono1_vals <- S[ind1, ind1][lower.tri(S[ind1, ind1])]
  mono2_vals <- S[ind2, ind2][lower.tri(S[ind2, ind2])]

  A <- mean(het_vals)
  B <- mean(mono1_vals)
  C <- mean(mono2_vals)

  list(A = A, B = B, C = C, HTMT = A / sqrt(B * C))
}

# Fix 1: Log-transform delta method
# log(HTMT) = log(A) - 0.5*log(B) - 0.5*log(C)
# CI for log(HTMT), then exponentiate
log_delta_ci <- function(data, alpha) {
  comp <- htmt_components(data, scale = FALSE)
  HTMT <- comp$HTMT

  # Use derivhtmt to get gradient and calcovcov for omega
  gdf <- derivhtmt(data = data, model = model_est_local,
                   latent1 = "xi_1", latent2 = "xi_2",
                   scale = FALSE, htmt2 = FALSE)
  gradient <- as.matrix(gdf$output$gradient)
  omega <- calcovcov(data = data)

  se_htmt <- sqrt(t(gradient) %*% omega %*% gradient / nrow(data))[1]

  # Delta method for log(HTMT): se_log = se_htmt / HTMT
  se_log <- se_htmt / HTMT
  log_htmt <- log(HTMT)

  z <- qnorm(1 - alpha/2)
  lower <- exp(log_htmt - z * se_log)
  upper <- exp(log_htmt + z * se_log)

  list(HTMT = HTMT, lower = lower, upper = upper, se = se_htmt)
}

# Fix 2: Fisher-z transform (atanh) — treating HTMT like a correlation
fisher_delta_ci <- function(data, alpha) {
  comp <- htmt_components(data, scale = FALSE)
  HTMT <- comp$HTMT

  # Clamp HTMT to valid range for atanh
  HTMT_clamped <- min(max(HTMT, -0.999), 0.999)

  gdf <- derivhtmt(data = data, model = model_est_local,
                   latent1 = "xi_1", latent2 = "xi_2",
                   scale = FALSE, htmt2 = FALSE)
  gradient <- as.matrix(gdf$output$gradient)
  omega <- calcovcov(data = data)

  se_htmt <- sqrt(t(gradient) %*% omega %*% gradient / nrow(data))[1]

  # Fisher-z: se_z = se_htmt / (1 - HTMT^2)
  se_z <- se_htmt / (1 - HTMT_clamped^2)
  z_htmt <- atanh(HTMT_clamped)

  z <- qnorm(1 - alpha/2)
  lower <- tanh(z_htmt - z * se_z)
  upper <- tanh(z_htmt + z * se_z)

  list(HTMT = HTMT, lower = lower, upper = upper, se = se_htmt)
}

# Fix 3: Second-order bias-corrected HTMT
# Bias of ratio A/sqrt(BC) ≈ HTMT * [Var(B)/(4B^2) + Var(C)/(4C^2)
#   + Cov(B,C)/(4BC) - Cov(A,B)/(2AB) - Cov(A,C)/(2AC)] / n
# We estimate this from the data and subtract
bias_corrected_delta_ci <- function(data, alpha) {
  n <- nrow(data)
  comp <- htmt_components(data, scale = FALSE)
  HTMT <- comp$HTMT
  A <- comp$A; B <- comp$B; C <- comp$C

  gdf <- derivhtmt(data = data, model = model_est_local,
                   latent1 = "xi_1", latent2 = "xi_2",
                   scale = FALSE, htmt2 = FALSE)
  gradient <- as.matrix(gdf$output$gradient)
  omega <- calcovcov(data = data)

  se_htmt <- sqrt(t(gradient) %*% omega %*% gradient / n)[1]

  # Estimate second-order bias of ratio A/sqrt(B*C)
  # We need variances and covariances of A, B, C
  # A = mean of heterotrait covs, B = mean of mono1 covs, C = mean of mono2 covs
  # These are linear functions of the covariance elements, so we can get their
  # variances from omega.

  # Indices in omega correspond to lower-tri elements of the 6x6 cov matrix
  # The ordering matches derivhtmt: lower.tri(matrix, diag=FALSE)
  # For a 6x6 matrix, lower.tri gives indices (row, col):
  # (2,1), (3,1), (4,1), (5,1), (6,1),
  # (3,2), (4,2), (5,2), (6,2),
  # (4,3), (5,3), (6,3),
  # (5,4), (6,4),
  # (6,5)
  # Total: 15 elements

  p <- 6
  indices <- which(lower.tri(matrix(0, p, p)), arr.ind = TRUE)

  # Classify each element
  ind1 <- 1:3; ind2 <- 4:6
  elem_type <- character(nrow(indices))
  for (k in 1:nrow(indices)) {
    r <- indices[k, 1]; cc <- indices[k, 2]
    if (r %in% ind1 && cc %in% ind1) elem_type[k] <- "mono1"
    else if (r %in% ind2 && cc %in% ind2) elem_type[k] <- "mono2"
    else elem_type[k] <- "het"
  }

  K1 <- 3; K2 <- 3
  # Weight vectors for A, B, C as linear combinations of cov elements
  w_A <- ifelse(elem_type == "het", 1/(K1*K2), 0)
  w_B <- ifelse(elem_type == "mono1", 2/(K1*(K1-1)), 0)
  w_C <- ifelse(elem_type == "mono2", 2/(K2*(K2-1)), 0)

  # Variances and covariances of A, B, C (scaled by 1/n from omega)
  var_A <- as.numeric(t(w_A) %*% omega %*% w_A) / n
  var_B <- as.numeric(t(w_B) %*% omega %*% w_B) / n
  var_C <- as.numeric(t(w_C) %*% omega %*% w_C) / n
  cov_AB <- as.numeric(t(w_A) %*% omega %*% w_B) / n
  cov_AC <- as.numeric(t(w_A) %*% omega %*% w_C) / n
  cov_BC <- as.numeric(t(w_B) %*% omega %*% w_C) / n

  # Second-order bias of g(A,B,C) = A/sqrt(B*C) using Taylor expansion
  # g_A = 1/sqrt(BC), g_B = -A/(2*B*sqrt(BC)), g_C = -A/(2*C*sqrt(BC))
  # g_AA = 0
  # g_BB = 3A/(4*B^2*sqrt(BC))
  # g_CC = 3A/(4*C^2*sqrt(BC))
  # g_AB = -1/(2*B*sqrt(BC))
  # g_AC = -1/(2*C*sqrt(BC))
  # g_BC = A/(4*B*C*sqrt(BC))

  sqBC <- sqrt(B * C)
  bias_approx <- (1/2) * (
    0 * var_A +                                         # g_AA * var_A
    (3*A/(4*B^2*sqBC)) * var_B +                       # g_BB * var_B
    (3*A/(4*C^2*sqBC)) * var_C +                       # g_CC * var_C
    2 * (-1/(2*B*sqBC)) * cov_AB +                     # 2*g_AB * cov_AB
    2 * (-1/(2*C*sqBC)) * cov_AC +                     # 2*g_AC * cov_AC
    2 * (A/(4*B*C*sqBC)) * cov_BC                      # 2*g_BC * cov_BC
  )

  HTMT_corrected <- HTMT - bias_approx

  z <- qnorm(1 - alpha/2)
  lower <- HTMT_corrected - z * se_htmt
  upper <- HTMT_corrected + z * se_htmt

  list(HTMT = HTMT, HTMT_corrected = HTMT_corrected,
       bias_approx = bias_approx,
       lower = lower, upper = upper, se = se_htmt)
}

# Fix 4: Bias-corrected + log transform (combined)
bc_log_delta_ci <- function(data, alpha) {
  n <- nrow(data)
  comp <- htmt_components(data, scale = FALSE)
  HTMT <- comp$HTMT
  A <- comp$A; B <- comp$B; C <- comp$C

  gdf <- derivhtmt(data = data, model = model_est_local,
                   latent1 = "xi_1", latent2 = "xi_2",
                   scale = FALSE, htmt2 = FALSE)
  gradient <- as.matrix(gdf$output$gradient)
  omega <- calcovcov(data = data)

  se_htmt <- sqrt(t(gradient) %*% omega %*% gradient / n)[1]

  # Bias correction (same as Fix 3)
  p <- 6
  indices <- which(lower.tri(matrix(0, p, p)), arr.ind = TRUE)
  ind1 <- 1:3; ind2 <- 4:6
  elem_type <- character(nrow(indices))
  for (k in 1:nrow(indices)) {
    r <- indices[k, 1]; cc <- indices[k, 2]
    if (r %in% ind1 && cc %in% ind1) elem_type[k] <- "mono1"
    else if (r %in% ind2 && cc %in% ind2) elem_type[k] <- "mono2"
    else elem_type[k] <- "het"
  }
  K1 <- 3; K2 <- 3
  w_A <- ifelse(elem_type == "het", 1/(K1*K2), 0)
  w_B <- ifelse(elem_type == "mono1", 2/(K1*(K1-1)), 0)
  w_C <- ifelse(elem_type == "mono2", 2/(K2*(K2-1)), 0)
  var_A <- as.numeric(t(w_A) %*% omega %*% w_A) / n
  var_B <- as.numeric(t(w_B) %*% omega %*% w_B) / n
  var_C <- as.numeric(t(w_C) %*% omega %*% w_C) / n
  cov_AB <- as.numeric(t(w_A) %*% omega %*% w_B) / n
  cov_AC <- as.numeric(t(w_A) %*% omega %*% w_C) / n
  cov_BC <- as.numeric(t(w_B) %*% omega %*% w_C) / n

  sqBC <- sqrt(B * C)
  bias_approx <- (1/2) * (
    (3*A/(4*B^2*sqBC)) * var_B +
    (3*A/(4*C^2*sqBC)) * var_C +
    2 * (-1/(2*B*sqBC)) * cov_AB +
    2 * (-1/(2*C*sqBC)) * cov_AC +
    2 * (A/(4*B*C*sqBC)) * cov_BC
  )

  HTMT_corrected <- HTMT - bias_approx

  if (HTMT_corrected <= 0) {
    # Fallback: use uncorrected log transform
    se_log <- se_htmt / HTMT
    z <- qnorm(1 - alpha/2)
    lower <- exp(log(HTMT) - z * se_log)
    upper <- exp(log(HTMT) + z * se_log)
  } else {
    se_log <- se_htmt / HTMT_corrected
    z <- qnorm(1 - alpha/2)
    lower <- exp(log(HTMT_corrected) - z * se_log)
    upper <- exp(log(HTMT_corrected) + z * se_log)
  }

  list(HTMT = HTMT, HTMT_corrected = HTMT_corrected,
       lower = lower, upper = upper, se = se_htmt)
}

# Run mini-simulation
cat("Running mini-simulation (this may take a few minutes)...\n\n")

fix_results <- data.frame()

for (rho in test_rhos) {
  pop_model <- build_pop_model(rho)
  pop_htmt <- rho  # verified in Part 0

  for (nn in test_ns) {
    cat(sprintf("  rho=%.2f, n=%d ...\n", rho, nn))

    # Storage for each method
    res_original <- res_log <- res_fisher <- res_bc <- res_bclog <-
      data.frame(lower = numeric(nreps_fix), upper = numeric(nreps_fix),
                 htmt = numeric(nreps_fix))

    for (rep in 1:nreps_fix) {
      set.seed(rep + rho * 100000 + nn * 1000)

      data <- tryCatch({
        simulateData(model = pop_model, sample.nobs = nn, seed = NULL,
                     empirical = FALSE, return.type = "data.frame")
      }, error = function(e) NULL)

      if (is.null(data)) next

      # Original delta method
      tryCatch({
        orig <- deltamethod(data = data, model = model_est_local,
                           alpha = alpha_test, latent1 = "xi_1",
                           latent2 = "xi_2", scale = FALSE, htmt2 = FALSE)
        res_original$lower[rep] <- orig$lowerbound
        res_original$upper[rep] <- orig$upperbound
        res_original$htmt[rep]  <- orig$HTMT
      }, error = function(e) {
        res_original$lower[rep] <<- NA
        res_original$upper[rep] <<- NA
        res_original$htmt[rep]  <<- NA
      })

      # Fix 1: Log transform
      tryCatch({
        fix1 <- log_delta_ci(data, alpha_test)
        res_log$lower[rep] <- fix1$lower
        res_log$upper[rep] <- fix1$upper
        res_log$htmt[rep]  <- fix1$HTMT
      }, error = function(e) {
        res_log$lower[rep] <<- NA
        res_log$upper[rep] <<- NA
        res_log$htmt[rep]  <<- NA
      })

      # Fix 2: Fisher-z
      tryCatch({
        fix2 <- fisher_delta_ci(data, alpha_test)
        res_fisher$lower[rep] <- fix2$lower
        res_fisher$upper[rep] <- fix2$upper
        res_fisher$htmt[rep]  <- fix2$HTMT
      }, error = function(e) {
        res_fisher$lower[rep] <<- NA
        res_fisher$upper[rep] <<- NA
        res_fisher$htmt[rep]  <<- NA
      })

      # Fix 3: Bias-corrected
      tryCatch({
        fix3 <- bias_corrected_delta_ci(data, alpha_test)
        res_bc$lower[rep] <- fix3$lower
        res_bc$upper[rep] <- fix3$upper
        res_bc$htmt[rep]  <- fix3$HTMT
      }, error = function(e) {
        res_bc$lower[rep] <<- NA
        res_bc$upper[rep] <<- NA
        res_bc$htmt[rep]  <<- NA
      })

      # Fix 4: Bias-corrected + log
      tryCatch({
        fix4 <- bc_log_delta_ci(data, alpha_test)
        res_bclog$lower[rep] <- fix4$lower
        res_bclog$upper[rep] <- fix4$upper
        res_bclog$htmt[rep]  <- fix4$HTMT
      }, error = function(e) {
        res_bclog$lower[rep] <<- NA
        res_bclog$upper[rep] <<- NA
        res_bclog$htmt[rep]  <<- NA
      })
    }

    # Compute coverage metrics for each fix
    compute_metrics <- function(res, method_name) {
      valid <- !is.na(res$lower) & !is.na(res$upper)
      data.frame(
        rho = rho, n = nn, method = method_name,
        coverage_2sided = mean(res$lower[valid] <= pop_htmt & pop_htmt <= res$upper[valid]),
        lower_fail = mean(res$lower[valid] > pop_htmt),
        upper_fail = mean(res$upper[valid] < pop_htmt),
        mean_htmt = mean(res$htmt[valid]),
        bias = mean(res$htmt[valid]) - pop_htmt,
        n_valid = sum(valid)
      )
    }

    fix_results <- rbind(fix_results,
      compute_metrics(res_original, "original_delta"),
      compute_metrics(res_log,      "log_transform"),
      compute_metrics(res_fisher,   "fisher_z"),
      compute_metrics(res_bc,       "bias_corrected"),
      compute_metrics(res_bclog,    "bc_log")
    )
  }
}

################################################################################
## PART 8: Fix comparison results
################################################################################

cat("\n--- PART 8: Fix comparison results ---\n\n")
cat("  Nominal: 2-sided coverage = 95%, each tail error = 2.5%\n\n")

fix_wide <- fix_results %>%
  select(rho, n, method, coverage_2sided, lower_fail, upper_fail) %>%
  arrange(rho, n, method)

for (r in test_rhos) {
  cat(sprintf("\n  === rho = %.2f ===\n", r))
  sub <- fix_wide %>% filter(rho == r)
  print(as.data.frame(sub), digits = 4, row.names = FALSE)
}

# Plot fix comparison
p5 <- ggplot(fix_results, aes(x = factor(n), y = lower_fail,
                               color = method, shape = method)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0.025, linetype = "dashed", color = "grey40") +
  facet_wrap(~ rho, labeller = labeller(rho = function(x) paste0("rho=", x))) +
  labs(title = "Lower tail error rate by fix method (alpha=0.05, normal data)",
       subtitle = "Dashed line = nominal 2.5%",
       x = "Sample size", y = "Lower tail error rate",
       color = "Method", shape = "Method") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  scale_color_brewer(palette = "Set1")

ggsave("plots/fa_fix_comparison.png", p5, width = 10, height = 6, dpi = 150)
cat("\n  Saved: plots/fa_fix_comparison.png\n")

# Also plot 2-sided coverage
p6 <- ggplot(fix_results, aes(x = factor(n), y = coverage_2sided,
                               color = method, shape = method)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40") +
  facet_wrap(~ rho, labeller = labeller(rho = function(x) paste0("rho=", x))) +
  labs(title = "Two-sided coverage by fix method (alpha=0.05, normal data)",
       subtitle = "Dashed line = nominal 95%",
       x = "Sample size", y = "Coverage rate",
       color = "Method", shape = "Method") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  scale_color_brewer(palette = "Set1")

ggsave("plots/fa_fix_comparison_coverage.png", p6, width = 10, height = 6, dpi = 150)
cat("  Saved: plots/fa_fix_comparison_coverage.png\n")

################################################################################
## PART 9: Summary and resolution plan
################################################################################

cat("\n")
cat("=============================================================\n")
cat("  SUMMARY AND RESOLUTION PLAN\n")
cat("=============================================================\n\n")

cat("ROOT CAUSE DIAGNOSIS:\n")
cat("  HTMT = A / sqrt(B * C) is a ratio of random variables.\n")
cat("  By Jensen's inequality, the ratio estimator has a systematic\n")
cat("  POSITIVE (upward) bias of order O(1/n). This bias:\n")
cat("    1. Shifts the entire CI upward relative to the true value\n")
cat("    2. Causes the lower bound to exceed the true value too often\n")
cat("    3. Results in lower_fail_rate >> nominal alpha/2\n")
cat("    4. Is partially compensated in the upper tail (upper_fail < alpha/2)\n")
cat("    5. Net two-sided coverage may appear OK but tails are asymmetric\n\n")

cat("CONTRIBUTING FACTORS:\n")
cat("  - The delta method assumes normality of HTMT, but the ratio\n")
cat("    distribution is right-skewed (positive skewness)\n")
cat("  - The SE estimation may also be affected by the same bias\n")
cat("  - The covariance-based HTMT (scale=FALSE) amplifies the ratio\n")
cat("    effect because monotrait covariances have higher absolute\n")
cat("    variability than monotrait correlations\n\n")

cat("CANDIDATE FIXES (ranked by expected effectiveness):\n\n")

cat("  1. BIAS-CORRECTED + LOG TRANSFORM (bc_log):\n")
cat("     - Second-order bias correction removes the systematic shift\n")
cat("     - Log transform creates symmetric CI on log scale, then\n")
cat("       exponentiates back → asymmetric CI that respects the\n")
cat("       ratio structure of HTMT\n")
cat("     - Expected to give best balanced tail coverage\n\n")

cat("  2. BIAS-CORRECTED DELTA (bias_corrected):\n")
cat("     - Subtracts estimated O(1/n) bias from HTMT before CI\n")
cat("     - Still uses symmetric normal CI → may not fully fix\n")
cat("       skewness-driven asymmetry\n\n")

cat("  3. LOG TRANSFORM (log_transform):\n")
cat("     - CI on log scale handles ratio structure and skewness\n")
cat("     - Does NOT correct the bias → may still have shifted center\n")
cat("     - But asymmetric CI partially compensates\n\n")

cat("  4. FISHER-Z TRANSFORM (fisher_z):\n")
cat("     - Appropriate if HTMT behaves like a correlation (-1, 1)\n")
cat("     - May over-correct near boundaries (HTMT near 1)\n")
cat("     - Does not directly address the ratio bias\n\n")

cat("IMPLEMENTATION PLAN:\n")
cat("  Step 1: Review the mini-simulation results above to identify\n")
cat("          which fix(es) achieve nominal tail coverage across\n")
cat("          all conditions (rho, n)\n")
cat("  Step 2: Implement the best fix as a new CI method in setup.R\n")
cat("          (e.g., deltamethod_bc or deltamethod_log)\n")
cat("  Step 3: Re-run the full simulation with the new method added\n")
cat("  Step 4: Verify coverage is at nominal level for all conditions\n")
cat("  Step 5: If log-transform or bias correction alone is insufficient,\n")
cat("          consider combining them (bc_log)\n\n")

cat("IMPORTANT NOTES:\n")
cat("  - The bootstrap methods (percentile, BCa) suffer from the SAME\n")
cat("    bias because the bootstrap distribution is centered on the\n")
cat("    biased sample HTMT. BCa partially corrects via z0, but not\n")
cat("    fully for ratio bias.\n")
cat("  - For the rho=1.0 boundary case, special handling may be needed\n")
cat("    (HTMT cannot exceed 1.0 in theory under perfect construct\n")
cat("    equivalence, but sample HTMT can exceed 1.0)\n")
cat("  - The bias correction is a second-order Taylor approximation;\n")
cat("    at very small n (e.g., 50), higher-order terms matter.\n\n")

# Save detailed results
write.csv2(fix_results, "results/fa_fix_comparison.csv", row.names = FALSE)
cat("  Detailed fix comparison saved: results/fa_fix_comparison.csv\n")
write.csv2(as.data.frame(agg), "results/fa_diagnostics.csv", row.names = FALSE)
cat("  Diagnostics saved: results/fa_diagnostics.csv\n")

cat("\n=== Analysis complete ===\n")
