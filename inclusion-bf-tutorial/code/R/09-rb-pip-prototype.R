# Rao-Blackwellized posterior inclusion probabilities: prototype study.
#
# bgms updates each edge indicator with a joint birth-death Metropolis move
# (omrf_model.cpp, update_edge_indicator): propose toggling the indicator
# together with a weight, accept with probability alpha = min(1, exp(la)).
# The one-step Rao-Blackwellized PIP estimator averages the EXPECTED next
# state instead of the realized one,
#     J_t = gamma_t + (1 - 2 gamma_t) * alpha_t,
# using a quantity the sampler already computes and discards. A chain
# saturated at inclusion still evaluates death moves, so J_t < 1: the
# estimator is interior at the boundary, yields finite inclusion Bayes
# factors without epsilon guards, and its chain is CONTINUOUS, so ordinary
# MCSE/ESS machinery applies to it.
#
# This prototype replicates the bgms move structure on a conjugate
# two-predictor spike-and-slab regression (correlated predictors), where
# the exact posterior inclusion probability follows from enumerating the
# four models, and full Rao-Blackwellization (the conditional inclusion
# probability given the other coefficient) is available in closed form.
# Estimators compared, per regime and replication:
#   raw       plain average of the 0/1 indicator draws
#   rb1       one-step Rao-Blackwell (acceptance probabilities)
#   rbfull    full Rao-Blackwell (closed-form conditional PIP)
#   smooth    Jeffreys-smoothed two-state stationary probability (post hoc)
#
# Output: tables/rb-pip-prototype.csv

set.seed(123)

n_obs <- 100
sigma2 <- 1
tau2 <- 1
rho_x <- 0.7
n_chains <- 4L
n_rep <- 300L
n_iter <- 2000L
sd_prop <- 0.4

# Fixed design; scale columns so sufficient statistics are simple.
x1 <- stats::rnorm(n_obs)
x2 <- rho_x * x1 + sqrt(1 - rho_x^2) * stats::rnorm(n_obs)
g11 <- sum(x1^2); g22 <- sum(x2^2); g12 <- sum(x1 * x2)

# Exact log marginal likelihood of model S (subset of {1, 2}) under the
# conjugate slab beta_S ~ N(0, tau2 I): y ~ N(0, sigma2 I + tau2 X_S X_S').
log_marginal <- function(y, cols) {
  X <- cbind(x1, x2)[, cols, drop = FALSE]
  V <- diag(sigma2, n_obs) + tau2 * X %*% t(X)
  as.numeric(mvtnorm::dmvnorm(y, sigma = V, log = TRUE))
}

true_pip1 <- function(y) {
  lm_ <- c(none = log_marginal(y, integer(0)), m1 = log_marginal(y, 1),
           m2 = log_marginal(y, 2), m12 = log_marginal(y, 1:2))
  w <- exp(lm_ - max(lm_))
  (w["m1"] + w["m12"]) / sum(w)
}

# Regimes: coefficients chosen so the true PIP of coefficient 1 spans the
# bulk, the threshold zone, and the boundary.
regimes <- list(
  bulk = c(b1 = 0.10, b2 = 0.15),
  threshold = c(b1 = 0.17, b2 = 0.15),
  decisive = c(b1 = 0.35, b2 = 0.15),
  boundary = c(b1 = 0.42, b2 = 0.15),
  deep_boundary = c(b1 = 0.50, b2 = 0.15)
)

# Closed-form conditional log inclusion odds for coefficient j given the
# other coefficient's current contribution: conjugate scalar update.
cond_log_odds <- function(s_adj, gjj) {
  A <- gjj / sigma2 + 1 / tau2
  b <- s_adj / sigma2
  0.5 * b^2 / A + 0.5 * log(1 / (tau2 * A))
}

run_regime <- function(regime_name, betas) {
  y <- as.numeric(cbind(x1, x2) %*% betas + stats::rnorm(n_obs))
  s1 <- sum(x1 * y); s2 <- sum(x2 * y)
  truth <- true_pip1(y)

  k <- n_rep * n_chains
  gam1 <- stats::rbinom(k, 1L, 0.5); beta1 <- rep(0, k)
  gam2 <- stats::rbinom(k, 1L, 0.5); beta2 <- rep(0, k)
  beta1[gam1 == 1L] <- stats::rnorm(sum(gam1), 0, sqrt(tau2))
  beta2[gam2 == 1L] <- stats::rnorm(sum(gam2), 0, sqrt(tau2))

  sum_raw <- sum_rb1 <- sum_rbfull <- rep(0, k)
  raw_draws <- matrix(0L, nrow = k, ncol = n_iter)

  # Vectorized log-likelihood difference for moving beta_j: old -> new,
  # holding the other coefficient fixed. Sufficient statistics only.
  ll_diff_1 <- function(new1, old1) {
    d <- new1 - old1
    (d * s1 - d * g12 * gam2 * beta2 - 0.5 * (new1^2 - old1^2) * g11) / sigma2
  }
  ll_diff_2 <- function(new2, old2) {
    d <- new2 - old2
    (d * s2 - d * g12 * gam1 * beta1 - 0.5 * (new2^2 - old2^2) * g22) / sigma2
  }

  for (t in seq_len(n_iter)) {
    ## --- coefficient 1: birth-death toggle (bgms move structure) ---------
    prop <- ifelse(gam1 == 0L, stats::rnorm(k, beta1, sd_prop), 0)
    la <- ll_diff_1(prop, beta1)
    slab_term <- stats::dnorm(ifelse(gam1 == 0L, prop, beta1), 0,
                              sqrt(tau2), log = TRUE)
    prop_term <- stats::dnorm(ifelse(gam1 == 0L, prop, beta1),
                              0, sd_prop, log = TRUE)
    la <- la + ifelse(gam1 == 0L, slab_term - prop_term,
                      prop_term - slab_term)
    alpha <- pmin(1, exp(la))
    sum_rb1 <- sum_rb1 + gam1 + (1 - 2 * gam1) * alpha
    flip <- stats::runif(k) < alpha
    beta1 <- ifelse(flip, prop, beta1)
    gam1 <- ifelse(flip, 1L - gam1, gam1)

    ## within-model random walk when included
    inc <- gam1 == 1L
    prop_b <- beta1 + stats::rnorm(k, 0, sd_prop)
    la_b <- ll_diff_1(prop_b, beta1) +
      stats::dnorm(prop_b, 0, sqrt(tau2), log = TRUE) -
      stats::dnorm(beta1, 0, sqrt(tau2), log = TRUE)
    acc_b <- inc & (log(stats::runif(k)) < la_b)
    beta1[acc_b] <- prop_b[acc_b]

    ## full Rao-Blackwell: conditional PIP given (gam2, beta2)
    lo <- cond_log_odds(s1 - g12 * gam2 * beta2, g11)
    sum_rbfull <- sum_rbfull + 1 / (1 + exp(-lo))

    sum_raw <- sum_raw + gam1
    raw_draws[, t] <- gam1

    ## --- coefficient 2: same moves (nuisance) ----------------------------
    prop <- ifelse(gam2 == 0L, stats::rnorm(k, beta2, sd_prop), 0)
    la <- ll_diff_2(prop, beta2)
    slab_term <- stats::dnorm(ifelse(gam2 == 0L, prop, beta2), 0,
                              sqrt(tau2), log = TRUE)
    prop_term <- stats::dnorm(ifelse(gam2 == 0L, prop, beta2),
                              0, sd_prop, log = TRUE)
    la <- la + ifelse(gam2 == 0L, slab_term - prop_term,
                      prop_term - slab_term)
    flip <- stats::runif(k) < pmin(1, exp(la))
    beta2 <- ifelse(flip, prop, beta2)
    gam2 <- ifelse(flip, 1L - gam2, gam2)
    inc <- gam2 == 1L
    prop_b <- beta2 + stats::rnorm(k, 0, sd_prop)
    la_b <- ll_diff_2(prop_b, beta2) +
      stats::dnorm(prop_b, 0, sqrt(tau2), log = TRUE) -
      stats::dnorm(beta2, 0, sqrt(tau2), log = TRUE)
    acc_b <- inc & (log(stats::runif(k)) < la_b)
    beta2[acc_b] <- prop_b[acc_b]
  }

  ## per-replication estimates (pool the 4 chains of each replication)
  rep_id <- rep(seq_len(n_rep), each = n_chains)
  est_raw <- tapply(sum_raw / n_iter, rep_id, mean)
  est_rb1 <- tapply(sum_rb1 / n_iter, rep_id, mean)
  est_rbfull <- tapply(sum_rbfull / n_iter, rep_id, mean)

  est_smooth <- vapply(seq_len(n_rep), function(r) {
    rows <- seq.int((r - 1L) * n_chains + 1L, r * n_chains)
    n01 <- n10 <- n0 <- n1 <- 0
    for (row in rows) {
      x <- raw_draws[row, ]
      xf <- x[-n_iter]; xt <- x[-1L]
      n0 <- n0 + sum(xf == 0L); n1 <- n1 + sum(xf == 1L)
      n01 <- n01 + sum(xf == 0L & xt == 1L)
      n10 <- n10 + sum(xf == 1L & xt == 0L)
    }
    a_hat <- (n01 + 0.5) / (n0 + 1); b_hat <- (n10 + 0.5) / (n1 + 1)
    a_hat / (a_hat + b_hat)
  }, numeric(1))

  ## among replications where the raw average saturates (exactly 1), what
  ## Bayes factor does each remaining estimator report, against the truth?
  sat <- est_raw == 1
  log10bf <- function(p) log10(p / (1 - p))
  data.frame(
    regime = regime_name, true_pip = truth,
    true_log10bf = log10bf(truth),
    prop_raw_saturated = mean(sat),
    rmse_raw = sqrt(mean((est_raw - truth)^2)),
    rmse_rb1 = sqrt(mean((est_rb1 - truth)^2)),
    rmse_rbfull = sqrt(mean((est_rbfull - truth)^2)),
    rmse_smooth = sqrt(mean((est_smooth - truth)^2)),
    med_log10bf_rb1_sat = if (any(sat)) stats::median(log10bf(est_rb1[sat]))
                          else NA_real_,
    med_log10bf_rbfull_sat = if (any(sat))
                               stats::median(log10bf(est_rbfull[sat]))
                             else NA_real_,
    med_log10bf_smooth_sat = if (any(sat))
                               stats::median(log10bf(est_smooth[sat]))
                             else NA_real_,
    sd_log10bf_rb1 = stats::sd(log10bf(pmin(est_rb1, 1 - 1e-15))),
    prop_rb1_saturated = mean(est_rb1 >= 1 - 1e-15)
  )
}

results <- NULL
for (nm in names(regimes)) {
  res <- run_regime(nm, regimes[[nm]])
  results <- rbind(results, res)
  cat(sprintf("%s: true PIP %.6f done\n", nm, res$true_pip))
}
utils::write.csv(results, "tables/rb-pip-prototype.csv", row.names = FALSE)
print(results, digits = 4)
