# R/scripts/cone_interior_diagnostics.R
#
# SUPERSEDED by cone_interior_diagnostics_expanded.R (which produces the data of
# figure fig:cone-interior); this version covers q <= 20, sigma = 1 only, and
# feeds only the legacy fig_cone_interior.R.
#
# Demonstrate the cone-interior collapse of the un-modified spike-and-slab
# prior, independent of the Laplace approximation to Z(Gamma). For each
# variant (discrete, continuous) and a sweep over q, samples (K, gamma)
# jointly from the base measure with gamma_ij ~ Bernoulli(gamma_prob), and
# reports:
#
#   1. Pr[K positive-definite] under the base measure (no PD constraint).
#   2. Pr[lambda_min(K) > epsilon] for epsilon in {0.05, 0.1, 0.2}.
#   3. Pr[lambda_min(K) < 0.05 | PD] (boundary concentration on PD-conditional).
#   4. Quantiles of lambda_min on PD draws (effective prior bulk location).
#   5. Effective vs nominal diagonal mean (PD-conditioned distortion).
#
# Parameters fixed at alpha = 3, beta = 3, gamma_prob = 0.5 (matching the
# Z project's joint-sampling diagnostic). The discrete variant uses
# sigma = 1; the continuous variant uses v0 = 0.1, v1 = 1.
#
# Runs in parallel via mclapply on 12 cores.

suppressPackageStartupMessages({
  library(Rcpp)
  library(RcppArmadillo)
  library(parallel)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
CPP_SRC      <- file.path(PROJECT_ROOT, "R", "src", "Z_prior_diag.cpp")
OUT_RDS      <- file.path(PROJECT_ROOT, "R", "results", "cone_interior_diagnostics.rds")

Rcpp::sourceCpp(CPP_SRC)

weighted_quantile <- function(x, w, probs) {
  ord <- order(x)
  x <- x[ord]; w <- w[ord]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1L))
}

run_one_cell <- function(args) {
  q        <- args$q
  variant  <- args$variant
  alpha    <- args$alpha
  beta     <- args$beta
  sigma    <- args$sigma
  v0       <- args$v0
  v1       <- args$v1
  gamma_p  <- args$gamma_prob
  n_draws  <- args$n_draws
  seed     <- args$seed
  delta_grid <- args$delta_grid

  variant_int <- if (variant == "discrete") 0L else 1L
  res <- z_prior_diag_cpp(q, alpha, beta, variant_int,
                          sigma, v0, v1, gamma_p, n_draws, seed)

  lmin    <- res$lmin
  log_det <- res$log_det
  pd      <- res$pd == 1L
  pd_rate <- mean(pd)
  lmin_pd <- lmin[pd]
  log_det_pd <- log_det[pd]
  diag_pd <- res$diag_mean[pd]
  n_pd    <- length(lmin_pd)

  eps_grid <- c(0, 0.05, 0.1, 0.2, 0.5)
  pr_above <- vapply(eps_grid, function(e) mean(lmin > e), numeric(1L))

  # Tilt-only diagnostics: re-weight PD draws by |K|^delta = exp(delta * log_det).
  tilt_summary <- lapply(delta_grid, function(d) {
    if (n_pd < 10L) {
      return(list(delta = d, median = NA_real_, q05 = NA_real_, q95 = NA_real_,
                  pr_lmin_lt_eps = setNames(rep(NA_real_, length(eps_grid)),
                                            sprintf("eps_%g", eps_grid))))
    }
    log_w <- d * log_det_pd
    log_w <- log_w - max(log_w)             # stabilise
    w <- exp(log_w)
    qs <- weighted_quantile(lmin_pd, w, c(0.05, 0.50, 0.95))
    pr_above_w <- vapply(eps_grid, function(e) sum(w * (lmin_pd <= e)) / sum(w), numeric(1L))
    list(delta = d,
         q05    = qs[1L],
         median = qs[2L],
         q95    = qs[3L],
         pr_lmin_le_eps = setNames(pr_above_w, sprintf("eps_%g", eps_grid)))
  })

  list(
    q              = q,
    variant        = variant,
    alpha          = alpha,
    beta           = beta,
    sigma          = sigma,
    v0             = v0,
    v1             = v1,
    gamma_prob     = gamma_p,
    n_draws        = n_draws,
    n_pd           = n_pd,
    pd_rate        = pd_rate,
    eps_grid       = eps_grid,
    pr_above       = pr_above,
    pr_lmin_lt05_given_pd = if (n_pd > 0L) mean(lmin_pd < 0.05) else NA_real_,
    lmin_pd_mean   = if (n_pd > 0L) mean(lmin_pd) else NA_real_,
    lmin_pd_median = if (n_pd > 0L) stats::median(lmin_pd) else NA_real_,
    lmin_pd_q05    = if (n_pd > 10L) stats::quantile(lmin_pd, 0.05, names = FALSE) else NA_real_,
    lmin_pd_q95    = if (n_pd > 10L) stats::quantile(lmin_pd, 0.95, names = FALSE) else NA_real_,
    diag_mean_pd   = if (n_pd > 0L) mean(diag_pd) else NA_real_,
    diag_mean_all  = mean(res$diag_mean),
    tilt_summary   = tilt_summary
  )
}

# ----- driver -----

q_grid       <- c(5L, 10L, 15L, 20L)
variants     <- c("discrete", "continuous")
ALPHA        <- 3.0
BETA         <- 3.0
SIGMA        <- 1.0
V0           <- 0.1
V1           <- 1.0

# Edge inclusion probability scaled with q (sparse-asymptotic regime matching
# the Z project's prior_filter_rate_sweep.R). At q=20 this gives an expected
# |E| ~ 19 edges, comparable to a realistic Bayesian-GGM analysis.
p_inc_for_q  <- c(`5` = 0.50, `10` = 0.20, `15` = 0.13, `20` = 0.10)
n_draws_by_q <- c(`5` = 5e5L,  `10` = 1e6L, `15` = 2e6L,  `20` = 4e6L)
DELTA_GRID   <- c(0, 0.5, 1, 2)

jobs <- list()
k <- 0L
for (q in q_grid) {
  for (variant in variants) {
    k <- k + 1L
    jobs[[length(jobs) + 1L]] <- list(
      q          = q,
      variant    = variant,
      alpha      = ALPHA,
      beta       = BETA,
      sigma      = SIGMA,
      v0         = V0,
      v1         = V1,
      gamma_prob = p_inc_for_q[as.character(q)],
      n_draws    = n_draws_by_q[as.character(q)],
      delta_grid = DELTA_GRID,
      seed       = 42L + 1000L * k
    )
  }
}

cat(sprintf("Running %d cells on 12 cores (gamma_prob scaled with q: %s)...\n",
            length(jobs),
            paste(sprintf("q=%s:%.2f", names(p_inc_for_q), p_inc_for_q), collapse = ", ")))
t0 <- Sys.time()
results <- parallel::mclapply(jobs, run_one_cell, mc.cores = 12L)
t1 <- Sys.time()
cat(sprintf("Done in %.1f s.\n", as.numeric(t1 - t0, units = "secs")))

# Summary table (single row per cell)
tab <- data.frame(
  variant         = vapply(results, function(r) r$variant, character(1L)),
  q               = vapply(results, function(r) r$q, integer(1L)),
  n_draws         = vapply(results, function(r) r$n_draws, integer(1L)),
  n_pd            = vapply(results, function(r) r$n_pd, integer(1L)),
  pd_rate         = vapply(results, function(r) r$pd_rate, numeric(1L)),
  pr_lt05_given_pd = vapply(results, function(r) r$pr_lmin_lt05_given_pd, numeric(1L)),
  lmin_pd_median  = vapply(results, function(r) r$lmin_pd_median, numeric(1L)),
  lmin_pd_q05     = vapply(results, function(r) r$lmin_pd_q05, numeric(1L)),
  lmin_pd_q95     = vapply(results, function(r) r$lmin_pd_q95, numeric(1L)),
  diag_mean_pd    = vapply(results, function(r) r$diag_mean_pd, numeric(1L)),
  diag_mean_all   = vapply(results, function(r) r$diag_mean_all, numeric(1L)),
  stringsAsFactors = FALSE
)
tab <- tab[order(tab$variant, tab$q), ]

# Pr[lmin > eps] wide table (one row per cell, columns per eps)
eps_grid <- c(0, 0.05, 0.1, 0.2, 0.5)
pr_above_mat <- t(vapply(results, function(r) r$pr_above, numeric(length(eps_grid))))
colnames(pr_above_mat) <- sprintf("pr_lmin_above_%g", eps_grid)
pr_tab <- cbind(
  data.frame(
    variant = vapply(results, function(r) r$variant, character(1L)),
    q       = vapply(results, function(r) r$q, integer(1L)),
    stringsAsFactors = FALSE
  ),
  as.data.frame(pr_above_mat)
)
pr_tab <- pr_tab[order(pr_tab$variant, pr_tab$q), ]

cat("\n=== Summary table (unconditional PD rate + boundary concentration) ===\n")
print(tab, row.names = FALSE, digits = 4)

cat("\n=== Pr[lambda_min > epsilon] ===\n")
print(pr_tab, row.names = FALSE, digits = 4)

cat("\n=== Effective vs nominal diagonal mean (PD-conditioned distortion) ===\n")
distort <- data.frame(
  variant         = tab$variant,
  q               = tab$q,
  nominal_mean    = ALPHA / BETA,
  diag_mean_pd    = tab$diag_mean_pd,
  distortion_ratio = tab$diag_mean_pd / (ALPHA / BETA),
  stringsAsFactors = FALSE
)
print(distort, row.names = FALSE, digits = 4)

# Tilt-only: how the median lambda_min shifts as delta grows (no floor)
tilt_rows <- list()
for (r in results) {
  for (ts in r$tilt_summary) {
    tilt_rows[[length(tilt_rows) + 1L]] <- data.frame(
      variant       = r$variant,
      q             = r$q,
      delta         = ts$delta,
      lmin_median   = ts$median,
      lmin_q05      = ts$q05,
      lmin_q95      = ts$q95,
      pr_lmin_le_005 = ts$pr_lmin_le_eps["eps_0.05"],
      pr_lmin_le_01  = ts$pr_lmin_le_eps["eps_0.1"],
      stringsAsFactors = FALSE
    )
  }
}
tilt_tab <- do.call(rbind, tilt_rows)
rownames(tilt_tab) <- NULL

cat("\n=== Tilt-only effect: boundary concentration as a function of delta ===\n")
cat("(weighted by |K|^delta on PD draws; no floor)\n")
print(tilt_tab[order(tilt_tab$variant, tilt_tab$q, tilt_tab$delta), ],
      row.names = FALSE, digits = 4)

saveRDS(list(jobs = jobs, results = results,
             summary = tab, pr_table = pr_tab, distortion = distort,
             tilt_only = tilt_tab),
        OUT_RDS)
cat(sprintf("\nSaved diagnostics to %s\n", OUT_RDS))
