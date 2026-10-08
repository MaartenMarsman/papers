# R/scripts/c_sensitivity_sweep.R
#
# Sensitivity of the dimension-adaptive coefficient c in delta = c * log(q)
# to the slab variance sigma and the diagonal-prior shape alpha. For each
# (alpha, sigma) in {1, 3} x {0.3, 1, 3} we run a (q, delta) grid and
# extract c(target) for several target medians of lambda_min.
#
# Sampler choice (correctness identical across samplers):
#   alpha = 1, sigma >= 1: G-Wishart-slab perfect sampler (fast slab accept).
#   alpha = 1, sigma <  1: rejection sampler with PD + tilt rejection.
#   alpha = 3            : rejection sampler (G-Wishart envelope needs alpha=1).
#
# Target: at least 500 accepted draws per cell, capped at n_proposals_max.
#
# Output: R/results/c_sensitivity.rds

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
SRC_GW <- file.path(PROJECT_ROOT, "R", "src", "Z_prior_gwishart_slab.cpp")
SRC_RJ <- file.path(PROJECT_ROOT, "R", "src", "Z_prior_tilt_resample.cpp")
OUT_RDS <- file.path(PROJECT_ROOT, "R", "results", "c_sensitivity.rds")

sourceCpp(SRC_GW)
sourceCpp(SRC_RJ)

set.seed(20260516L)

V0 <- 0.1
N_TARGET <- 500L

# Sweep design
alphas <- c(1, 3)
sigmas <- c(0.3, 1, 3)
deltas <- c(0, 0.5, 1, 2)
qs     <- c(5L, 10L, 15L, 20L, 30L, 50L)

p_inc_for_q <- c(`5` = 0.50, `10` = 0.20, `15` = 0.13,
                 `20` = 0.10, `30` = 0.06, `50` = 0.03)

# Per-cell proposal cap. Bounded by ~5 min wall time at q=50 (eig_sym
# dominates). Cells that don't reach n_target = 500 by their cap are kept
# but flagged for the analysis to drop or downweight.
n_cap_table <- expand.grid(alpha = alphas, sigma = sigmas, q = qs)
n_cap_table$cap <- with(n_cap_table, {
  q_cap <- ifelse(q <= 10, 5e7,
           ifelse(q == 15, 3e7,
           ifelse(q == 20, 1.5e7,
           ifelse(q == 30, 8e6, 3e6))))
  # Hard cap regardless of (alpha, sigma)
  pmin(q_cap, 5e7)
})

cells <- list()
seed_counter <- 1L
for (alpha in alphas) {
  for (sigma in sigmas) {
    for (q in qs) {
      cap <- n_cap_table$cap[n_cap_table$alpha == alpha &
                              n_cap_table$sigma == sigma & n_cap_table$q == q]
      for (d in deltas) {
        cells[[length(cells) + 1L]] <- list(
          alpha = alpha, beta = alpha,   # E[theta_ll] = alpha/beta = 1
          sigma = sigma, q = q, delta = d,
          n_target = N_TARGET,
          n_proposal_cap = as.integer(cap),
          p_inc = p_inc_for_q[as.character(q)],
          seed = 5000L + seed_counter
        )
        seed_counter <- seed_counter + 1L
      }
    }
  }
}
cat(sprintf("Cells: %d\n", length(cells)))

run_cell <- function(cell) {
  use_gwishart <- (cell$alpha == 1) && (cell$sigma >= 1)
  variant_int <- 0L  # discrete

  t0 <- Sys.time()
  if (use_gwishart) {
    res <- gwishart_slab_resample_cpp(
      q = cell$q, delta = cell$delta, beta = cell$beta, sigma = cell$sigma,
      gamma_prob = cell$p_inc,
      n_target = cell$n_target,
      n_proposal_cap = cell$n_proposal_cap,
      max_inner_tries = 10000L,
      seed = cell$seed
    )
    lmin <- res$lmin
    n_pd <- res$n_proposals     # G-Wishart envelope produces only PD draws
    n_accepted <- res$n_accepted
    n_proposals <- res$n_proposals
  } else {
    res <- z_prior_tilt_resample_cpp(
      q = cell$q, alpha = cell$alpha, beta = cell$beta, delta = cell$delta,
      variant = variant_int, sigma = cell$sigma,
      v0 = V0, v1 = cell$sigma,
      gamma_prob = cell$p_inc,
      n_proposals = cell$n_proposal_cap,
      seed = cell$seed
    )
    lmin <- res$lmin
    n_pd <- res$n_pd
    n_accepted <- res$n_accepted
    n_proposals <- res$n_proposals
  }
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  med <- if (length(lmin) >= 5L) median(lmin) else NA_real_
  q05 <- if (length(lmin) >= 20L) as.numeric(quantile(lmin, 0.05)) else NA_real_
  q95 <- if (length(lmin) >= 20L) as.numeric(quantile(lmin, 0.95)) else NA_real_

  list(alpha = cell$alpha, beta = cell$beta, sigma = cell$sigma,
       q = cell$q, delta = cell$delta,
       sampler = if (use_gwishart) "gwishart" else "rejection",
       n_pd = n_pd, n_accepted = n_accepted, n_proposals = n_proposals,
       median = med, q05 = q05, q95 = q95,
       elapsed_s = dt)
}

t0 <- Sys.time()
results <- mclapply(cells, run_cell, mc.cores = 12L, mc.preschedule = FALSE)
cat(sprintf("Sweep done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

summary_df <- do.call(rbind, lapply(results, function(r) {
  data.frame(alpha = r$alpha, beta = r$beta, sigma = r$sigma,
             q = r$q, delta = r$delta, sampler = r$sampler,
             n_pd = r$n_pd, n_accepted = r$n_accepted,
             n_proposals = r$n_proposals,
             median = r$median, q05 = r$q05, q95 = r$q95,
             elapsed_s = r$elapsed_s,
             stringsAsFactors = FALSE)
}))

saveRDS(list(summary = summary_df, results = results),
        OUT_RDS)
cat(sprintf("Wrote %s\n", OUT_RDS))
print(summary_df[order(summary_df$alpha, summary_df$sigma, summary_df$q, summary_df$delta),
                 c("alpha", "sigma", "q", "delta", "sampler",
                   "n_accepted", "n_proposals", "median")])
