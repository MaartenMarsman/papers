# R/scripts/cone_interior_tilt_resample.R
#
# Direct (rejection-sampled) draws from the delta-tilted spike-and-slab prior
# conditional on K positive-definite. Replaces the importance-reweighted tilt
# summary used by the original cone_interior_diagnostics sweep, which loses
# effective sample size as delta grows and leaves gaps in Figure 2.
#
# For each (variant, q, delta) at sigma = 1, draw K from a proposal with
# diagonals Gamma(alpha + delta, beta) and accept with probability
# (det K / prod K_ll)^delta on the PD draws. Accepted draws are i.i.d. from
# the delta-tilted PD-conditional prior. Each accepted draw is one effective
# sample.
#
# Output: R/results/cone_interior_tilt_resample.rds

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
SRC_CPP <- file.path(PROJECT_ROOT, "R", "src", "Z_prior_tilt_resample.cpp")
OUT_RDS <- file.path(PROJECT_ROOT, "R", "results", "cone_interior_tilt_resample.rds")

sourceCpp(SRC_CPP)

set.seed(20260513L)

# Sweep design at sigma = 1.
qs       <- c(5L, 10L, 15L, 20L, 30L, 50L)
deltas   <- c(0, 0.5, 1, 2)
variants <- c("discrete", "continuous")

# Match the q-dependent inclusion probabilities from the original sweep.
p_inc_for_q <- c(`5` = 0.50, `10` = 0.20, `15` = 0.13,
                 `20` = 0.10, `30` = 0.06, `50` = 0.03)

# Per-q proposal budget. High-q cells get more proposals because PD rate is
# astronomically low.
n_proposals_for_q <- c(`5`  = 500000L,   `10` = 1000000L,
                       `15` = 2000000L,  `20` = 5000000L,
                       `30` = 15000000L, `50` = 15000000L)

ALPHA <- 3
BETA  <- 3
V0    <- 0.1
V1    <- 1.0     # continuous slab sd (matches sigma = 1 for discrete)
SIGMA <- 1.0     # discrete slab sd

cells <- list()
seed_counter <- 1L
for (variant in variants) {
  for (q in qs) {
    for (d in deltas) {
      cells[[length(cells) + 1L]] <- list(
        variant = variant,
        q       = q,
        delta   = d,
        sigma   = SIGMA,
        n_proposals = n_proposals_for_q[as.character(q)],
        p_inc   = p_inc_for_q[as.character(q)],
        seed    = 1000L + seed_counter
      )
      seed_counter <- seed_counter + 1L
    }
  }
}
cat(sprintf("Sweep has %d cells (sigma = 1 only).\n", length(cells)))

run_cell <- function(cell) {
  variant_int <- if (cell$variant == "discrete") 0L else 1L
  res <- z_prior_tilt_resample_cpp(
    q = cell$q,
    alpha = ALPHA, beta = BETA, delta = cell$delta,
    variant = variant_int,
    sigma = SIGMA,
    v0 = V0, v1 = V1,
    gamma_prob = cell$p_inc,
    n_proposals = cell$n_proposals,
    seed = cell$seed
  )
  lmin <- res$lmin
  if (length(lmin) >= 1L) {
    med <- median(lmin)
    q05 <- as.numeric(quantile(lmin, 0.05, names = FALSE))
    q95 <- as.numeric(quantile(lmin, 0.95, names = FALSE))
  } else {
    med <- NA_real_; q05 <- NA_real_; q95 <- NA_real_
  }
  list(
    variant   = cell$variant,
    q         = cell$q,
    delta     = cell$delta,
    sigma     = cell$sigma,
    n_proposals = cell$n_proposals,
    p_inc     = cell$p_inc,
    n_pd      = res$n_pd,
    n_accepted = res$n_accepted,
    median    = med,
    q05       = q05,
    q95       = q95,
    lmin      = lmin       # accepted draws of lambda_min
  )
}

# Parallel sweep.
t0 <- Sys.time()
results <- mclapply(cells, run_cell, mc.cores = 12L, mc.preschedule = FALSE)
cat(sprintf("Sweep done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Tidy summary table (drop the raw lmin draws from the table).
summary_df <- do.call(rbind, lapply(results, function(r) {
  data.frame(variant = r$variant, q = r$q, delta = r$delta, sigma = r$sigma,
             n_proposals = r$n_proposals, n_pd = r$n_pd,
             n_accepted = r$n_accepted,
             median = r$median, q05 = r$q05, q95 = r$q95,
             stringsAsFactors = FALSE)
}))

saveRDS(list(summary = summary_df, results = results),
        OUT_RDS)
cat(sprintf("Wrote %s\n", OUT_RDS))
print(summary_df)
