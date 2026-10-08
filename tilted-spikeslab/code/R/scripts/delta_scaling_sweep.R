# R/scripts/delta_scaling_sweep.R
#
# Rejection-sampled sweep across delta(q) scaling rules. We compare four rules
# on the same axes:
#   rule "const1":  delta = 1
#   rule "const2":  delta = 2
#   rule "logq":    delta = log(q)
#   rule "linear":  delta = (q - 2) / 2     (G-Wishart-style)
#
# For each (variant, q, rule) we draw from the delta-tilted PD-conditional
# prior using the rejection sampler from R/src/Z_prior_tilt_resample.cpp.
# Cells with const1 and const2 reuse the existing draws from
# R/results/cone_interior_tilt_resample.rds; logq and linear cells are run
# here.
#
# Output: R/results/delta_scaling.rds

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
SRC_CPP <- file.path(PROJECT_ROOT, "R", "src", "Z_prior_tilt_resample.cpp")
EXISTING_RDS <- file.path(PROJECT_ROOT, "R", "results", "cone_interior_tilt_resample.rds")
OUT_RDS <- file.path(PROJECT_ROOT, "R", "results", "delta_scaling.rds")

sourceCpp(SRC_CPP)

set.seed(20260514L)

qs       <- c(5L, 10L, 15L, 20L, 30L, 50L)
variants <- c("discrete", "continuous")

delta_for <- function(rule, q) {
  switch(rule,
         const1 = 1,
         const2 = 2,
         logq   = log(q),
         linear = (q - 2) / 2,
         stop("unknown rule"))
}

p_inc_for_q <- c(`5` = 0.50, `10` = 0.20, `15` = 0.13,
                 `20` = 0.10, `30` = 0.06, `50` = 0.03)

# Per-q proposal budget. Logq and linear rules push delta higher, which shifts
# the proposal toward diagonally-dominant K; PD rate is typically much higher
# and we need fewer proposals than for the const sweeps.
n_proposals_for_q <- c(`5`  = 200000L,   `10` = 500000L,
                       `15` = 1000000L,  `20` = 2000000L,
                       `30` = 5000000L,  `50` = 10000000L)

ALPHA <- 3
BETA  <- 3
V0    <- 0.1
V1    <- 1.0
SIGMA <- 1.0

# Build the cells to run (new rules only).
new_rules <- c("logq", "linear")
cells <- list()
seed_counter <- 1L
for (variant in variants) {
  for (q in qs) {
    for (rule in new_rules) {
      d <- delta_for(rule, q)
      cells[[length(cells) + 1L]] <- list(
        variant = variant, q = q, rule = rule, delta = d,
        sigma = SIGMA,
        n_proposals = n_proposals_for_q[as.character(q)],
        p_inc = p_inc_for_q[as.character(q)],
        seed = 2000L + seed_counter
      )
      seed_counter <- seed_counter + 1L
    }
  }
}
cat(sprintf("New cells to run: %d  (rules: %s)\n",
            length(cells), paste(new_rules, collapse = ", ")))

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
  list(variant = cell$variant, q = cell$q, rule = cell$rule, delta = cell$delta,
       sigma = cell$sigma, n_proposals = cell$n_proposals,
       p_inc = cell$p_inc, n_pd = res$n_pd, n_accepted = res$n_accepted,
       median = med, q05 = q05, q95 = q95, lmin = lmin)
}

t0 <- Sys.time()
new_results <- mclapply(cells, run_cell, mc.cores = 12L, mc.preschedule = FALSE)
cat(sprintf("New sweep done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Pull const1 and const2 cells from the existing dataset.
existing <- readRDS(EXISTING_RDS)
const_results <- list()
for (r in existing$results) {
  if (abs(r$sigma - 1) > 1e-12) next
  if (!isTRUE(r$delta %in% c(1, 2))) next
  rule <- if (r$delta == 1) "const1" else "const2"
  const_results[[length(const_results) + 1L]] <- list(
    variant = r$variant, q = r$q, rule = rule, delta = r$delta,
    sigma = r$sigma, n_proposals = r$n_proposals,
    p_inc = r$p_inc, n_pd = r$n_pd, n_accepted = r$n_accepted,
    median = r$median, q05 = r$q05, q95 = r$q95, lmin = r$lmin
  )
}

all_results <- c(const_results, new_results)
summary_df <- do.call(rbind, lapply(all_results, function(r) {
  data.frame(variant = r$variant, q = r$q, rule = r$rule, delta = r$delta,
             sigma = r$sigma, n_proposals = r$n_proposals,
             n_pd = r$n_pd, n_accepted = r$n_accepted,
             median = r$median, q05 = r$q05, q95 = r$q95,
             stringsAsFactors = FALSE)
}))
saveRDS(list(summary = summary_df, results = all_results), OUT_RDS)
cat(sprintf("Wrote %s\n", OUT_RDS))
print(summary_df[order(summary_df$variant, summary_df$rule, summary_df$q),
                 c("variant","rule","q","delta","n_accepted","median","q05","q95")])
