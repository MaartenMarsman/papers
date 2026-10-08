# R/scripts/delta_scaling_gwishart.R
#
# Scaling analysis on non-empty graphs at alpha = 1 using the G-Wishart-slab
# perfect rejection sampler. For each (q, rule) we draw joint (Gamma, K) from
# the conditional spike-and-slab spec with Gamma ~ Bern(p_inc(q)) and
# K | Gamma from the delta-tilted exponential-diagonal slab prior.
#
# Rules: const0, const1, const2, logq, linear ( (q - 2) / 2, G-Wishart-style ).
# Output: R/results/delta_scaling_gwishart.rds

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
SRC_CPP <- file.path(PROJECT_ROOT, "R", "src", "Z_prior_gwishart_slab.cpp")
OUT_RDS <- file.path(PROJECT_ROOT, "R", "results", "delta_scaling_gwishart.rds")

sourceCpp(SRC_CPP)

set.seed(20260514L)

ALPHA <- 1     # required regime for the G-Wishart envelope rejection
BETA  <- 1     # E[theta_ll] = (1 + delta) / beta, ranges 1..4 across delta
SIGMA <- 1.0

qs       <- c(5L, 10L, 15L, 20L, 30L, 50L)
p_inc_for_q <- c(`5` = 0.50, `10` = 0.20, `15` = 0.13,
                 `20` = 0.10, `30` = 0.06, `50` = 0.03)

# Per-q proposal budgets, scaled by expected acceptance.
n_proposal_cap_for_q <- c(`5`  = 100000L,    `10` = 500000L,
                          `15` = 2000000L,   `20` = 5000000L,
                          `30` = 30000000L,  `50` = 50000000L)
n_target_for_q       <- c(`5`  = 1000L,      `10` = 1000L,
                          `15` = 1000L,      `20` = 1000L,
                          `30` = 500L,       `50` = 200L)
max_inner_tries <- 10000L

delta_rules <- list(
  const0 = function(q) 0,
  const1 = function(q) 1,
  const2 = function(q) 2,
  logq   = function(q) log(q),
  linear = function(q) (q - 2) / 2
)

cells <- list()
seed_counter <- 1L
for (q in qs) {
  for (rule_name in names(delta_rules)) {
    d <- delta_rules[[rule_name]](q)
    cells[[length(cells) + 1L]] <- list(
      q = q, rule = rule_name, delta = d,
      n_target = n_target_for_q[as.character(q)],
      n_proposal_cap = n_proposal_cap_for_q[as.character(q)],
      p_inc = p_inc_for_q[as.character(q)],
      seed = 3000L + seed_counter
    )
    seed_counter <- seed_counter + 1L
  }
}
cat(sprintf("Cells to run: %d\n", length(cells)))

run_cell <- function(cell) {
  res <- gwishart_slab_resample_cpp(
    q = cell$q,
    delta = cell$delta,
    beta = BETA, sigma = SIGMA,
    gamma_prob = cell$p_inc,
    n_target = cell$n_target,
    n_proposal_cap = cell$n_proposal_cap,
    max_inner_tries = max_inner_tries,
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
  list(q = cell$q, rule = cell$rule, delta = cell$delta,
       n_target = cell$n_target, n_proposal_cap = cell$n_proposal_cap,
       p_inc = cell$p_inc,
       n_proposals = res$n_proposals, n_accepted = res$n_accepted,
       n_graphs_tried = res$n_graphs_tried, n_graphs_kept = res$n_graphs_kept,
       median = med, q05 = q05, q95 = q95,
       mean_inner_tries = if (length(res$inner_tries) > 0L) mean(res$inner_tries) else NA_real_,
       mean_n_edges = if (length(res$n_edges) > 0L) mean(res$n_edges) else NA_real_,
       lmin = lmin)
}

t0 <- Sys.time()
results <- mclapply(cells, run_cell, mc.cores = 12L, mc.preschedule = FALSE)
cat(sprintf("Sweep done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

summary_df <- do.call(rbind, lapply(results, function(r) {
  data.frame(q = r$q, rule = r$rule, delta = r$delta,
             n_proposals = r$n_proposals, n_accepted = r$n_accepted,
             n_graphs_tried = r$n_graphs_tried, n_graphs_kept = r$n_graphs_kept,
             mean_inner_tries = r$mean_inner_tries,
             mean_n_edges = r$mean_n_edges,
             median = r$median, q05 = r$q05, q95 = r$q95,
             stringsAsFactors = FALSE)
}))
saveRDS(list(summary = summary_df, results = results,
             params = list(alpha = ALPHA, beta = BETA, sigma = SIGMA)),
        OUT_RDS)
cat(sprintf("Wrote %s\n", OUT_RDS))
print(summary_df[, c("q","rule","delta","n_proposals","n_accepted","median","q05","q95",
                     "mean_n_edges","mean_inner_tries")])
