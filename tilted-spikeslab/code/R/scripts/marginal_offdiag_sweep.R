# R/scripts/marginal_offdiag_sweep.R
#
# Marginal of an off-diagonal entry theta_{ij} | gamma_{ij} = 1 under the
# tilted spike-and-slab prior. Uses the alpha = 1 G-Wishart-slab perfect
# sampler. For each (q, delta) cell we accept enough joint (Gamma, K) draws
# at gamma_{1,2} = 1 to estimate the marginal moments of K[1,2].
#
# Output: R/results/marginal_offdiag.rds

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
SRC_CPP <- file.path(PROJECT_ROOT, "R", "src", "Z_prior_gwishart_slab_k12.cpp")
OUT_RDS <- file.path(PROJECT_ROOT, "R", "results", "marginal_offdiag.rds")

sourceCpp(SRC_CPP)

set.seed(20260515L)

ALPHA <- 1     # G-Wishart-compatible regime (perfect sampler applies)
BETA  <- 1
SIGMA <- 1.0

qs       <- c(5L, 10L, 15L, 20L, 30L)

# Rerun control (added 2026-10-08). The inner cap MAX_INNER_TRIES is the number
# of proposals per sampled graph before the graph is abandoned; the original
# runs used 1e4, which discarded a biased subset of graphs at q >= 20 (see the
# README, known issue 1). With 1e6 a graph is abandoned only when its slab
# acceptance rate is below about 1e-6, and the kept fraction
# n_graphs_kept / n_graphs_tried is recorded per cell. Set the environment
# variable SPIKESLAB_QS (for example "20,30") to rerun only those dimensions;
# the new cells then replace the matching cells of the archived results file
# and all other cells are kept as archived.
MAX_INNER_TRIES <- 1000000L
QS_RUN <- if (nzchar(Sys.getenv("SPIKESLAB_QS"))) {
  as.integer(strsplit(Sys.getenv("SPIKESLAB_QS"), ",")[[1]])
} else NULL
if (!is.null(QS_RUN)) qs <- qs[qs %in% QS_RUN]
deltas   <- c(0, 1, 2)
p_inc_for_q <- c(`5` = 0.50, `10` = 0.20, `15` = 0.13,
                 `20` = 0.10, `30` = 0.06)

# Per-q proposal budgets, scaled by acceptance. With condition_on_edge_01,
# every accepted draw is on-edge (n_on_edge = n_accepted), so n_target is
# the number of usable on-edge draws directly.
n_proposal_cap_for_q <- c(`5`  = 500000L,    `10` = 2000000L,
                          `15` = 10000000L,  `20` = 200000000L,
                          `30` = 1000000000L)
n_target_for_q       <- c(`5`  = 5000L, `10` = 5000L,
                          `15` = 3000L, `20` = 3000L,
                          `30` = 2000L)

# Draws of K[1,2] conditional on gamma_{1,2} = 1: the recorder forces the edge
# (1,2) into every proposed graph (condition_on_edge_01 = TRUE below) and draws
# the remaining edges iid Bernoulli, so every accepted draw is an exact draw
# from the tilted prior given gamma_{1,2} = 1.
#
# Caveat (2026-10-08 code review): the perfect sampler drops a proposed graph
# when no K is accepted within max_inner_tries = 1e4 proposals; at q = 20 with
# delta >= 1 and at q = 30 this discards a substantial fraction of graphs and
# biases the kept graphs toward high slab acceptance. The q <= 15 cells are
# unaffected. The C++ recorder now reports n_graphs_tried / n_graphs_kept.

cells <- list(); seed_counter <- 1L
for (q in qs) {
  for (d in deltas) {
    cells[[length(cells) + 1L]] <- list(
      q = q, delta = d,
      n_target = n_target_for_q[as.character(q)],
      n_proposal_cap = n_proposal_cap_for_q[as.character(q)],
      p_inc = p_inc_for_q[as.character(q)],
      seed = 4000L + seed_counter
    )
    seed_counter <- seed_counter + 1L
  }
}
cat(sprintf("Cells: %d\n", length(cells)))

run_one <- function(cell) {
  res <- gwishart_slab_record_k12_cpp(
    q = cell$q, delta = cell$delta, beta = BETA, sigma = SIGMA,
    gamma_prob = cell$p_inc,
    n_target = cell$n_target,
    n_proposal_cap = cell$n_proposal_cap,
    max_inner_tries = MAX_INNER_TRIES,
    seed = cell$seed,
    condition_on_edge_01 = TRUE
  )
  # Filter for gamma_{1,2} = 1 draws
  on_edge <- res$gamma12 == 1
  k12_on <- res$k12[on_edge]
  k12_off <- res$k12[!on_edge]
  a_on  <- res$a[on_edge]
  c_on  <- res$c[on_edge]
  m_on  <- res$m[on_edge]
  list(
    q = cell$q, delta = cell$delta, p_inc = cell$p_inc,
    n_proposals = res$n_proposals, n_accepted = res$n_accepted,
    n_graphs_tried = res$n_graphs_tried, n_graphs_kept = res$n_graphs_kept,
    n_inv_fail = res$n_inv_fail, max_inner_tries = MAX_INNER_TRIES,
    n_on_edge = sum(on_edge),
    k12_on_edge = k12_on,
    k12_off_edge = k12_off,
    a_on = a_on, c_on = c_on, m_on = m_on
  )
}

t0 <- Sys.time()
results <- mclapply(cells, run_one, mc.cores = 12L, mc.preschedule = FALSE)
if (!is.null(QS_RUN) && file.exists(OUT_RDS)) {
  old <- readRDS(OUT_RDS)$results
  keep <- Filter(function(r) !(r$q %in% QS_RUN), old)
  cat(sprintf("Merging %d rerun cells with %d archived cells\n", length(results), length(keep)))
  results <- c(keep, results)
  ord <- order(vapply(results, function(r) r$q, numeric(1)),
               vapply(results, function(r) r$delta, numeric(1)))
  results <- results[ord]
}
cat(sprintf("Done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

summary_df <- do.call(rbind, lapply(results, function(r) {
  on <- r$k12_on_edge
  if (length(on) >= 2L) {
    mean_on <- mean(on); var_on <- var(on)
    q05 <- as.numeric(quantile(on, 0.05))
    q95 <- as.numeric(quantile(on, 0.95))
  } else {
    mean_on <- NA_real_; var_on <- NA_real_; q05 <- NA_real_; q95 <- NA_real_
  }
  data.frame(q = r$q, delta = r$delta, p_inc = r$p_inc,
             n_accepted = r$n_accepted, n_on_edge = r$n_on_edge,
             kept_frac = if (is.null(r$n_graphs_tried)) NA_real_ else r$n_graphs_kept / r$n_graphs_tried,
             inv_fail = if (is.null(r$n_inv_fail)) NA_real_ else r$n_inv_fail,
             mean = mean_on, var = var_on, q05 = q05, q95 = q95,
             stringsAsFactors = FALSE)
}))

saveRDS(list(summary = summary_df, results = results,
             params = list(alpha = ALPHA, beta = BETA, sigma = SIGMA)),
        OUT_RDS)
cat(sprintf("Wrote %s\n", OUT_RDS))
print(summary_df)
