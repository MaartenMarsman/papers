# R/scripts/partial_correlation_sweep.R
#
# Induced prior on the partial correlation rho_{ij} = -K[i,j]/sqrt(K[i,i]K[j,j])
# under the tilted spike-and-slab prior, on gamma_{ij} = 1 draws. This is the
# "prior you get" on the interpretable dependence parameter, to be contrasted
# with the "prior you specify" (the N(0, sigma^2) slab on K[i,j]).
#
# Uses the alpha = 1 G-Wishart-slab perfect sampler with diagonals recorded
# (Z_prior_gwishart_slab_rho.cpp). Standardized frame: sigma = beta = 1,
# so the geometry knob eta = sigma * beta = 1.
#
# For each (q, delta) cell we accept joint (Gamma, K) draws with gamma_{1,2}=1
# and record (k11, k22, k12), from which rho_{1,2} is formed.
#
# Output: R/results/partial_correlation.rds
#
# Run: source("R/scripts/partial_correlation_sweep.R")  (long-running)

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
SRC_CPP <- file.path(PROJECT_ROOT, "R", "src", "Z_prior_gwishart_slab_rho.cpp")
OUT_RDS <- file.path(PROJECT_ROOT, "R", "results", "partial_correlation.rds")

sourceCpp(SRC_CPP)

set.seed(20260724L)

ALPHA <- 1     # G-Wishart-compatible regime (perfect sampler applies)
BETA  <- 1
SIGMA <- 1.0

qs <- c(5L, 10L, 15L, 20L, 30L)

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

# Fixed tilts {0, 1, 2} plus the dimension-adaptive default delta = 0.5*log q.
# Keeping delta = 0 (untilted) and the adaptive value in every q lets the figure
# contrast "drift with q" (delta = 0) against "stabilizes with q" (adaptive).
delta_grid_for_q <- function(q) {
  sort(unique(round(c(0, 1, 2, 0.5 * log(q)), 3)))
}

p_inc_for_q <- c(`5` = 0.50, `10` = 0.20, `15` = 0.13,
                 `20` = 0.10, `30` = 0.06)

n_proposal_cap_for_q <- c(`5`  = 500000L,    `10` = 2000000L,
                          `15` = 10000000L,  `20` = 200000000L,
                          `30` = 1000000000L)
n_target_for_q       <- c(`5`  = 5000L, `10` = 5000L,
                          `15` = 3000L, `20` = 3000L,
                          `30` = 2000L)

cells <- list(); seed_counter <- 1L
for (q in qs) {
  for (d in delta_grid_for_q(q)) {
    cells[[length(cells) + 1L]] <- list(
      q = q, delta = d,
      n_target = n_target_for_q[as.character(q)],
      n_proposal_cap = n_proposal_cap_for_q[as.character(q)],
      p_inc = p_inc_for_q[as.character(q)],
      seed = 7000L + seed_counter
    )
    seed_counter <- seed_counter + 1L
  }
}
cat(sprintf("Cells: %d\n", length(cells)))

run_one <- function(cell) {
  res <- gwishart_slab_record_rho_cpp(
    q = cell$q, delta = cell$delta, beta = BETA, sigma = SIGMA,
    gamma_prob = cell$p_inc,
    n_target = cell$n_target,
    n_proposal_cap = cell$n_proposal_cap,
    max_inner_tries = MAX_INNER_TRIES,
    seed = cell$seed,
    condition_on_edge_01 = TRUE
  )
  on_edge <- res$gamma12 == 1
  k11 <- res$k11[on_edge]
  k22 <- res$k22[on_edge]
  k12 <- res$k12[on_edge]
  rho <- -k12 / sqrt(k11 * k22)
  ok  <- is.finite(rho)
  list(
    q = cell$q, delta = cell$delta, p_inc = cell$p_inc,
    n_proposals = res$n_proposals, n_accepted = res$n_accepted,
    n_graphs_tried = res$n_graphs_tried, n_graphs_kept = res$n_graphs_kept,
    n_inv_fail = res$n_inv_fail, max_inner_tries = MAX_INNER_TRIES,
    n_on_edge = sum(on_edge),
    rho = rho[ok],
    k12 = k12[ok], k11 = k11[ok], k22 = k22[ok]
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
  rho <- r$rho
  if (length(rho) >= 2L) {
    data.frame(
      q = r$q, delta = r$delta, p_inc = r$p_inc,
      n = length(rho),
      kept_frac = if (is.null(r$n_graphs_tried)) NA_real_ else r$n_graphs_kept / r$n_graphs_tried,
      inv_fail = if (is.null(r$n_inv_fail)) NA_real_ else r$n_inv_fail,
      rho_sd   = sd(rho),
      rho_mean = mean(rho),
      k12_sd   = sd(r$k12),
      p_abs_gt_0.3 = mean(abs(rho) > 0.3),
      p_abs_gt_0.5 = mean(abs(rho) > 0.5),
      q05 = as.numeric(quantile(rho, 0.05)),
      q95 = as.numeric(quantile(rho, 0.95)),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(q = r$q, delta = r$delta, p_inc = r$p_inc, n = length(rho),
               kept_frac = NA, inv_fail = NA, rho_sd = NA, rho_mean = NA, k12_sd = NA,
               p_abs_gt_0.3 = NA, p_abs_gt_0.5 = NA, q05 = NA, q95 = NA,
               stringsAsFactors = FALSE)
  }
}))

saveRDS(list(summary = summary_df, results = results,
             params = list(alpha = ALPHA, beta = BETA, sigma = SIGMA)),
        OUT_RDS)
cat(sprintf("Wrote %s\n", OUT_RDS))
print(summary_df, row.names = FALSE)
