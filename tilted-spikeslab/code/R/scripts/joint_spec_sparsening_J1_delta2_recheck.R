# R/scripts/joint_spec_sparsening_J1_delta2_recheck.R
#
# Recheck of the brute-force (importance-sampling) normalizing constants at
# delta = 2, q = 5, with ten times the draws of joint_spec_sparsening_J1.R
# (2e6 instead of 2e5 per graph) and two independent seeds. It quantifies the
# between-seed spread of the enumeration-based expected edge count E_J[m] at
# delta = 2, where the importance weights are heavy-tailed. This spread is the
# reason the manuscript's sparsening table is generated from chain runs
# (joint_spec_sparsening_J1_chain_table.R) rather than from the enumeration.
#
# The original recheck (2026-07-28) was run interactively and only its console
# output survived as R/results/j1_delta2_recheck.log; this script reproduces
# that computation with the same seeds (91001, 91002).
#
# Output: R/results/j1_delta2_recheck.rds (+ console log)
# Runtime: about 15 minutes on 12 cores.

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
sourceCpp(file.path(PROJECT_ROOT, "R", "src", "monotonicity_check.cpp"))
OUT_RDS <- file.path(PROJECT_ROOT, "R", "results", "j1_delta2_recheck.rds")

q       <- 5L
delta   <- 2
etas    <- c(1, 2)
p_incs  <- c(0.25, 0.5)
seeds   <- c(91001L, 91002L)
n_draws <- 2e6

edge_list <- function(q) t(utils::combn(q, 2))
graph_from_mask <- function(mask, edges, q) {
  G <- matrix(0L, q, q)
  for (k in seq_len(nrow(edges))) {
    if (bitwAnd(mask, bitwShiftL(1L, k - 1L)) != 0L) {
      G[edges[k, 1], edges[k, 2]] <- 1L
      G[edges[k, 2], edges[k, 1]] <- 1L
    }
  }
  G
}
edge_count <- function(mask) sum(bitwAnd(mask, bitwShiftL(1L, 0:30)) != 0L)

edges  <- edge_list(q)
M      <- nrow(edges)
masks  <- 0:(2^M - 1L)
m_of   <- vapply(masks, edge_count, integer(1))
graphs <- lapply(masks, graph_from_mask, edges = edges, q = q)

# Expected edge count under p_J(Gamma) proportional to w(Gamma) Z(Gamma).
expected_edges <- function(z_hat, p_inc) {
  logw <- m_of * log(p_inc) + (M - m_of) * log(1 - p_inc)
  pj <- exp(logw) * z_hat
  pj <- pj / sum(pj)
  sum(m_of * pj)
}

rows <- list()
t_all <- Sys.time()
for (eta in etas) for (seed in seeds) {
  t0 <- Sys.time()
  z <- mclapply(seq_along(masks), function(i) {
    r <- z_graph_aligned_cpp(graphs[[i]], beta = eta, sigma = 1,
                             delta = delta, n_draws = n_draws, seed = seed)
    c(r$z_hat, r$z_se)
  }, mc.cores = 12L)
  z <- do.call(rbind, z)
  for (p in p_incs) {
    E_J <- expected_edges(z[, 1], p)
    cat(sprintf("eta=%.0f seed=%d p=%.2f: E_J = %.4f\n", eta, seed, p, E_J))
    rows[[length(rows) + 1L]] <- data.frame(q = q, delta = delta, eta = eta,
                                            seed = seed, p_inc = p, E_J = E_J,
                                            n_draws = n_draws)
  }
  cat(sprintf("   (eta=%.0f seed=%d: %d graphs in %.1f s)\n", eta, seed,
              length(masks), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
res <- do.call(rbind, rows)

# Between-seed half-differences per (eta, p_inc): the Monte Carlo spread.
spread <- do.call(rbind, lapply(split(res, interaction(res$eta, res$p_inc)), function(d) {
  data.frame(eta = d$eta[1], p_inc = d$p_inc[1],
             E_J_seed1 = d$E_J[d$seed == seeds[1]],
             E_J_seed2 = d$E_J[d$seed == seeds[2]],
             half_diff = abs(diff(d$E_J)) / 2)
}))
cat("\nBetween-seed half-differences of E_J[m] at delta = 2, q = 5, 2e6 draws per graph:\n")
print(spread, row.names = FALSE, digits = 4)
cat(sprintf("\nTotal time %.1f min\n", as.numeric(difftime(Sys.time(), t_all, units = "mins"))))

saveRDS(list(results = res, spread = spread,
             settings = list(q = q, delta = delta, etas = etas, p_incs = p_incs,
                             seeds = seeds, n_draws = n_draws)),
        OUT_RDS)
cat(sprintf("Wrote %s\n", OUT_RDS))
