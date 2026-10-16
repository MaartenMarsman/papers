# Shared diagnostic helpers for split-Rhat and the two-state summary of
# indicator chains (used by R/threshold-study-helpers.R).

# Classic split-Rhat on a [n_iter x n_chains] draws matrix, replicating the
# bgms implementation: halve each chain, pooled W and B over the 2m
# sub-chains. Returns NA when every sub-chain is constant and equal (bgms
# reports NA for constant chains), Inf when sub-chains are constant but
# unequal.
split_rhat <- function(draws) {
  n_iter_local <- nrow(draws)
  half <- n_iter_local %/% 2L
  sub <- cbind(draws[seq_len(half), , drop = FALSE],
               draws[seq.int(n_iter_local - half + 1L, n_iter_local), ,
                     drop = FALSE])
  n <- nrow(sub)
  chain_means <- colMeans(sub)
  chain_vars <- colSums(sub^2)
  chain_vars <- (chain_vars - n * chain_means^2) / (n - 1)
  chain_vars[chain_vars < 0] <- 0
  w <- mean(chain_vars)
  b <- n * stats::var(chain_means)
  if (w == 0) {
    if (b == 0) return(NA_real_)
    return(Inf)
  }
  sqrt(((n - 1) / n * w + b / n) / w)
}

# Two-state Markov chain summary for one indicator chain: Jeffreys-smoothed
# transition probabilities give a PIP estimate and Markov-chain variance that
# remain interior for saturated chains (all 0 or all 1), where the empirical
# variance p(1 - p) and hence the usual ESS collapse to zero. Flips are the
# boundary-robust information measure: well-defined even when zero.
chain_pip_stats <- function(x) {
  n_draws <- length(x)
  x_from <- x[-n_draws]
  x_to <- x[-1L]
  n0 <- sum(x_from == 0L)
  n1 <- n_draws - 1L - n0
  n01 <- sum(x_from == 0L & x_to == 1L)
  n10 <- sum(x_from == 1L & x_to == 0L)
  a_hat <- (n01 + 0.5) / (n0 + 1)
  b_hat <- (n10 + 0.5) / (n1 + 1)
  p_hat <- a_hat / (a_hat + b_hat)
  ess_smooth <- n_draws * (a_hat + b_hat) / (2 - a_hat - b_hat)
  c(p_hat = p_hat, var_p = p_hat * (1 - p_hat) / ess_smooth,
    flips = n01 + n10, ess_smooth = ess_smooth)
}

# Chi-square-style agreement statistic across chains from per-chain smoothed
# PIPs and variances (df = n_chains - 1). A binary analogue of Rhat;
# miscalibrated in the extreme few-flip regime, so read it next
# to the flip counts, not instead of them.
indicator_agreement_stats <- function(per_chain) {
  weights <- 1 / per_chain["var_p", ]
  p_pool <- sum(weights * per_chain["p_hat", ]) / sum(weights)
  sum((per_chain["p_hat", ] - p_pool)^2 / per_chain["var_p", ])
}

# Per-edge diagnostic table for a bgms fit: one row per edge with PIP,
# manual split-Rhat for indicator and unconditional pairwise chains,
# flip counts, per-chain PIP range, smoothed ESS, and agreement statistic.
edge_diagnostics <- function(fit, n_iter, n_chains) {
  indicators <- bgms::extract_indicators(fit)
  pairwise <- bgms::extract_pairwise_interactions(fit)
  stopifnot(nrow(indicators) == n_iter * n_chains)
  chain_id <- rep(seq_len(n_chains), each = n_iter)
  edge_names <- colnames(indicators)
  rows <- vector("list", ncol(indicators))
  for (e in seq_along(rows)) {
    x <- indicators[, e]
    per_chain <- vapply(seq_len(n_chains),
                        function(ch) chain_pip_stats(x[chain_id == ch]),
                        c(p_hat = 0, var_p = 0, flips = 0, ess_smooth = 0))
    rows[[e]] <- data.frame(
      edge = edge_names[e],
      pip = mean(x),
      rhat_indicator = split_rhat(matrix(x, nrow = n_iter)),
      rhat_pairwise = split_rhat(matrix(pairwise[, e], nrow = n_iter)),
      total_flips = sum(per_chain["flips", ]),
      pip_chain_min = min(per_chain["p_hat", ]),
      pip_chain_max = max(per_chain["p_hat", ]),
      ess_smooth_total = sum(per_chain["ess_smooth", ]),
      agreement_stat = indicator_agreement_stats(per_chain)
    )
  }
  do.call(rbind, rows)
}
