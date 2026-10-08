# c(0)-compensated joint specification, enumeration-based version.
# SUPERSEDED by the compensated arm of joint_spec_sparsening_J1_chain_table.R,
# which the manuscript quotes (the enumeration constants are unreliable at
# delta = 2).
#
# The J1 finding that the realized average per-edge odds factor tracks c(0)
# suggests a zero-cost corrective: state the inclusion odds inflated by 1/c(0),
#   p_c / (1 - p_c) = [p / (1 - p)] / c(0),
# and run the joint specification with p_c. This script recomputes the exact
# enumeration under the compensated prior and reports how closely the intended
# expected edge count is recovered, per replicate seed (to separate MC noise).
#
# Output: R/results/joint_spec_compensation.rds (+ printed table)

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
zt <- readRDS(file.path(PROJECT_ROOT, "R", "results", "joint_spec_sparsening_J1_z.rds"))

c0_fun <- function(delta, eta) {
  s_int <- function(s) {
    t <- (delta + 1) * s
    exp(delta * log(t) - t - lgamma(delta + 1) + log(delta + 1)) /
      sqrt(1 + 2 * eta^2 / t)
  }
  stats::integrate(s_int, 0, Inf, rel.tol = 1e-11)$value
}

rows <- list()
for (q in c(4L, 5L)) {
  deltas <- sort(unique(zt$delta[zt$q == q]))
  for (eta in c(0.5, 1, 2)) for (delta in deltas) {
    c0 <- c0_fun(delta, eta)
    M <- q * (q - 1) / 2
    for (p in c(0.25, 0.5)) {
      odds_c <- (p / (1 - p)) / c0
      pc <- odds_c / (1 + odds_c)
      per_seed <- vapply(sort(unique(zt$seed)), function(s) {
        d <- zt[zt$q == q & zt$eta == eta & abs(zt$delta - delta) < 1e-9 & zt$seed == s, ]
        d <- d[order(d$mask), ]
        pj <- exp(d$m * log(pc) + (M - d$m) * log(1 - pc)) * d$z_hat
        sum(d$m * pj / sum(pj))
      }, numeric(1))
      rows[[length(rows) + 1L]] <- data.frame(
        q = q, eta = eta, delta = delta, p_target = p, c0 = c0, p_comp = pc,
        E_intended = M * p, E_comp = mean(per_seed),
        rel_err = mean(per_seed) / (M * p) - 1,
        rel_err_seed_halfdiff = abs(diff(per_seed)) / 2 / (M * p))
    }
  }
}
res <- do.call(rbind, rows)
saveRDS(res, file.path(PROJECT_ROOT, "R", "results", "joint_spec_compensation.rds"))

cat("== compensated joint specification: recovery of the intended E[m] ==\n")
print(res[order(res$q, res$eta, res$delta, res$p_target), ],
      row.names = FALSE, digits = 3)
cat(sprintf("\nrel err range (all cells): %+.3f to %+.3f; all positive: %s\n",
            min(res$rel_err), max(res$rel_err), all(res$rel_err > 0)))
