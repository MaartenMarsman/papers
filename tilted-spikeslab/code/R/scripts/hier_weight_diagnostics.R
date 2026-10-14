# R/scripts/hier_weight_diagnostics.R
# Purpose: separate the three quantities that Section 5.3 (H4) conflated in the
#          "underflow" description of the hierarchical ratio estimator at scale:
#          (a) the fraction of auxiliary draws whose Cholesky completion overflows
#              (log weight -Inf, weight exactly zero),
#          (b) the fraction of draws with finite log weight at or below -700 (the
#              threshold the archived "usable" fraction used), and
#          (c) the importance-weight effective sample size, ESS = (sum w)^2 / sum w^2
#              with w = exp(lw - max lw) and zero-weight draws counted in the
#              denominator N, reported as ESS / N.
#          Same graphs, edge, eta = 1, delta = 0.5 log q as hier_ratio_accuracy_at_scale_A2.R.
# Run:     Rscript R/scripts/hier_weight_diagnostics.R   (single core, a few minutes)
# Output:  R/results/hier_weight_diagnostics.rds (+ log)

suppressPackageStartupMessages(library(Rcpp))
PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
sourceCpp(file.path(PROJECT_ROOT, "R", "src", "spikeslab_sampler.cpp"))

er_graph <- function(q, seed) {
  set.seed(seed)
  G <- matrix(0L, q, q)
  p <- 4 / (q - 1)
  for (a in 1:(q - 1)) for (b in (a + 1):q)
    if (runif(1) < p) { G[a, b] <- 1L; G[b, a] <- 1L }
  G
}
pick_edge <- function(G) {
  absent <- which(upper.tri(G) & G == 0, arr.ind = TRUE)
  as.integer(absent[nrow(absent) %/% 2, ])
}
diag_one <- function(lw) {
  N <- length(lw); fin <- is.finite(lw)
  w <- numeric(N); w[fin] <- exp(lw[fin] - max(lw[fin]))
  c(n = N, frac_finite = mean(fin), frac_le_m700 = mean(fin & lw <= -700),
    frac_usable_archived = mean(lw > -700),
    med_finite = median(lw[fin]), q10_finite = unname(quantile(lw[fin], 0.1)),
    ess_over_n = (sum(w)^2 / sum(w^2)) / N, max_share = max(w) / sum(w))
}

rows <- list()
for (q in c(10L, 20L, 50L, 100L)) {
  G <- er_graph(q, 4500L + q); e <- pick_edge(G); delta <- 0.5 * log(q)
  for (N in c(1000L, 10000L)) {
    set.seed(7700L + q)
    d <- akm_ratio_cpp(G, e[1], e[2], delta, 1, 1, N, coupled = TRUE, return_draws = TRUE)
    for (arm in c("minus", "plus")) {
      lw <- if (arm == "minus") d$draws_minus else d$draws_plus
      rows[[length(rows) + 1L]] <- data.frame(q = q, N = N, arm = arm, t(diag_one(lw)))
    }
  }
}
tab <- do.call(rbind, rows)
cat("== completion overflow, threshold exclusion, and weight ESS (coupled estimator, eta = 1, delta = 0.5 log q) ==\n")
print(tab, row.names = FALSE, digits = 3)
saveRDS(tab, file.path(PROJECT_ROOT, "R", "results", "hier_weight_diagnostics.rds"))
cat("Done. Results saved to R/results/hier_weight_diagnostics.rds\n")
