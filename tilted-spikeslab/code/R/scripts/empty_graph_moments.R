# R/scripts/empty_graph_moments.R
#
# SUPERSEDED (pre-pivot, 2026-05): empty-graph moment check with a Gamma(alpha,
# beta) diagonal and a small eigenvalue floor eps. The current manuscript uses
# the exponential diagonal and no floor, and does not cite this check.
#
# Empty-graph moment verification.
# On the empty graph the tilted prior reduces to iid Gamma(alpha + delta, beta)
# on the diagonal entries (truncated by the small numerical-stability
# indicator). The script verifies the closed-form moments by direct simulation
# and reports the empirical-vs-closed-form table for E[theta_ll], Var[theta_ll],
# E[log theta_ll], E[trace], E[||Theta||_F^2], E[log|Theta|].

suppressPackageStartupMessages({
  library(parallel)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
OUT_RDS      <- file.path(PROJECT_ROOT, "R", "results", "empty_graph_moments.rds")

# Direct sampler from the tilted empty-graph prior:
# theta_ll iid Gamma(alpha + delta, rate = beta), truncated to (eps, infty).
# The tilt |Theta|^delta with Theta = diag(theta_ll) factors entry-by-entry to
# theta_ll^delta, which combined with the Gamma(alpha, beta) base measure gives
# Gamma(alpha + delta, beta) marginally.
sample_empty <- function(q, alpha, beta, delta, eps, n_draws, seed) {
  set.seed(seed)
  M <- matrix(stats::rgamma(n_draws * q, shape = alpha + delta, rate = beta),
              nrow = n_draws)
  # Apply numerical-stability truncation eps. Rejection on rows where any
  # entry <= eps.
  if (eps > 0) {
    keep <- apply(M > eps, 1, all)
    M <- M[keep, , drop = FALSE]
  }
  M
}

closed_form_moments <- function(q, alpha, beta, delta) {
  shape <- alpha + delta
  list(
    E_theta_ll      = shape / beta,
    Var_theta_ll    = shape / beta^2,
    E_log_theta_ll  = digamma(shape) - log(beta),
    E_trace         = q * shape / beta,
    Var_trace       = q * shape / beta^2,
    E_frob_sq       = q * shape * (shape + 1) / beta^2,
    E_log_det       = q * (digamma(shape) - log(beta))
  )
}

empirical_moments <- function(M) {
  list(
    E_theta_ll      = mean(M),
    Var_theta_ll    = stats::var(as.vector(M)),
    E_log_theta_ll  = mean(log(M)),
    E_trace         = mean(rowSums(M)),
    Var_trace       = stats::var(rowSums(M)),
    E_frob_sq       = mean(rowSums(M^2)),
    E_log_det       = mean(rowSums(log(M)))
  )
}

# Sweep over (q, delta) with alpha = beta = 3 throughout
ALPHA   <- 3.0
BETA    <- 3.0
EPS     <- 1e-3
N_DRAWS <- 5e5

q_grid     <- c(5L, 20L, 50L)
delta_grid <- c(0, 0.5, 1, 2)

jobs <- list()
k <- 0L
for (q in q_grid) {
  for (d in delta_grid) {
    k <- k + 1L
    jobs[[length(jobs) + 1L]] <- list(
      q = q, alpha = ALPHA, beta = BETA, delta = d, eps = EPS,
      n_draws = N_DRAWS, seed = 42L + 1000L * k
    )
  }
}

cat(sprintf("Running %d empty-graph cells on 12 cores...\n", length(jobs)))
t0 <- Sys.time()
results <- parallel::mclapply(jobs, function(j) {
  M <- sample_empty(j$q, j$alpha, j$beta, j$delta, j$eps, j$n_draws, j$seed)
  cf <- closed_form_moments(j$q, j$alpha, j$beta, j$delta)
  em <- empirical_moments(M)
  list(q = j$q, alpha = j$alpha, beta = j$beta, delta = j$delta, eps = j$eps,
       n_kept = nrow(M),
       closed = cf, empirical = em)
}, mc.cores = 12L)
t1 <- Sys.time()
cat(sprintf("Done in %.2f s.\n", as.numeric(t1 - t0, units = "secs")))

# Build a tidy long-format table
rows <- list()
metrics <- c("E_theta_ll", "Var_theta_ll", "E_log_theta_ll",
             "E_trace", "Var_trace", "E_frob_sq", "E_log_det")
for (r in results) {
  for (m in metrics) {
    rows[[length(rows) + 1L]] <- data.frame(
      q       = r$q,
      delta   = r$delta,
      metric  = m,
      closed_form = r$closed[[m]],
      empirical   = r$empirical[[m]],
      rel_err     = (r$empirical[[m]] - r$closed[[m]]) / r$closed[[m]],
      stringsAsFactors = FALSE
    )
  }
}
tab <- do.call(rbind, rows)

cat("\n=== Empty-graph moment verification ===\n")
print(tab[order(tab$metric, tab$q, tab$delta), ], row.names = FALSE, digits = 4)

saveRDS(list(jobs = jobs, results = results, summary = tab,
             alpha = ALPHA, beta = BETA, eps = EPS, n_draws = N_DRAWS),
        OUT_RDS)
cat(sprintf("\nSaved empty-graph diagnostics to %s\n", OUT_RDS))
