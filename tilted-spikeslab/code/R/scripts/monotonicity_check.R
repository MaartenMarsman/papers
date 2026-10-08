# R/scripts/monotonicity_check.R
# Purpose: exhaustive small-graph verification of the edge-monotonicity theorem
#          c(0) <= Z(Gamma + e; delta) / Z(Gamma; delta) < 1,
#          with equality at the lower end exactly when the endpoints of e are
#          disconnected in Gamma. All graphs at q = 3, 4 (full hyperparameter
#          grid) and q = 5 (two settings). alpha = 1, standardized scale
#          sigma = 1, beta = eta.
# Run:     Rscript R/scripts/monotonicity_check.R
# Output:  R/results/monotonicity_check.rds (+ console summary)

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
set.seed(20260718)

suppressMessages({
  library(Rcpp)
  library(parallel)
})

sourceCpp(file.path(PROJECT_ROOT, "R", "src", "monotonicity_check.cpp"))

# --- parameters ---
settings <- data.frame(
  eta   = c(0.5, 1.0, 1.0, 2.0, 1.0),
  delta = c(1.6, 0.8, 1.6, 1.6, 0.0)
)
qs       <- c(3, 4, 5)
n_draws  <- c(`3` = 2e6, `4` = 1e6, `5` = 5e5)
q5_keep  <- c(2, 3)          # settings rows to run at q = 5 (runtime)
n_c0     <- 3e7              # draws for the 2x2 constant c(0)
n_cores  <- 12L
seed_mc  <- 42L              # shared across graphs -> common random numbers

# --- helpers ---
edge_list <- function(q) t(combn(q, 2))

graph_from_mask <- function(mask, edges, q) {
  G <- matrix(0L, q, q)
  on <- which(bitwAnd(mask, bitwShiftL(1L, seq_len(nrow(edges)) - 1L)) != 0L)
  for (e in on) {
    G[edges[e, 1], edges[e, 2]] <- 1L
    G[edges[e, 2], edges[e, 1]] <- 1L
  }
  G
}

components <- function(G) {
  q <- nrow(G); comp <- rep(0L, q); cur <- 0L
  for (s in seq_len(q)) {
    if (comp[s] > 0L) next
    cur <- cur + 1L; stack <- s
    while (length(stack)) {
      v <- stack[[1]]; stack <- stack[-1]
      if (comp[v] > 0L) next
      comp[v] <- cur
      stack <- c(stack, which(G[v, ] == 1L & comp == 0L))
    }
  }
  comp
}

# c(0) = Z(2x2, one edge) / Z(2x2, empty); the empty value is analytic:
# E[k^delta] under Exp(beta) is Gamma(1 + delta) / beta^delta.
c0_constant <- function(eta, delta) {
  G2 <- matrix(c(0L, 1L, 1L, 0L), 2, 2)
  z2 <- z_graph_aligned_cpp(G2, beta = eta, sigma = 1, delta = delta,
                            n_draws = n_c0, seed = seed_mc)
  z_empty <- (gamma(1 + delta) / eta^delta)^2
  list(c0 = z2$z_hat / z_empty, se = z2$z_se / z_empty)
}

# --- computation ---
results <- list()

for (q in qs) {
  edges <- edge_list(q)
  E <- nrow(edges)
  masks <- 0:(2^E - 1)
  N <- n_draws[[as.character(q)]]
  rows <- if (q == 5) q5_keep else seq_len(nrow(settings))

  for (r in rows) {
    eta <- settings$eta[r]; delta <- settings$delta[r]
    t0 <- Sys.time()

    z <- mclapply(masks, function(m) {
      G <- graph_from_mask(m, edges, q)
      z_graph_aligned_cpp(G, beta = eta, sigma = 1, delta = delta,
                          n_draws = N, seed = seed_mc)
    }, mc.cores = n_cores)
    z_hat <- vapply(z, `[[`, 0, "z_hat")
    z_se  <- vapply(z, `[[`, 0, "z_se")

    c0 <- c0_constant(eta, delta)

    # analytic empty-graph value as an MC calibration check
    z_empty_true <- (gamma(1 + delta) / eta^delta)^q
    empty_relerr <- z_hat[1] / z_empty_true - 1

    # every (graph, absent edge) pair
    checks <- do.call(rbind, lapply(masks, function(m) {
      G <- graph_from_mask(m, edges, q)
      comp <- components(G)
      absent <- which(bitwAnd(m, bitwShiftL(1L, seq_len(E) - 1L)) == 0L)
      if (!length(absent)) return(NULL)
      do.call(rbind, lapply(absent, function(e) {
        m_plus <- bitwOr(m, bitwShiftL(1L, e - 1L))
        ratio <- z_hat[m_plus + 1] / z_hat[m + 1]
        # conservative SE (ignores the positive coupling correlation)
        se <- ratio * sqrt((z_se[m_plus + 1] / z_hat[m_plus + 1])^2 +
                           (z_se[m + 1] / z_hat[m + 1])^2)
        data.frame(mask = m, edge = e,
                   disconnected = comp[edges[e, 1]] != comp[edges[e, 2]],
                   ratio = ratio, se = se)
      }))
    }))

    tol <- 4 * (checks$se + c0$se)
    checks$below_lower <- checks$ratio < c0$c0 - tol
    checks$above_one <- checks$ratio >= 1 + tol
    checks$disc_dev <- ifelse(checks$disconnected,
                              abs(checks$ratio - c0$c0), NA)

    res <- list(q = q, eta = eta, delta = delta, n_draws = N,
                c0 = c0$c0, c0_se = c0$se,
                empty_relerr = empty_relerr,
                n_pairs = nrow(checks),
                n_below_lower = sum(checks$below_lower),
                n_at_or_above_one = sum(checks$above_one),
                max_ratio = max(checks$ratio),
                min_ratio = min(checks$ratio),
                max_disc_dev = suppressWarnings(max(checks$disc_dev, na.rm = TRUE)),
                max_disc_dev_tol = suppressWarnings(
                  max(tol[checks$disconnected], na.rm = TRUE)),
                checks = checks,
                minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
    results[[sprintf("q%d_eta%g_delta%g", q, eta, delta)]] <- res

    cat(sprintf(
      "q=%d eta=%.2g delta=%.2g | c0=%.4f(%.1e) emptyerr=%+.2e | pairs=%d below=%d >=1: %d | min=%.4f max=%.6f | discdev=%.2e (tol %.2e) | %.1f min\n",
      q, eta, delta, c0$c0, c0$se, empty_relerr,
      res$n_pairs, res$n_below_lower, res$n_at_or_above_one,
      res$min_ratio, res$max_ratio, res$max_disc_dev, res$max_disc_dev_tol,
      res$minutes))
  }
}

# --- save ---
dir.create("R/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(results, file.path(PROJECT_ROOT, "R", "results", "monotonicity_check.rds"))

ok <- all(vapply(results, function(r)
  r$n_below_lower == 0 && r$n_at_or_above_one == 0 &&
    r$max_disc_dev <= r$max_disc_dev_tol, TRUE))
cat(sprintf("\nOVERALL: %s\n", if (ok) "PASS" else "FAIL"))
cat("Done. Results saved to R/results/monotonicity_check.rds\n")
