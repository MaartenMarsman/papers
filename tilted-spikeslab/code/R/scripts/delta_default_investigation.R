# R/scripts/delta_default_investigation.R
#
# EXPLORATORY, not used in the manuscript: uses the HMC sampler of the bgms
# package as the prior sampler. Superseded by delta_scaling_regime_sweep.R,
# which answers the same question with the row-block Gibbs sampler. Its output
# is not archived.
#
# Informed default for the tilt exponent delta.
#
# Question: the earlier rule delta = 0.5 log q was read off a sweep that held
# the expected degree roughly constant (edge-inclusion probability decreasing
# with q). That is one design choice. How does the delta needed to hold a fixed
# boundary margin depend on the graph regime, and what default should bgms use?
#
# Two things are computed, both in the standardized frame (exponential diagonal,
# sigma = 1, eta = sigma*beta = 1, so beta = 1; bgms takes eta_bgms = 2*beta):
#
# 1. delta*(q, regime): the tilt that holds the median smallest eigenvalue of K
#    at target levels {0.25, 0.5}, across four graph-density regimes:
#      degree2  : expected degree ~ 2  (p_inc = 2/(q-1)); the earlier schedule
#      degree4  : expected degree ~ 4  (p_inc = 4/(q-1))
#      density25: constant density 0.25 (degree grows with q)
#      complete : p_inc = 1
#    Sampled from the tilted PD-conditional prior with bgms (HMC), which is PD by
#    construction, over a delta-grid; delta* is found by monotone interpolation.
#
# 2. Slack law check: the diagonal slack psi_l = 1/(K^{-1})_ll should be
#    Gamma(delta+1, beta) exactly, regardless of q or graph (Appendix B). We
#    confirm this from the HMC draws at a few (q, regime, delta) cells, including
#    the complete graph at q = 50 where importance sampling from the base measure
#    cannot reach.
#
# Run:    Rscript R/scripts/delta_default_investigation.R
# Output: R/results/delta_default_investigation.rds (+ console summary)

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
set.seed(20260719L)
suppressPackageStartupMessages({
  library(bgms)
  library(parallel)
})

# --- frame ---
SIGMA    <- 1
BETA     <- 1
ETA_BGMS <- 2 * BETA          # bgms exponential_prior(eta) sets diagonal rate eta/2
QS       <- c(5L, 10L, 15L, 20L, 30L, 50L, 75L, 100L)
DELTAS   <- c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12)
TARGETS  <- c(0.25, 0.5)
NGRAPH   <- 4L
NS       <- 400L
NW       <- 500L
NCORES   <- 12L

regime_dens <- function(regime, q) {
  switch(regime,
    degree2   = min(1, 2 / (q - 1)),
    degree4   = min(1, 4 / (q - 1)),
    density25 = 0.25,
    complete  = 1,
    stop("unknown regime"))
}
REGIMES <- c("degree2", "degree4", "density25", "complete")

# --- helpers ---
mk_graph <- function(q, dens, seed) {
  set.seed(seed)
  A <- matrix(0L, q, q)
  if (dens >= 1) { A[] <- 1L; diag(A) <- 1L; return(A) }
  A[upper.tri(A)] <- rbinom(choose(q, 2), 1L, dens)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 1L
  A
}

# Reconstruct K from a bgms prior-sample row and return (lambda_min, slacks).
reconstruct <- function(r, q) {
  di <- do.call(rbind, lapply(strsplit(sub("^K_", "", r$diag_names),    "_"), as.integer))
  oi <- do.call(rbind, lapply(strsplit(sub("^K_", "", r$offdiag_names), "_"), as.integer))
  S <- nrow(r$K_diag)
  lmin  <- numeric(S)
  slack <- vector("list", S)
  for (s in seq_len(S)) {
    K <- matrix(0, q, q)
    K[cbind(di[, 1], di[, 2])] <- r$K_diag[s, ]
    K[cbind(oi[, 1], oi[, 2])] <- r$K_offdiag[s, ]
    K[cbind(oi[, 2], oi[, 1])] <- r$K_offdiag[s, ]
    ev <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
    lmin[s] <- min(ev)
    Ki <- chol2inv(chol(K))
    slack[[s]] <- 1 / diag(Ki)
  }
  list(lmin = lmin, slack = unlist(slack))
}

# One (regime, q) cell: sample at every delta over NGRAPH graphs, return the
# median lambda_min per delta and pooled slacks at each delta.
run_cell <- function(regime, q) {
  dens <- regime_dens(regime, q)
  out_lmin  <- matrix(NA_real_, length(DELTAS), 3,
                      dimnames = list(NULL, c("lmin_med", "lmin_q05", "lmin_q95")))
  slack_by_delta <- vector("list", length(DELTAS))
  for (di in seq_along(DELTAS)) {
    delta <- DELTAS[di]
    lm <- numeric(0); sl <- numeric(0)
    for (g in seq_len(NGRAPH)) {
      A <- mk_graph(q, dens, 9000L + q * 131L + g * 7L)
      r <- sample_ggm_prior(p = q, n_samples = NS, n_warmup = NW,
             interaction_prior = normal_prior(scale = SIGMA),
             precision_scale_prior = exponential_prior(eta = ETA_BGMS),
             delta = delta, spec = "conditional", edge_indicators = A,
             step_size = 0.1, seed = 20260719L + di * 1000L + g, verbose = FALSE)
      rc <- reconstruct(r, q)
      lm <- c(lm, rc$lmin)
      if (delta %in% c(1, 2, 4)) sl <- c(sl, rc$slack)  # store slacks at a few deltas
    }
    out_lmin[di, ] <- c(median(lm), quantile(lm, 0.05), quantile(lm, 0.95))
    slack_by_delta[[di]] <- if (length(sl)) sl else NULL
  }
  # invert median lambda_min(delta) for each target (monotone increasing in delta)
  med <- out_lmin[, "lmin_med"]
  dstar <- sapply(TARGETS, function(tt) {
    if (all(med < tt)) return(NA_real_)          # even delta=max too small
    if (med[1] >= tt) return(0)                   # delta=0 already enough
    approx(x = med, y = DELTAS, xout = tt, ties = "ordered")$y
  })
  names(dstar) <- paste0("target", TARGETS)
  list(regime = regime, q = q, dens = dens,
       deltas = DELTAS, lmin = out_lmin, dstar = dstar,
       slack_by_delta = setNames(slack_by_delta, DELTAS))
}

# --- run all cells in parallel ---
cells <- expand.grid(regime = REGIMES, q = QS, stringsAsFactors = FALSE)
t0 <- Sys.time()
results <- mcmapply(function(regime, q) run_cell(regime, q),
                    cells$regime, cells$q,
                    SIMPLIFY = FALSE, mc.cores = NCORES, mc.preschedule = FALSE)
cat(sprintf("Sweep done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# --- assemble delta* table ---
dstar_df <- do.call(rbind, lapply(results, function(r)
  data.frame(regime = r$regime, q = r$q, dens = round(r$dens, 3),
             exp_degree = round(r$dens * (r$q - 1), 2),
             dstar_med025 = r$dstar["target0.25"],
             dstar_med050 = r$dstar["target0.5"],
             row.names = NULL)))
dstar_df <- dstar_df[order(dstar_df$regime, dstar_df$q), ]

# --- slack law summary: compare pooled slacks to Gamma(delta+1, beta) ---
slack_df <- do.call(rbind, lapply(results, function(r) {
  do.call(rbind, lapply(c("1", "2", "4"), function(dk) {
    sl <- r$slack_by_delta[[dk]]
    if (is.null(sl)) return(NULL)
    delta <- as.numeric(dk)
    data.frame(regime = r$regime, q = r$q, delta = delta,
      slack_med = median(sl), gamma_med = qgamma(0.5, delta + 1, rate = BETA),
      slack_q10 = quantile(sl, 0.1), gamma_q10 = qgamma(0.1, delta + 1, rate = BETA),
      slack_q90 = quantile(sl, 0.9), gamma_q90 = qgamma(0.9, delta + 1, rate = BETA),
      row.names = NULL)
  }))
}))

saveRDS(list(results = results, dstar = dstar_df, slack = slack_df,
             frame = list(sigma = SIGMA, beta = BETA, deltas = DELTAS,
                          targets = TARGETS, ngraph = NGRAPH, ns = NS, nw = NW)),
        file.path(PROJECT_ROOT, "R", "results", "delta_default_investigation.rds"))

cat("\n=== delta* to hold median lambda_min (by regime, q) ===\n")
print(dstar_df, digits = 3, row.names = FALSE)
cat("\n=== delta*/log(q) for target median 0.5 ===\n")
tab <- dstar_df
tab$ratio <- tab$dstar_med050 / log(tab$q)
print(tab[, c("regime", "q", "exp_degree", "dstar_med050", "ratio")], digits = 3, row.names = FALSE)
cat("\n=== slack law: pooled slack vs Gamma(delta+1, beta) ===\n")
print(slack_df, digits = 3, row.names = FALSE)
cat("\nDone. Saved R/results/delta_default_investigation.rds\n")
