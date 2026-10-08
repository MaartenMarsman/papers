# Regime-explicit delta-scaling sweep (redesign of delta_scaling_sweep.R).
# See notes/delta-scaling-regime-dependence.md.
#
# Design: factorial over
#   q            in {10, 20, 50}
#   regime axis  fixed expected degree dbar in {2, 4, 8}  (p = dbar/(q-1))
#                fixed density        rho  in {0.05, 0.1, 0.2}  (p = rho)
#   eta          in {0.5, 1, 2}          (standardized: sigma = 1, beta = eta)
#   delta grid   {0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24}
#   graph reps   2 (independent Erdos-Renyi draws per cell)
#
# Sampler: exact prior row-block Gibbs (the manuscript's within-model
# conditionals at n = 0): per node l, given the rest,
#   psi ~ Gamma(delta + 1, beta),   k_F ~ N(0, (sigma^-2 I + 2 beta M_F)^-1),
#   k_ll = psi + k_F' M_F k_F,      M_F = ((K_{-l})^{-1})_{F,F},
# with F the neighbors of l. Positive definiteness is automatic.
#
# Validation (runs first, aborts on failure): q = 5 fixed graph, Gibbs vs
# importance resampling from the untilted base prior with weights |K|^delta 1{PD}.
#
# Metrics per (cell, delta): median/IQR of lambda_min; median diagonal;
# realized slab variance ratio mean(k_edge^2)/sigma^2; SD of the partial
# correlation over edges and draws; edge count m.
#
# Output: R/results/delta_scaling_regime_sweep.rds (+ .log)

suppressPackageStartupMessages(library(parallel))

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
OUT_RDS <- file.path(PROJECT_ROOT, "R", "results", "delta_scaling_regime_sweep.rds")

SIGMA <- 1

## ---------- prior row-block Gibbs ----------

gibbs_prior_chain <- function(adj, delta, beta, n_burn, n_keep, seed) {
  set.seed(seed)
  q <- nrow(adj)
  K <- diag((delta + 1) / beta, q)
  nbrs <- lapply(seq_len(q), function(l) which(adj[l, ] == 1L))
  draws_lmin <- numeric(n_keep)
  draws_diag_med <- numeric(n_keep)
  edge_idx <- which(upper.tri(adj) & adj == 1L, arr.ind = TRUE)
  m <- nrow(edge_idx)
  sum_e2 <- 0; sum_rho <- 0; sum_rho2 <- 0; n_rho <- 0
  n_sweeps <- n_burn + n_keep
  for (it in seq_len(n_sweeps)) {
    for (l in seq_len(q)) {
      Fl <- nbrs[[l]]
      psi <- stats::rgamma(1L, shape = delta + 1, rate = beta)
      if (length(Fl) == 0L) {
        K[l, ] <- 0; K[, l] <- 0; K[l, l] <- psi
      } else {
        Km <- K[-l, -l, drop = FALSE]
        Fm <- match(Fl, seq_len(q)[-l])
        rhs <- matrix(0, q - 1L, length(Fm))
        rhs[cbind(Fm, seq_along(Fm))] <- 1
        MF <- solve(Km, rhs)[Fm, , drop = FALSE]
        MF <- (MF + t(MF)) / 2
        Q <- diag(1 / SIGMA^2, length(Fm)) + 2 * beta * MF
        R <- chol(Q)
        kF <- backsolve(R, stats::rnorm(length(Fm)))
        row_l <- numeric(q)
        row_l[Fl] <- kF
        K[l, ] <- row_l; K[, l] <- row_l
        K[l, l] <- psi + drop(crossprod(kF, MF %*% kF))
      }
    }
    if (it > n_burn) {
      j <- it - n_burn
      ev <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
      draws_lmin[j] <- min(ev)
      dg <- diag(K)
      draws_diag_med[j] <- stats::median(dg)
      if (m > 0L) {
        e2 <- K[edge_idx]^2
        rho <- -K[edge_idx] / sqrt(dg[edge_idx[, 1]] * dg[edge_idx[, 2]])
        sum_e2 <- sum_e2 + sum(e2)
        sum_rho <- sum_rho + sum(rho)
        sum_rho2 <- sum_rho2 + sum(rho^2)
        n_rho <- n_rho + m
      }
    }
  }
  var_ratio <- if (m > 0L) (sum_e2 / n_rho) / SIGMA^2 else NA_real_
  rho_sd <- if (m > 0L) sqrt(sum_rho2 / n_rho - (sum_rho / n_rho)^2) else NA_real_
  list(lmin_med = stats::median(draws_lmin),
       lmin_q25 = stats::quantile(draws_lmin, 0.25, names = FALSE),
       lmin_q75 = stats::quantile(draws_lmin, 0.75, names = FALSE),
       diag_med = stats::median(draws_diag_med),
       var_ratio = var_ratio, rho_sd = rho_sd, m = m,
       lmin_first_half = stats::median(draws_lmin[seq_len(n_keep %/% 2)]),
       lmin_second_half = stats::median(draws_lmin[(n_keep %/% 2 + 1L):n_keep]))
}

## ---------- validation: Gibbs vs importance resampling at q = 5 ----------

weighted_median <- function(x, w) {
  o <- order(x)
  cw <- cumsum(w[o]) / sum(w)
  x[o][which(cw >= 0.5)[1]]
}

validate <- function() {
  set.seed(11)
  q <- 5L
  adj <- matrix(0L, q, q)
  edges <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(2, 5))
  adj[edges] <- 1L; adj[edges[, c(2, 1)]] <- 1L
  ok <- TRUE
  for (delta in c(0.5, 2)) {
    for (eta in c(1)) {
      beta <- eta
      n_is <- 4e5
      diag_draws <- matrix(stats::rgamma(n_is * q, 1, beta), n_is, q)
      slab_draws <- matrix(stats::rnorm(n_is * nrow(edges), 0, SIGMA), n_is, nrow(edges))
      lmin <- numeric(n_is); logw <- numeric(n_is)
      for (i in seq_len(n_is)) {
        K <- diag(diag_draws[i, ])
        K[edges] <- slab_draws[i, ]; K[edges[, c(2, 1)]] <- slab_draws[i, ]
        ev <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
        lmin[i] <- min(ev)
        logw[i] <- if (min(ev) > 0) delta * sum(log(ev)) else -Inf
      }
      w <- exp(logw - max(logw[is.finite(logw)]))
      keep <- is.finite(logw)
      ref_med <- weighted_median(lmin[keep], w[keep])
      ess <- sum(w[keep])^2 / sum(w[keep]^2)
      gb <- gibbs_prior_chain(adj, delta, beta, n_burn = 500, n_keep = 4000, seed = 22)
      rel <- abs(gb$lmin_med / ref_med - 1)
      cat(sprintf("validate delta=%.1f eta=%.1f: IS median lmin = %.4f (ESS %.0f), Gibbs = %.4f, rel diff = %.3f\n",
                  delta, eta, ref_med, ess, gb$lmin_med, rel))
      if (rel > 0.05) ok <- FALSE
    }
  }
  ok
}

cat("== validation ==\n")
if (!validate()) stop("Gibbs sampler failed validation against importance resampling", call. = FALSE)
cat("validation PASSED\n\n")

## ---------- factorial sweep ----------

qs <- c(10L, 20L, 50L)
regimes <- rbind(
  data.frame(axis = "degree", level = c(2, 4, 8)),
  data.frame(axis = "density", level = c(0.05, 0.10, 0.20))
)
etas <- c(0.5, 1, 2)
deltas <- c(0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24)
n_reps <- 2L

cells <- list()
sc <- 0L
for (q in qs) for (r in seq_len(nrow(regimes))) for (eta in etas)
  for (delta in deltas) for (rep in seq_len(n_reps)) {
    sc <- sc + 1L
    p_edge <- if (regimes$axis[r] == "degree") min(regimes$level[r] / (q - 1), 1) else regimes$level[r]
    cells[[sc]] <- list(q = q, axis = regimes$axis[r], level = regimes$level[r],
                        p_edge = p_edge, eta = eta, delta = delta, rep = rep,
                        seed = 40000L + sc)
  }
cat(sprintf("running %d cells on 12 cores\n", length(cells)))

run_cell <- function(cell) {
  set.seed(cell$seed)
  q <- cell$q
  adj <- matrix(0L, q, q)
  ut <- upper.tri(adj)
  adj[ut] <- stats::rbinom(sum(ut), 1L, cell$p_edge)
  adj <- adj + t(adj)
  res <- gibbs_prior_chain(adj, cell$delta, cell$eta, n_burn = 300, n_keep = 500,
                           seed = cell$seed + 1L)
  data.frame(q = q, axis = cell$axis, level = cell$level, p_edge = cell$p_edge,
             eta = cell$eta, delta = cell$delta, rep = cell$rep,
             m = res$m, dbar_realized = 2 * res$m / q,
             lmin_med = res$lmin_med, lmin_q25 = res$lmin_q25, lmin_q75 = res$lmin_q75,
             diag_med = res$diag_med, var_ratio = res$var_ratio, rho_sd = res$rho_sd,
             lmin_first_half = res$lmin_first_half, lmin_second_half = res$lmin_second_half)
}

t0 <- Sys.time()
out <- mclapply(cells, run_cell, mc.cores = 12L)
bad <- vapply(out, function(x) inherits(x, "try-error") || is.null(x), logical(1))
if (any(bad)) cat(sprintf("WARNING: %d failed cells\n", sum(bad)))
res <- do.call(rbind, out[!bad])
saveRDS(res, OUT_RDS)
cat(sprintf("done in %.1f min; %d rows -> %s\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(res), OUT_RDS))
