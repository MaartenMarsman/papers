# R/scripts/tilt_LO_check.R
#
# Numerical verification of the tilt-LO formula derived in notes/tilt-mapping.md.
# Compares log Z_LO^tilt against the brute-force MC values from Phase 1.4
# (R/results/Z_prior_brute_q234.rds) on the discrete variant.
#
# Tested at (alpha, beta) = (3, 3), sigma = 1, delta in {0, 0.5, 1}, epsilon = 0
# (the floor is empirically negligible at LO for our parameter range; see
# notes/floor-and-NLO.md).

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
BRUTE_RDS    <- file.path(PROJECT_ROOT, "R", "results", "Z_prior_brute_q234.rds")
OUT_RDS      <- file.path(PROJECT_ROOT, "R", "results", "tilt_LO_check.rds")

# ----- tilt-LO formula -----
#
# log Z_LO^tilt(G; alpha, beta, sigma, delta) =
#   sum_v log Gamma((nu_v + 2(alpha+delta))/2) - q log Gamma(alpha)
#   - |E| log beta + (|E|/2) log pi - (|E|/2) log(2 pi sigma^2)
#   + 0.5 sum_e log(2 beta / H_e^tilt)
# with H_e^tilt = 2 beta + (M_k^tilt) / (2 sigma^2 beta)
#                 - 4 beta (alpha - 1) / M_v^tilt
# and M_v^tilt = nu_v + 2(alpha + delta) - 1.

log_Z_LO_tilt <- function(G, alpha, beta, sigma, delta) {
  q <- nrow(G)
  edge_mat <- which(G == 1L & upper.tri(G), arr.ind = TRUE)
  n_edges <- nrow(edge_mat)

  nu <- integer(q)
  if (q > 1L) {
    for (l in 1:(q - 1L)) {
      for (m in (l + 1L):q) {
        if (G[l, m] == 1L) nu[l] <- nu[l] + 1L
      }
    }
  }

  M_tilt <- nu + 2 * (alpha + delta) - 1
  sigma2 <- sigma * sigma

  log_C0 <- sum(lgamma((nu + 2 * (alpha + delta)) / 2)) -
            q * lgamma(alpha) -
            n_edges * log(beta) -
            q * delta * log(beta) +
            (n_edges / 2) * log(pi) -
            (n_edges / 2) * log(2 * pi * sigma2)

  log_slab <- 0
  if (n_edges > 0L) {
    H_e <- numeric(n_edges)
    for (e in seq_len(n_edges)) {
      k <- edge_mat[e, 1L]; v <- edge_mat[e, 2L]
      H_e[e] <- 2 * beta + M_tilt[k] / (2 * sigma2 * beta) -
                4 * beta * (alpha - 1) / M_tilt[v]
    }
    if (any(H_e <= 0)) {
      return(list(log_z = NaN, H_e = H_e, M_tilt = M_tilt, note = "non-positive H_e"))
    }
    log_slab <- 0.5 * sum(log(2 * beta / H_e))
  }

  list(log_z = log_C0 + log_slab,
       log_C0 = log_C0,
       log_slab = log_slab,
       M_tilt = M_tilt,
       H_e = if (n_edges > 0L) H_e else NULL)
}

# ----- graph builder (matches Z_prior_brute_force.R) -----
build_graphs <- function() {
  g2_empty <- matrix(0L, 2, 2)
  g2_full  <- matrix(0L, 2, 2); g2_full[1,2] <- g2_full[2,1] <- 1L

  g3_empty <- matrix(0L, 3, 3)
  g3_path  <- matrix(0L, 3, 3)
  g3_path[1,2] <- g3_path[2,1] <- 1L
  g3_path[2,3] <- g3_path[3,2] <- 1L
  g3_tri   <- matrix(1L, 3, 3); diag(g3_tri) <- 0L

  g4_empty <- matrix(0L, 4, 4)
  g4_path  <- matrix(0L, 4, 4)
  g4_path[1,2] <- g4_path[2,1] <- 1L
  g4_path[2,3] <- g4_path[3,2] <- 1L
  g4_path[3,4] <- g4_path[4,3] <- 1L
  g4_cyc   <- g4_path; g4_cyc[1,4] <- g4_cyc[4,1] <- 1L
  g4_full  <- matrix(1L, 4, 4); diag(g4_full) <- 0L

  list(
    q2_empty = g2_empty, q2_full = g2_full,
    q3_empty = g3_empty, q3_path = g3_path, q3_tri = g3_tri,
    q4_empty = g4_empty, q4_path = g4_path, q4_cyc = g4_cyc, q4_full = g4_full
  )
}

# ----- closed-form empty-graph (tilt-LO is exact here) -----
log_z_closed_empty <- function(q, alpha, beta, delta) {
  q * (lgamma(alpha + delta) - lgamma(alpha) - delta * log(beta))
}

# ----- run the comparison at epsilon = 0 (LO is for the un-floored prior) -----
graphs <- build_graphs()
brute  <- readRDS(BRUTE_RDS)
brute_disc <- brute[brute$variant == "discrete" & brute$epsilon == 0, ]

rows <- list()
for (gname in names(graphs)) {
  G <- graphs[[gname]]
  for (delta in c(0, 0.5, 1.0)) {
    lo <- log_Z_LO_tilt(G, alpha = 3, beta = 3, sigma = 1, delta = delta)
    brow <- brute_disc[brute_disc$graph == gname & brute_disc$delta == delta, ]
    if (nrow(brow) != 1L) {
      stop(sprintf("expected 1 brute row for %s delta=%g, got %d", gname, delta, nrow(brow)))
    }
    rows[[length(rows) + 1L]] <- data.frame(
      graph    = gname,
      q        = nrow(G),
      delta    = delta,
      log_z_LO = lo$log_z,
      log_z_brute = brow$log_z,
      miss_LO  = lo$log_z - brow$log_z,
      z_brute  = brow$z_hat,
      z_se_brute = brow$z_se,
      stringsAsFactors = FALSE
    )
  }
}
tab <- do.call(rbind, rows)

# Cross-check on empty graphs: closed-form should equal LO exactly
empty_rows <- tab[grepl("_empty$", tab$graph), ]
empty_rows$log_z_closed <- mapply(function(q, d) log_z_closed_empty(q, 3, 3, d),
                                   empty_rows$q, empty_rows$delta)
empty_rows$LO_vs_closed <- empty_rows$log_z_LO - empty_rows$log_z_closed

cat("=== tilt-LO vs brute-force (epsilon = 0, discrete variant) ===\n\n")
print(tab[order(tab$q, tab$graph, tab$delta), ], row.names = FALSE, digits = 5)

cat("\n=== empty-graph: tilt-LO vs closed-form (should match exactly) ===\n\n")
print(empty_rows[order(empty_rows$q, empty_rows$delta),
                 c("graph", "delta", "log_z_LO", "log_z_closed", "LO_vs_closed")],
      row.names = FALSE, digits = 6)

cat("\n=== LO miss summary (|miss| < 0.5 nats expected for q <= 4 sparse graphs) ===\n\n")
tab$abs_miss <- abs(tab$miss_LO)
agg <- aggregate(abs_miss ~ q + delta, data = tab, FUN = mean)
print(agg, row.names = FALSE, digits = 4)

saveRDS(list(tab = tab, empty = empty_rows), OUT_RDS)
cat(sprintf("\nSaved tilt-LO check to %s\n", OUT_RDS))
