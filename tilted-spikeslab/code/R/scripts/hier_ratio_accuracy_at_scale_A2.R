# Section 5.4, component A2: estimator precision and coupling AT SCALE.
#
# q in {10, 20, 50, 100}, Erdos-Renyi graphs of expected degree 4 (fixed
# seed), one representative absent edge per graph (fixed), eta = 1,
# delta = 0.5 log q, N in {1e2, 1e3, 1e4}, coupled and independent arms,
# R = 100 replicate estimates per cell.
#
# Reported: SD of log J-hat (the operative precision; bias was undetectable
# at q = 4 and is re-anchored here at q in {10, 20} against a 20 x N = 1e6
# coupled reference), coupling correlation, and the SD ratio
# independent/coupled.
#
# Output: R/results/hier_ratio_accuracy_at_scale_A2.rds (+ log)

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})
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

qs <- c(10L, 20L, 50L, 100L)
Ns <- c(100L, 1000L, 10000L)
R_reps <- 100L

cells <- list()
for (q in qs) for (N in Ns) for (arm in c("coupled", "independent"))
  cells[[length(cells) + 1L]] <- list(q = q, N = N, arm = arm)

pick_edge <- function(G) {
  q <- nrow(G)
  absent <- which(upper.tri(G) & G == 0, arr.ind = TRUE)
  as.integer(absent[nrow(absent) %/% 2, ])
}

run_cell <- function(cell) {
  q <- cell$q
  G <- er_graph(q, 4500L + q)
  e <- pick_edge(G)
  delta <- 0.5 * log(q)
  set.seed(4600L + q + cell$N %% 97L + (cell$arm == "coupled"))
  reps <- t(vapply(seq_len(R_reps), function(r) {
    est <- akm_ratio_cpp(G, e[1], e[2], delta, 1, 1, cell$N,
                         coupled = (cell$arm == "coupled"))
    c(est$log_J, est$mean_minus, est$mean_plus)
  }, numeric(3)))
  # log-weight diagnostics: the collapse mechanism (one diagnostic call);
  # draws are LOG weights, -Inf when the completion overflows
  dd <- akm_ratio_cpp(G, e[1], e[2], delta, 1, 1, min(cell$N, 1000L),
                      coupled = TRUE, return_draws = TRUE)
  lw <- dd$draws_minus
  data.frame(q = q, N = cell$N, arm = cell$arm,
             mean_logJ = mean(reps[, 1]), sd_logJ = sd(reps[, 1]),
             corr_means = suppressWarnings(cor(reps[, 2], reps[, 3])),
             lw_med = median(lw[is.finite(lw)]),
             frac_usable = mean(lw > -700))
}

t0 <- Sys.time()
res <- do.call(rbind, mclapply(cells, run_cell, mc.cores = 12L,
                               mc.preschedule = FALSE))
cat(sprintf("main grid done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

## ---- bias anchor at q in {10, 20}: 20 coupled reps at N = 1e6 ----
anchors <- list()
for (q in c(10L, 20L)) {
  G <- er_graph(q, 4500L + q)
  e <- pick_edge(G)
  set.seed(4700L + q)
  ref <- vapply(1:20, function(r)
    akm_ratio_cpp(G, e[1], e[2], 0.5 * log(q), 1, 1, 1000000L)$log_J,
    numeric(1))
  anchors[[length(anchors) + 1L]] <- data.frame(
    q = q, logJ_ref = mean(ref), ref_se = sd(ref) / sqrt(20))
}
anchors <- do.call(rbind, anchors)

res$bias_vs_ref <- NA_real_
for (i in seq_len(nrow(anchors))) {
  sel <- res$q == anchors$q[i]
  res$bias_vs_ref[sel] <- res$mean_logJ[sel] - anchors$logJ_ref[i]
}
saveRDS(list(res = res, anchors = anchors),
        file.path(PROJECT_ROOT, "R", "results",
                  "hier_ratio_accuracy_at_scale_A2.rds"))

cat("\n== SD of log J-hat and coupling, by q and N (delta = 0.5 log q, eta = 1) ==\n")
s <- res[order(res$q, res$N, res$arm), ]
print(s, row.names = FALSE, digits = 3)
cat("\n== SD ratio independent / coupled ==\n")
sc <- res[res$arm == "coupled", ]; si <- res[res$arm == "independent", ]
si <- si[match(paste(sc$q, sc$N), paste(si$q, si$N)), ]
print(data.frame(sc[, c("q", "N")], sd_ratio = si$sd_logJ / sc$sd_logJ,
                 corr_coupled = sc$corr_means),
      row.names = FALSE, digits = 3)
cat("\n== bias anchors (q = 10, 20; reference se in parentheses) ==\n")
print(anchors, row.names = FALSE, digits = 4)
cat(sprintf("\nmax |bias_vs_ref| (coupled, q <= 20): %.4f\n",
            max(abs(res$bias_vs_ref[res$arm == "coupled"]), na.rm = TRUE)))
