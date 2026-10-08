# Section 5.4, component A: bias and RMSE of the estimated log ratio, and the
# coupling correlation, over settings.
#
# Test pairs at q = 4: a cycle-closing addition (path 1-2, 2-3 plus edge
# (1,3)) and a disconnected-endpoint addition (same path plus edge (3,4),
# where the truth is 1/c(0) analytically). Factorial: eta in {0.5, 1, 2},
# delta in {0, 1, 2}, N in {1e2, 1e3, 1e4}, arms coupled and independent,
# R = 200 replicate estimates per cell. Truth for the cycle-closing pair:
# 20 coupled estimates at N = 1e6 (the estimator is validated against brute
# force in verify_akm_ratio.R);
# for the disconnected pair: analytic.
#
# Output: R/results/hier_ratio_accuracy_A.rds (+ log)

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})
PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
sourceCpp(file.path(PROJECT_ROOT, "R", "src", "spikeslab_sampler.cpp"))

c0_fun <- function(delta, eta) {
  s_int <- function(s) {
    t <- (delta + 1) * s
    exp(delta * log(t) - t - lgamma(delta + 1) + log(delta + 1)) /
      sqrt(1 + 2 * eta^2 / t)
  }
  stats::integrate(s_int, 0, Inf, rel.tol = 1e-11)$value
}

q <- 4L
G_path <- matrix(0L, q, q)
G_path[1, 2] <- G_path[2, 1] <- 1L
G_path[2, 3] <- G_path[3, 2] <- 1L
pairs <- list(cycle = c(1L, 3L), disconnected = c(3L, 4L))

etas <- c(0.5, 1, 2)
deltas <- c(0, 1, 2)
Ns <- c(100L, 1000L, 10000L)
R_reps <- 200L

cells <- list()
for (pair in names(pairs)) for (eta in etas) for (delta in deltas)
  cells[[length(cells) + 1L]] <- list(pair = pair, eta = eta, delta = delta)

run_cell <- function(cell) {
  ij <- pairs[[cell$pair]]
  set.seed(5501L)
  logJ_true <- if (cell$pair == "disconnected") {
    -log(c0_fun(cell$delta, cell$eta))
  } else {
    mean(vapply(1:20, function(r)
      akm_ratio_cpp(G_path, ij[1], ij[2], cell$delta, cell$eta, 1, 1000000L)$log_J,
      numeric(1)))
  }
  rows <- list()
  for (N in Ns) for (arm in c("coupled", "independent")) {
    reps <- t(vapply(seq_len(R_reps), function(r) {
      est <- akm_ratio_cpp(G_path, ij[1], ij[2], cell$delta, cell$eta, 1, N,
                           coupled = (arm == "coupled"))
      c(est$log_J, est$mean_minus, est$mean_plus)
    }, numeric(3)))
    bias <- mean(reps[, 1]) - logJ_true
    rows[[length(rows) + 1L]] <- data.frame(
      pair = cell$pair, eta = cell$eta, delta = cell$delta, N = N, arm = arm,
      logJ_true = logJ_true, bias = bias,
      bias_mcse = sd(reps[, 1]) / sqrt(R_reps),
      rmse = sqrt(mean((reps[, 1] - logJ_true)^2)),
      corr_means = cor(reps[, 2], reps[, 3]))
  }
  do.call(rbind, rows)
}

t0 <- Sys.time()
res <- do.call(rbind, mclapply(cells, run_cell, mc.cores = 12L,
                               mc.preschedule = FALSE))
cat(sprintf("done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
saveRDS(res, file.path(PROJECT_ROOT, "R", "results", "hier_ratio_accuracy_A.rds"))

cat("\n== RMSE of log J-hat (coupled), by cell ==\n")
s <- res[res$arm == "coupled", ]
print(s[order(s$pair, s$eta, s$delta, s$N),
        c("pair", "eta", "delta", "N", "logJ_true", "bias", "bias_mcse",
          "rmse", "corr_means")],
      row.names = FALSE, digits = 3)

cat("\n== variance reduction: RMSE(independent) / RMSE(coupled) ==\n")
sc <- res[res$arm == "coupled", ]
si <- res[res$arm == "independent", ]
key <- function(d) paste(d$pair, d$eta, d$delta, d$N)
si <- si[match(key(sc), key(si)), ]
vr <- data.frame(sc[, c("pair", "eta", "delta", "N")],
                 rmse_ratio = si$rmse / sc$rmse,
                 corr_coupled = sc$corr_means)
print(vr[order(vr$pair, vr$eta, vr$delta, vr$N), ], row.names = FALSE,
      digits = 3)

cat("\n== bias detectability: |bias|/mcse (coupled) ==\n")
s$z <- abs(s$bias) / s$bias_mcse
cat(sprintf("cells with |bias| > 3 mcse: %d of %d; max |bias| = %.4f (N = %d)\n",
            sum(s$z > 3), nrow(s), max(abs(s$bias)), s$N[which.max(abs(s$bias))]))
