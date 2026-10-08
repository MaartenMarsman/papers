# J2: sampled consequences of the joint specification at q = 5, 10, 20, 50
# (Section 5.3), using the verified C++ sampler in prior-only mode.
#
# Factors: q in {5, 10, 20, 50}; eta in {0.5, 1, 2}; delta in {0, 1, 2, 0.5 log q};
# p_inc in {0.1, 0.25, 0.5}; two independent chain seeds per cell.
# Estimands: E_J[m] / E_w[m]; the implied average per-edge odds factor against
# c(0); at q = 5 the agreement with the exact J1 enumeration doubles as the
# joint chain's implementation check at scale.
#
# Output: R/results/joint_spec_sparsening_J2.rds (+ log)

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

qs <- c(5L, 10L, 20L, 50L)
etas <- c(0.5, 1, 2)
p_incs <- c(0.1, 0.25, 0.5)
seeds <- c(8801L, 8802L)

cells <- list()
for (q in qs) for (eta in etas) for (delta in c(0, 1, 2, 0.5 * log(q)))
  for (p in p_incs) for (seed in seeds) {
    cells[[length(cells) + 1L]] <- list(q = q, eta = eta, delta = delta,
                                        p_inc = p, seed = seed)
  }
cat(sprintf("running %d chains on 12 cores\n", length(cells)))

run_cell <- function(cell) {
  set.seed(cell$seed)
  q <- cell$q
  ch <- tryCatch(
    sampler_chain_cpp(S = matrix(0, q, q), n_obs = 0, q = q,
                      delta = cell$delta, beta = cell$eta, sigma = 1,
                      p_inc = cell$p_inc, spec = 0L, log_z = numeric(0),
                      n_iter = 4000L, n_burn = 500L, refresh_every = 1L,
                      track_masks = FALSE),
    error = function(e) NULL)
  if (is.null(ch)) {
    return(data.frame(q = q, eta = cell$eta, delta = cell$delta,
                      p_inc = cell$p_inc, seed = cell$seed, E_m = NA_real_,
                      acc = NA_real_, drift = NA_real_))
  }
  data.frame(q = q, eta = cell$eta, delta = cell$delta, p_inc = cell$p_inc,
             seed = cell$seed, E_m = ch$E_m, acc = ch$acc_rate,
             drift = ch$max_sigma_drift)
}

t0 <- Sys.time()
out <- mclapply(cells, run_cell, mc.cores = 12L, mc.preschedule = FALSE)
res <- do.call(rbind, out[vapply(out, is.data.frame, logical(1))])
n_fail <- sum(is.na(res$E_m))
if (n_fail > 0) {
  cat("failed chains:\n")
  print(res[is.na(res$E_m), c("q", "eta", "delta", "p_inc", "seed")],
        row.names = FALSE)
}
res <- res[!is.na(res$E_m), ]
cat(sprintf("done in %.1f min; %d/%d chains ok; max drift %.2e\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            nrow(res), length(cells), max(res$drift)))

## ---- aggregate ----
key <- interaction(res$q, res$eta, res$delta, res$p_inc, drop = TRUE)
agg <- do.call(rbind, lapply(split(res, key), function(d) {
  M <- d$q[1] * (d$q[1] - 1) / 2
  p <- d$p_inc[1]
  E_J <- mean(d$E_m)
  pbar <- E_J / M
  halfdiff <- if (nrow(d) >= 2) abs(diff(d$E_m))[1] / 2 else NA_real_
  data.frame(q = d$q[1], eta = d$eta[1], delta = d$delta[1], p_inc = p,
             E_w = M * p, E_J = E_J, E_J_halfdiff = halfdiff,
             ratio = E_J / (M * p),
             odds_factor = (pbar / (1 - pbar)) / (p / (1 - p)),
             c0 = c0_fun(d$delta[1], d$eta[1]),
             acc = mean(d$acc))
}))
rownames(agg) <- NULL
saveRDS(agg, file.path(PROJECT_ROOT, "R", "results", "joint_spec_sparsening_J2.rds"))

## ---- J1 agreement at q = 5 ----
j1 <- readRDS(file.path(PROJECT_ROOT, "R", "results", "joint_spec_sparsening_J1.rds"))$summary
cat("\n== implementation check: J2 sampled vs J1 exact at q = 5 ==\n")
for (i in which(agg$q == 5)) {
  j <- j1[j1$q == 5 & j1$eta == agg$eta[i] & abs(j1$delta - agg$delta[i]) < 1e-9 &
            j1$p_inc == agg$p_inc[i], ]
  if (nrow(j) == 1) {
    cat(sprintf("eta=%.1f delta=%.2f p=%.2f: sampled E[m] = %.3f (+-%.3f), exact = %.3f, diff = %+.3f\n",
                agg$eta[i], agg$delta[i], agg$p_inc[i],
                agg$E_J[i], agg$E_J_halfdiff[i], j$E_J, agg$E_J[i] - j$E_J))
  }
}

## ---- drift of the average odds factor from c(0) with dimension ----
cat("\n== odds factor vs c(0) across dimension (p_inc = 0.1, sparse regime) ==\n")
s <- agg[agg$p_inc == 0.1, ]
s <- s[order(s$eta, s$delta, s$q), ]
print(s[, c("q", "eta", "delta", "ratio", "odds_factor", "c0", "acc")],
      row.names = FALSE, digits = 3)
cat("\n== odds factor vs c(0) (p_inc = 0.5, dense regime) ==\n")
s <- agg[agg$p_inc == 0.5, ]
s <- s[order(s$eta, s$delta, s$q), ]
print(s[, c("q", "eta", "delta", "ratio", "odds_factor", "c0", "acc")],
      row.names = FALSE, digits = 3)
