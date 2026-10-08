# J2b: extend the J2 sparse-regime (p_inc = 0.1) sampled study to q = 100.
# Same design as J2 (verified C++ chain, prior-only joint mode, two seeds,
# 4000 iterations, exact Sigma refresh each sweep).
#
# Output: R/results/joint_spec_sparsening_J2b_q100.rds (+ log)

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

q <- 100L
p_inc <- 0.1
cells <- list()
for (eta in c(0.5, 1, 2)) for (delta in c(0, 1, 2, 0.5 * log(q)))
  for (seed in c(8801L, 8802L))
    cells[[length(cells) + 1L]] <- list(eta = eta, delta = delta, seed = seed)
cat(sprintf("running %d chains on 12 cores\n", length(cells)))

run_cell <- function(cell) {
  set.seed(cell$seed)
  ch <- tryCatch(
    sampler_chain_cpp(S = matrix(0, q, q), n_obs = 0, q = q,
                      delta = cell$delta, beta = cell$eta, sigma = 1,
                      p_inc = p_inc, spec = 0L, log_z = numeric(0),
                      n_iter = 4000L, n_burn = 500L, refresh_every = 1L,
                      track_masks = FALSE),
    error = function(e) NULL)
  if (is.null(ch)) {
    return(data.frame(eta = cell$eta, delta = cell$delta, seed = cell$seed,
                      E_m = NA_real_, acc = NA_real_, drift = NA_real_))
  }
  data.frame(eta = cell$eta, delta = cell$delta, seed = cell$seed,
             E_m = ch$E_m, acc = ch$acc_rate, drift = ch$max_sigma_drift)
}

t0 <- Sys.time()
res <- do.call(rbind, mclapply(cells, run_cell, mc.cores = 12L,
                               mc.preschedule = FALSE))
n_fail <- sum(is.na(res$E_m))
if (n_fail > 0) {
  cat("failed chains:\n")
  print(res[is.na(res$E_m), c("eta", "delta", "seed")], row.names = FALSE)
}
res_ok <- res[!is.na(res$E_m), ]
cat(sprintf("done in %.1f min; %d/%d chains ok\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            nrow(res_ok), nrow(res)))

key <- interaction(res_ok$eta, res_ok$delta, drop = TRUE)
M <- q * (q - 1) / 2
agg <- do.call(rbind, lapply(split(res_ok, key), function(d) {
  E_J <- mean(d$E_m)
  pbar <- E_J / M
  data.frame(eta = d$eta[1], delta = d$delta[1], E_w = M * p_inc, E_J = E_J,
             halfdiff = if (nrow(d) >= 2) abs(diff(d$E_m))[1] / 2 else NA,
             ratio = E_J / (M * p_inc),
             odds_factor = (pbar / (1 - pbar)) / (p_inc / (1 - p_inc)),
             c0 = c0_fun(d$delta[1], d$eta[1]),
             max_drift = max(d$drift))
}))
rownames(agg) <- NULL
agg <- agg[order(agg$eta, agg$delta), ]
saveRDS(agg, file.path(PROJECT_ROOT, "R", "results",
                       "joint_spec_sparsening_J2b_q100.rds"))
cat("\n== q = 100, p_inc = 0.1: odds factor vs c(0) ==\n")
print(agg, row.names = FALSE, digits = 3)
