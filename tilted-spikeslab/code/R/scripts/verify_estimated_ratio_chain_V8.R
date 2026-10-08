# V8: the hierarchical plug-in chain with ESTIMATED ratios (coupled AKM
# estimator, N auxiliary draws per candidate) against the exact target at
# q = 4, as N grows. The exact hierarchical prior-level target over graphs is
# the Bernoulli law w(Gamma) itself; the exact-ratio chain (V4) provides the
# attainable noise floor.
#
# Settings match V4's hierarchical arm: eta = 1, delta = 1, p in {0.5, 0.25}.
# N in {100, 1000, 10000}; two chain seeds per cell.
#
# Output: R/results/verify_estimated_ratio_chain_V8.rds (+ log)

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})
PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
sourceCpp(file.path(PROJECT_ROOT, "R", "src", "spikeslab_sampler.cpp"))

q <- 4L
M <- 6L
m_of <- vapply(0:(2^M - 1L), function(mask)
  sum(bitwAnd(mask, 2^(0:(M - 1L))) != 0), numeric(1))
tv <- function(p, phat) 0.5 * sum(abs(p - phat))

cells <- list()
for (p_inc in c(0.5, 0.25)) for (N in c(100L, 1000L, 10000L))
  for (seed in c(9301L, 9302L)) {
    cells[[length(cells) + 1L]] <- list(p_inc = p_inc, N = N, seed = seed,
                                        n_iter = if (N == 10000L) 2e4L else 6e4L)
  }
cat(sprintf("running %d chains on 12 cores\n", length(cells)))

run_cell <- function(cell) {
  set.seed(cell$seed)
  t0 <- Sys.time()
  ch <- sampler_chain_cpp(S = matrix(0, q, q), n_obs = 0, q = q, delta = 1,
                          beta = 1, sigma = 1, p_inc = cell$p_inc, spec = 2L,
                          log_z = numeric(0), n_iter = cell$n_iter,
                          n_burn = 4000L, refresh_every = 50L,
                          track_masks = TRUE, n_aux = cell$N)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(p_inc = cell$p_inc, N = cell$N, seed = cell$seed,
       phat = ch$mask_visits / sum(ch$mask_visits), E_m = ch$E_m,
       acc = ch$acc_rate, secs = el, n_iter = cell$n_iter)
}
out <- mclapply(cells, run_cell, mc.cores = 12L, mc.preschedule = FALSE)

rows <- list()
for (p_inc in c(0.5, 0.25)) for (N in c(100L, 1000L, 10000L)) {
  sel <- Filter(function(o) o$p_inc == p_inc && o$N == N, out)
  logw <- m_of * log(p_inc) + (M - m_of) * log(1 - p_inc)
  p_exact <- exp(logw) / sum(exp(logw))
  pooled <- (sel[[1]]$phat + sel[[2]]$phat) / 2
  rows[[length(rows) + 1L]] <- data.frame(
    p_inc = p_inc, N = N, n_iter = sel[[1]]$n_iter,
    tv_exact = tv(p_exact, pooled),
    tv_between_seeds = tv(sel[[1]]$phat, sel[[2]]$phat),
    E_m_exact = sum(m_of * p_exact),
    E_m_hat = mean(c(sel[[1]]$E_m, sel[[2]]$E_m)),
    acc = mean(c(sel[[1]]$acc, sel[[2]]$acc)),
    secs_per_chain = mean(c(sel[[1]]$secs, sel[[2]]$secs)))
}
res <- do.call(rbind, rows)
cat("\n== V8: plug-in chain vs exact Bernoulli target (eta = 1, delta = 1) ==\n")
print(res, row.names = FALSE, digits = 3)
cat("\nReference (V4, exact-ratio chain): TV(exact) = 0.015 / 0.007 at\n")
cat("p = 0.5 / 0.25 with between-seed floors 0.034 / 0.014 (6e4 iterations).\n")

saveRDS(res, file.path(PROJECT_ROOT, "R", "results",
                       "verify_estimated_ratio_chain_V8.rds"))
cat("\nSaved R/results/verify_estimated_ratio_chain_V8.rds\n")
