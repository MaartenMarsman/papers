# Section 5.4, component C2: simulation-based calibration of the plug-in
# chain AT SCALE (q = 10), posterior tier, as a function of N.
#
# Same machinery as component C but with Bernoulli edge draws for the prior
# graph (no enumeration needed at any q) and the verified row-block Gibbs for
# K | graph. N in {1e2, 1e3}; A2 predicts SD(log J-hat) of about 0.8 and 0.26
# there, so visible miscalibration is expected at N = 1e2 and marginal
# behavior at N = 1e3.
#
# Setting: q = 10, eta = 1, delta = 0.5 log q, p_inc = 0.2, n = 20.
#
# Output: R/results/hier_sbc_at_scale_C2.rds (+ log)

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})
PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
sourceCpp(file.path(PROJECT_ROOT, "R", "src", "spikeslab_sampler.cpp"))
source(file.path(PROJECT_ROOT, "R", "scripts", "sampler_reference.R"))

q <- 10L
n_obs <- 20L
p_inc <- 0.2
delta <- 0.5 * log(q)
L <- 99L; thin <- 40L; burn <- 500L
n_iter <- burn + L * thin
R_reps <- 400L
S0 <- matrix(0, q, q)

prior_graph_draw <- function() {
  G <- matrix(0L, q, q)
  for (a in 1:(q - 1)) for (b in (a + 1):q)
    if (runif(1) < p_inc) { G[a, b] <- 1L; G[b, a] <- 1L }
  G
}
prior_K_draw <- function(G, n_sweeps = 400L) {
  K <- diag((delta + 1), q)
  for (s in seq_len(n_sweeps)) for (l in seq_len(q))
    K <- sample_row_block(K, G, l, S0, 0, delta = delta, beta = 1, sigma = 1)
  K
}
rank_of <- function(g, post)
  sum(post < g) + sample.int(sum(post == g) + 1L, 1L) - 1L

run_rep <- function(rep, N) {
  set.seed(990000L + N * 100L + rep)
  Gt <- prior_graph_draw()
  Kt <- prior_K_draw(Gt)
  Sig <- solve(Kt)
  Y <- matrix(rnorm(n_obs * q), n_obs, q) %*% chol((Sig + t(Sig)) / 2)
  ch <- sampler_chain_cpp(S = crossprod(Y), n_obs = n_obs, q = q,
                          delta = delta, beta = 1, sigma = 1, p_inc = p_inc,
                          spec = 2L, log_z = numeric(0), n_iter = n_iter,
                          n_burn = burn, refresh_every = 50L,
                          track_masks = FALSE, n_aux = N, stat_thin = thin)
  c(m = rank_of(sum(Gt) / 2, ch$rec_m),
    k11 = rank_of(Kt[1, 1], ch$rec_k11),
    k12 = rank_of(Kt[1, 2], ch$rec_k12),
    logdet = rank_of(as.numeric(determinant(Kt)$modulus), ch$rec_logdet))
}

analyze <- function(ranks) {
  bins <- tabulate(ranks %/% 5 + 1L, nbins = 20L)
  chisq <- sum((bins - R_reps / 20)^2 / (R_reps / 20))
  u <- sort((ranks + 0.5) / (L + 1))
  ks <- max(pmax(abs(u - seq_len(R_reps) / R_reps),
                 abs(u - (seq_len(R_reps) - 1) / R_reps)))
  c(chisq = chisq, p = pchisq(chisq, df = 19, lower.tail = FALSE), ks = ks)
}

ks_band <- 1.358 / sqrt(R_reps)
rows <- list()
for (N in c(100L, 1000L)) {
  t0 <- Sys.time()
  ranks <- do.call(rbind, mclapply(seq_len(R_reps), run_rep, N = N,
                                   mc.cores = 12L, mc.preschedule = TRUE))
  for (stat in colnames(ranks)) {
    a <- analyze(ranks[, stat])
    rows[[length(rows) + 1L]] <- data.frame(
      N = N, stat = stat, chisq = a["chisq"], p = a["p"], ks = a["ks"],
      ks_ok = a["ks"] < ks_band)
  }
  cat(sprintf("N = %5d done in %.1f min\n", N,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
res <- do.call(rbind, rows)
rownames(res) <- NULL
cat(sprintf("\n== plug-in SBC at q = 10 (R = %d, KS band %.4f) ==\n",
            R_reps, ks_band))
print(res, row.names = FALSE, digits = 3)
saveRDS(res, file.path(PROJECT_ROOT, "R", "results", "hier_sbc_at_scale_C2.rds"))
cat("Saved R/results/hier_sbc_at_scale_C2.rds\n")
