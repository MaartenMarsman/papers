# Section 5.4, component C: simulation-based calibration of the hierarchical
# PLUG-IN chain as a function of the auxiliary sample size N, posterior tier
# (n = 10), q = 4, eta = 1, delta = 1, p_inc = 0.5.
#
# Same SBC machinery as V7 (400 replications; ranks of four statistics among
# 99 thinned draws; exact Bernoulli graph prior draws; verified row-block
# Gibbs for K | graph), but the chain estimates every Z-ratio with the
# coupled estimator at N in {2, 10, 100, 1000}. Deliberately tiny N included
# to locate where calibration visibly breaks.
#
# Output: R/results/hier_sbc_vs_N_C.rds (+ log)

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})
PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
sourceCpp(file.path(PROJECT_ROOT, "R", "src", "spikeslab_sampler.cpp"))
source(file.path(PROJECT_ROOT, "R", "scripts", "sampler_reference.R"))

q <- 4L
M <- 6L
n_obs <- 10L
edges4 <- t(utils::combn(q, 2))
masks <- 0:(2^M - 1L)
graph_of_mask <- function(mask) {
  G <- matrix(0L, q, q)
  for (k in seq_len(M)) {
    if (bitwAnd(mask, bitwShiftL(1L, k - 1L)) != 0L) {
      G[edges4[k, 1], edges4[k, 2]] <- 1L
      G[edges4[k, 2], edges4[k, 1]] <- 1L
    }
  }
  G
}

L <- 99L; thin <- 40L; burn <- 500L
n_iter <- burn + L * thin
R_reps <- 400L
S0 <- matrix(0, q, q)

prior_K_draw <- function(G, n_sweeps = 500L) {
  K <- diag(2, q)
  for (s in seq_len(n_sweeps)) for (l in seq_len(q))
    K <- sample_row_block(K, G, l, S0, 0, delta = 1, beta = 1, sigma = 1)
  K
}
rank_of <- function(g, post)
  sum(post < g) + sample.int(sum(post == g) + 1L, 1L) - 1L

run_rep <- function(rep, N) {
  set.seed(880000L + N * 1000L + rep)
  mask_t <- sample(masks, 1L)                    # Bernoulli(0.5): uniform
  Gt <- graph_of_mask(mask_t)
  Kt <- prior_K_draw(Gt)
  Sig <- solve(Kt)
  Y <- matrix(rnorm(n_obs * q), n_obs, q) %*% chol((Sig + t(Sig)) / 2)
  ch <- sampler_chain_cpp(S = crossprod(Y), n_obs = n_obs, q = q, delta = 1,
                          beta = 1, sigma = 1, p_inc = 0.5, spec = 2L,
                          log_z = numeric(0), n_iter = n_iter, n_burn = burn,
                          refresh_every = 50L, track_masks = FALSE,
                          n_aux = N, stat_thin = thin)
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
for (N in c(2L, 10L, 100L, 1000L)) {
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
cat(sprintf("\n== plug-in SBC vs N (R = %d, KS band %.4f) ==\n", R_reps, ks_band))
print(res, row.names = FALSE, digits = 3)
saveRDS(res, file.path(PROJECT_ROOT, "R", "results", "hier_sbc_vs_N_C.rds"))
cat("Saved R/results/hier_sbc_vs_N_C.rds\n")
