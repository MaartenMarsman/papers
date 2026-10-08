# V7: simulation-based calibration of the full chain, both specifications,
# prior-only (n = 0, posterior = prior) and posterior (n = 10) tiers.
#
# Per replication: draw a graph from the specification's exact graph prior
# (joint: p_J prop. to w Z with brute-force Z at delta = 1, the reliable
# regime; hierarchical: Bernoulli w directly), draw K | graph by the verified
# row-block Gibbs sampler (500 sweeps), simulate data when n > 0, run the
# C++ chain (hierarchical with exact Z-ratios), and record the rank of each
# test statistic of the prior draw among L = 99 thinned chain draws
# (randomized tie-breaking). Under a correct sampler the ranks are uniform
# on {0, ..., 99}.
#
# Statistics: edge count m, k_11, k_12, log|K|. R = 400 replications per
# (specification, tier); chi-squared over 20 rank bins (expected 20 per bin)
# and the maximum ecdf deviation against the 95% Kolmogorov band.
#
# Setting: q = 4, eta = 1, delta = 1, p_inc = 0.5.
#
# Output: R/results/verify_sbc_V7.rds (+ log)

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})
PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
sourceCpp(file.path(PROJECT_ROOT, "R", "src", "spikeslab_sampler.cpp"))
source(file.path(PROJECT_ROOT, "R", "scripts", "sampler_reference.R"))

q <- 4L
M <- 6L
edges4 <- t(utils::combn(q, 2))
masks <- 0:(2^M - 1L)
m_of <- vapply(masks, function(mask)
  sum(bitwAnd(mask, 2^(0:(M - 1L))) != 0), numeric(1))
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

zt <- readRDS(file.path(PROJECT_ROOT, "R", "results", "joint_spec_sparsening_J1_z.rds"))
d <- zt[zt$q == 4 & zt$eta == 1 & abs(zt$delta - 1) < 1e-9, ]
s1 <- d[d$seed == min(d$seed), ]; s2 <- d[d$seed == max(d$seed), ]
s1 <- s1[order(s1$mask), ]; s2 <- s2[order(s2$mask), ]
lz <- log((s1$z_hat + s2$z_hat) / 2)

logw <- m_of * log(0.5) + (M - m_of) * log(0.5)
pJ <- exp(logw + lz); pJ <- pJ / sum(pJ)
pW <- exp(logw); pW <- pW / sum(pW)

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

run_rep <- function(rep, spec, n_obs) {
  seed <- 770000L + (n_obs > 0) * 100000L +
    (spec == "hierarchical") * 10000L + rep
  set.seed(seed)
  mask_t <- sample(masks, 1L, prob = if (spec == "joint") pJ else pW)
  Gt <- graph_of_mask(mask_t)
  Kt <- prior_K_draw(Gt)
  if (n_obs > 0) {
    Sig <- solve(Kt)
    Y <- matrix(rnorm(n_obs * q), n_obs, q) %*% chol((Sig + t(Sig)) / 2)
    S <- crossprod(Y)
  } else {
    S <- S0
  }
  ch <- sampler_chain_cpp(S = S, n_obs = n_obs, q = q, delta = 1, beta = 1,
                          sigma = 1, p_inc = 0.5,
                          spec = if (spec == "joint") 0L else 1L,
                          log_z = lz, n_iter = n_iter, n_burn = burn,
                          refresh_every = 50L, track_masks = FALSE,
                          n_aux = 0L, stat_thin = thin)
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
for (n_obs in c(0, 10)) for (spec in c("joint", "hierarchical")) {
  t0 <- Sys.time()
  ranks <- do.call(rbind, mclapply(seq_len(R_reps), run_rep, spec = spec,
                                   n_obs = n_obs, mc.cores = 12L,
                                   mc.preschedule = TRUE))
  for (stat in colnames(ranks)) {
    a <- analyze(ranks[, stat])
    rows[[length(rows) + 1L]] <- data.frame(
      tier = if (n_obs == 0) "prior" else "posterior", spec = spec,
      stat = stat, chisq = a["chisq"], p = a["p"], ks = a["ks"],
      ks_ok = a["ks"] < ks_band)
  }
  cat(sprintf("%s / %-13s done in %.1f min\n",
              if (n_obs == 0) "prior" else "posterior", spec,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
res <- do.call(rbind, rows)
rownames(res) <- NULL
cat(sprintf("\n== V7 SBC (R = %d, L = %d, 20 bins, KS band %.4f) ==\n",
            R_reps, L, ks_band))
print(res, row.names = FALSE, digits = 3)
cat(sprintf("\nmin p = %.4f over %d tests; all KS in band: %s\n",
            min(res$p), nrow(res), all(res$ks_ok)))
saveRDS(res, file.path(PROJECT_ROOT, "R", "results", "verify_sbc_V7.rds"))
cat("Saved R/results/verify_sbc_V7.rds\n")
