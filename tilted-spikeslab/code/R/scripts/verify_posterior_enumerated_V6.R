# V6: the full chain (C++ implementation) against the exactly enumerated
# posterior over graphs at q = 4, both specifications.
#
# Enumerated reference: for each of the 64 graphs, the marginal likelihood
# factor m_lik(Gamma) = E_{p(K | Gamma)}[L(Y | K)] is estimated by fixed-graph
# prior Gibbs draws (verified row-block sampler, R reference) with likelihood
# weights, two independent reference seeds; the posterior over graphs is then
#   joint:        p(Gamma | Y) prop. to w(Gamma) Z(Gamma) m_lik(Gamma)
#   hierarchical: p(Gamma | Y) prop. to w(Gamma) m_lik(Gamma)
# with Z(Gamma) from the J1 brute-force table (delta = 1: reliable regime).
#
# Chains: C++ sampler with data (n = 10), two seeds, 6e4 iterations; joint
# spec and hierarchical spec with exact Z-ratios. PASS criterion as in V4:
# TV(pooled, enumerated) at or below max(between-seed floor, reference noise).
#
# Settings: (eta = 1, delta = 1, p_inc = 0.5) and (eta = 2, delta = 1,
# p_inc = 0.25).
#
# Output: R/results/verify_posterior_enumerated_V6.rds (+ log)

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})
PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
sourceCpp(file.path(PROJECT_ROOT, "R", "src", "spikeslab_sampler.cpp"))
source(file.path(PROJECT_ROOT, "R", "scripts", "sampler_reference.R"))

q <- 4L
n_obs <- 10
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
tv <- function(p, phat) 0.5 * sum(abs(p - phat))
logsumexp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

## ---- data ----
set.seed(6001L)
K_true <- diag(1, q)
K_true[cbind(1:3, 2:4)] <- -0.4
K_true[cbind(2:4, 1:3)] <- -0.4
Y <- matrix(rnorm(n_obs * q), n_obs, q) %*% chol(solve(K_true))
S <- crossprod(Y)

## ---- z table (brute force, q = 4, delta = 1) ----
zt <- readRDS(file.path(PROJECT_ROOT, "R", "results", "joint_spec_sparsening_J1_z.rds"))
log_z_of <- function(eta) {
  d <- zt[zt$q == 4 & zt$eta == eta & abs(zt$delta - 1) < 1e-9, ]
  s1 <- d[d$seed == min(d$seed), ]; s2 <- d[d$seed == max(d$seed), ]
  s1 <- s1[order(s1$mask), ]; s2 <- s2[order(s2$mask), ]
  log((s1$z_hat + s2$z_hat) / 2)
}

## ---- reference: log m_lik per graph via prior Gibbs + likelihood weights ----
S0 <- matrix(0, q, q)
# n_keep = 60000L prior Gibbs draws per graph and reference seed. The archived
# run of 2026-07-28 used this value (a rerun on 2026-10-08 reproduced it to
# four decimals); its log printed a stale "(of 15000)" string, corrected below.
ref_cell <- function(mask, eta, seed, n_keep = 60000L, n_burn = 500L) {
  set.seed(seed)
  G <- graph_of_mask(mask)
  K <- diag(2 / eta, q)
  lw <- numeric(n_keep)
  for (it in seq_len(n_burn + n_keep)) {
    for (l in seq_len(q))
      K <- sample_row_block(K, G, l, S0, 0, delta = 1, beta = eta, sigma = 1)
    if (it > n_burn) {
      ev <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
      lw[it - n_burn] <- (n_obs / 2) * sum(log(ev)) - 0.5 * sum(S * K)
    }
  }
  ess <- exp(2 * logsumexp(lw) - logsumexp(2 * lw))
  c(log_mlik = logsumexp(lw) - log(n_keep), ess = ess)
}

settings <- list(list(eta = 1, p_inc = 0.5), list(eta = 2, p_inc = 0.25))
ref_seeds <- c(6101L, 6102L)

results <- list()
for (st in settings) {
  eta <- st$eta; p_inc <- st$p_inc
  cat(sprintf("\n== setting eta = %.0f, delta = 1, p_inc = %.2f ==\n", eta, p_inc))
  jobs <- expand.grid(mask = masks, seed = ref_seeds)
  ref <- mclapply(seq_len(nrow(jobs)), function(r)
    ref_cell(jobs$mask[r], eta, jobs$seed[r]),
    mc.cores = 12L, mc.preschedule = TRUE)
  ref <- do.call(rbind, ref)
  lm1 <- ref[jobs$seed == ref_seeds[1], "log_mlik"]
  lm2 <- ref[jobs$seed == ref_seeds[2], "log_mlik"]
  cat(sprintf("reference min ESS = %.0f (of 60000)\n", min(ref[, "ess"])))

  logw <- m_of * log(p_inc) + (M - m_of) * log(1 - p_inc)
  lz <- log_z_of(eta)
  targets <- list(
    joint = list(l1 = logw + lz + lm1, l2 = logw + lz + lm2),
    hierarchical = list(l1 = logw + lm1, l2 = logw + lm2))

  for (spec_name in names(targets)) {
    tg <- targets[[spec_name]]
    p1 <- exp(tg$l1 - logsumexp(tg$l1)); p2 <- exp(tg$l2 - logsumexp(tg$l2))
    p_ref <- (p1 + p2) / 2
    ref_noise <- tv(p1, p2)
    spec_code <- if (spec_name == "joint") 0L else 1L
    chains <- lapply(c(6201L, 6202L), function(seed) {
      set.seed(seed)
      sampler_chain_cpp(S = S, n_obs = n_obs, q = q, delta = 1, beta = eta,
                        sigma = 1, p_inc = p_inc, spec = spec_code,
                        log_z = lz, n_iter = 6e4L, n_burn = 4e3L,
                        refresh_every = 50L, track_masks = TRUE)
    })
    phat <- lapply(chains, function(ch) ch$mask_visits / sum(ch$mask_visits))
    pooled <- (phat[[1]] + phat[[2]]) / 2
    row <- data.frame(
      eta = eta, p_inc = p_inc, spec = spec_name,
      tv_exact = tv(p_ref, pooled),
      tv_between_seeds = tv(phat[[1]], phat[[2]]),
      tv_reference_noise = ref_noise,
      E_m_ref = sum(m_of * p_ref),
      E_m_hat = mean(vapply(chains, function(ch) ch$E_m, numeric(1))),
      min_ess_ref = min(ref[, "ess"]))
    results[[length(results) + 1L]] <- row
    cat(sprintf("%-13s TV(enum) = %.4f  [seeds %.4f, ref noise %.4f]  E[m] %.3f vs %.3f\n",
                spec_name, row$tv_exact, row$tv_between_seeds,
                row$tv_reference_noise, row$E_m_hat, row$E_m_ref))
  }
}
res <- do.call(rbind, results)
saveRDS(list(res = res, S = S, K_true = K_true),
        file.path(PROJECT_ROOT, "R", "results", "verify_posterior_enumerated_V6.rds"))
cat("\nSaved R/results/verify_posterior_enumerated_V6.rds\n")
