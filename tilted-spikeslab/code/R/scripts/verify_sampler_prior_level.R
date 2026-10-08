# Prior-level verification of the R reference sampler: check V4 (complete chain
# at the prior level against exact enumeration) and check V5 (within-model
# update at the posterior level against importance sampling), in the numbering
# of the manuscript's verification appendix.
#
# The full chain (row-block sweep + edge toggles) is run WITHOUT data at q = 4
# and its graph marginal is compared against exact enumeration:
#   - joint specification: target p_J(Gamma) prop. to w(Gamma) Z(Gamma), with
#     Z from the J1 brute-force table (R/results/joint_spec_sparsening_J1_z.rds);
#   - hierarchical specification with EXACT log Z ratios from the same table:
#     target graph marginal is w(Gamma) itself (uniform at p_inc = 0.5).
# Two independent chain seeds per setting give the Monte Carlo noise floor.
#
# V4-lite: posterior within-model check at q = 3 (path graph, n = 10): moments
# from the within-only Gibbs chain against importance sampling with likelihood
# weights over exact prior draws.
#
# Output: R/results/verify_sampler_prior_level.rds (+ log)

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
source(file.path(PROJECT_ROOT, "R", "scripts", "sampler_reference.R"))

zt <- readRDS(file.path(PROJECT_ROOT, "R", "results", "joint_spec_sparsening_J1_z.rds"))

z_table <- function(q, eta, delta) {
  d <- zt[zt$q == q & zt$eta == eta & abs(zt$delta - delta) < 1e-9, ]
  s1 <- d[d$seed == min(d$seed), ]; s2 <- d[d$seed == max(d$seed), ]
  s1 <- s1[order(s1$mask), ]; s2 <- s2[order(s2$mask), ]
  (s1$z_hat + s2$z_hat) / 2
}

tv <- function(p, phat) 0.5 * sum(abs(p - phat))

run_setting <- function(label, spec, q, eta, delta, p_inc, n_iter, n_burn) {
  z <- z_table(q, eta, delta)
  M <- q * (q - 1) / 2
  m_of <- vapply(0:(2^M - 1L), function(mask)
    sum(bitwAnd(mask, 2^(0:(M - 1L))) != 0), numeric(1))
  logw <- m_of * log(p_inc) + (M - m_of) * log(1 - p_inc)
  p_exact <- if (spec == "joint") {
    pj <- exp(logw) * z; pj / sum(pj)
  } else {
    w <- exp(logw); w / sum(w)
  }
  chains <- lapply(c(9101L, 9102L), function(seed) {
    run_chain(q = q, S_mat = matrix(0, q, q), n = 0, delta = delta, beta = eta,
              sigma = 1, p_inc = p_inc, spec = spec, n_iter = n_iter,
              n_burn = n_burn, seed = seed,
              log_z = if (spec == "hierarchical") log(z) else NULL)
  })
  phat <- lapply(chains, function(ch) ch$visits / sum(ch$visits))
  pooled <- (phat[[1]] + phat[[2]]) / 2
  out <- data.frame(
    label = label, spec = spec, eta = eta, delta = delta, p_inc = p_inc,
    tv_exact = tv(p_exact, pooled),
    tv_seed1 = tv(p_exact, phat[[1]]), tv_seed2 = tv(p_exact, phat[[2]]),
    tv_between_seeds = tv(phat[[1]], phat[[2]]),
    max_abs_dev = max(abs(pooled - p_exact)),
    acc_rate = mean(vapply(chains, function(ch) ch$acc_rate, numeric(1))),
    E_m_exact = sum(m_of * p_exact), E_m_hat = sum(m_of * pooled))
  cat(sprintf("%-28s TV(exact) = %.4f  [seeds: %.4f / %.4f; between = %.4f]  E[m] %.3f vs %.3f  acc = %.2f\n",
              label, out$tv_exact, out$tv_seed1, out$tv_seed2,
              out$tv_between_seeds, out$E_m_hat, out$E_m_exact, out$acc_rate))
  out
}

cat("== full-chain prior-level checks at q = 4 (n_iter = 6e4) ==\n")
settings <- list(
  list("joint eta=1 delta=1 p=.5",  "joint",        4L, 1, 1, 0.5),
  list("joint eta=1 delta=0 p=.5",  "joint",        4L, 1, 0, 0.5),
  list("joint eta=2 delta=1 p=.25", "joint",        4L, 2, 1, 0.25),
  list("hier  eta=1 delta=1 p=.5",  "hierarchical", 4L, 1, 1, 0.5),
  list("hier  eta=1 delta=1 p=.25", "hierarchical", 4L, 1, 1, 0.25)
)
res <- do.call(rbind, lapply(settings, function(s)
  run_setting(s[[1]], s[[2]], s[[3]], s[[4]], s[[5]], s[[6]],
              n_iter = 6e4, n_burn = 4e3)))

## ---- V4-lite: posterior within-model check at q = 3 ----

cat("\n== V4-lite: posterior moments, within-only Gibbs vs importance sampling (q = 3, path graph, n = 10) ==\n")
set.seed(3301)
q <- 3L
G <- matrix(0L, q, q); G[1, 2] <- G[2, 1] <- 1L; G[2, 3] <- G[3, 2] <- 1L
K_true <- matrix(c(1.5, 0.4, 0, 0.4, 1.5, -0.5, 0, -0.5, 1.5), 3, 3)
n_obs <- 10L
Y <- matrix(stats::rnorm(n_obs * q), n_obs, q) %*% chol(solve(K_true))
S_mat <- crossprod(Y)
delta <- 1; eta <- 1

# posterior Gibbs (within-only, graph fixed)
K <- diag((delta + 1) / eta, q)
n_it <- 4e4; n_bu <- 2e3
post_draws <- matrix(0, n_it - n_bu, 4)
for (it in seq_len(n_it)) {
  for (l in seq_len(q)) K <- sample_row_block(K, G, l, S_mat, n_obs, delta, eta, 1)
  if (it > n_bu) {
    ev <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
    post_draws[it - n_bu, ] <- c(K[1, 2], K[2, 3], K[2, 2], min(ev))
  }
}
gibbs_mom <- colMeans(post_draws)

# importance sampling: prior draws from the (validated) prior Gibbs, likelihood weights
K <- diag((delta + 1) / eta, q)
n_pr <- 1.2e5
pri_stats <- matrix(0, n_pr, 4); pri_logw <- numeric(n_pr)
S0 <- matrix(0, q, q)
for (it in seq_len(n_pr)) {
  for (l in seq_len(q)) K <- sample_row_block(K, G, l, S0, 0, delta, eta, 1)
  ev <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  pri_stats[it, ] <- c(K[1, 2], K[2, 3], K[2, 2], min(ev))
  pri_logw[it] <- (n_obs / 2) * sum(log(ev)) - 0.5 * sum(S_mat * K)
}
w <- exp(pri_logw - max(pri_logw))
ess <- sum(w)^2 / sum(w^2)
is_mom <- colSums(pri_stats * w) / sum(w)
labels <- c("E[k12]", "E[k23]", "E[k22]", "E[lambda_min]")
for (a in 1:4) {
  cat(sprintf("%-14s Gibbs = %8.4f   IS = %8.4f   rel diff = %.4f\n",
              labels[a], gibbs_mom[a], is_mom[a],
              abs(gibbs_mom[a] / is_mom[a] - 1)))
}
cat(sprintf("IS effective sample size: %.0f of %d\n", ess, n_pr))

saveRDS(list(prior_level = res,
             posterior_check = data.frame(stat = labels, gibbs = gibbs_mom,
                                          importance = is_mom, ess = ess)),
        file.path(PROJECT_ROOT, "R", "results", "verify_sampler_prior_level.rds"))
cat("\nSaved R/results/verify_sampler_prior_level.rds\n")
