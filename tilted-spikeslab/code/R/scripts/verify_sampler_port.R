# Port verification: the C++ sampler against exact enumeration (same targets
# as the R-reference verification) plus the Sigma-drift diagnostic and a speed
# reading. PASS criterion: TV(exact) at or below the between-seed noise floor,
# as for the R reference, and negligible Sigma drift.
#
# Output: R/results/verify_sampler_port.rds (+ log)

suppressPackageStartupMessages(library(Rcpp))
PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
sourceCpp(file.path(PROJECT_ROOT, "R", "src", "spikeslab_sampler.cpp"))

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
  spec_code <- if (spec == "joint") 0L else 1L
  chains <- lapply(c(9201L, 9202L), function(seed) {
    set.seed(seed)
    sampler_chain_cpp(S = matrix(0, q, q), n_obs = 0, q = q, delta = delta,
                      beta = eta, sigma = 1, p_inc = p_inc, spec = spec_code,
                      log_z = log(z), n_iter = n_iter, n_burn = n_burn,
                      refresh_every = 50L, track_masks = TRUE)
  })
  phat <- lapply(chains, function(ch) ch$mask_visits / sum(ch$mask_visits))
  pooled <- (phat[[1]] + phat[[2]]) / 2
  drift <- max(vapply(chains, function(ch) ch$max_sigma_drift, numeric(1)))
  out <- data.frame(
    label = label, spec = spec, eta = eta, delta = delta, p_inc = p_inc,
    tv_exact = tv(p_exact, pooled),
    tv_between_seeds = tv(phat[[1]], phat[[2]]),
    E_m_exact = sum(m_of * p_exact),
    E_m_hat = mean(vapply(chains, function(ch) ch$E_m, numeric(1))),
    acc_rate = mean(vapply(chains, function(ch) ch$acc_rate, numeric(1))),
    max_sigma_drift = drift)
  cat(sprintf("%-28s TV(exact) = %.4f  [between-seed = %.4f]  E[m] %.3f vs %.3f  acc = %.2f  drift = %.2e\n",
              label, out$tv_exact, out$tv_between_seeds,
              out$E_m_hat, out$E_m_exact, out$acc_rate, drift))
  out
}

cat("== C++ port vs exact enumeration at q = 4 (n_iter = 6e4) ==\n")
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

cat("\n== speed reading (prior-only joint chain) ==\n")
for (q in c(20L, 50L)) {
  set.seed(1)
  t0 <- Sys.time()
  ch <- sampler_chain_cpp(S = matrix(0, q, q), n_obs = 0, q = q, delta = 1,
                          beta = 1, sigma = 1, p_inc = 0.1, spec = 0L,
                          log_z = numeric(0), n_iter = 2000L, n_burn = 500L,
                          refresh_every = 50L, track_masks = FALSE)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("q = %2d: %.0f iterations/s (E[m] = %.1f, acc = %.2f, drift = %.2e)\n",
              q, 2000 / el, ch$E_m, ch$acc_rate, ch$max_sigma_drift))
}

saveRDS(res, file.path(PROJECT_ROOT, "R", "results", "verify_sampler_port.rds"))
cat("\nSaved R/results/verify_sampler_port.rds\n")
