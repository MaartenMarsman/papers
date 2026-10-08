# R/scripts/standardized_rate_geometry.R
# Purpose: Characterise how the STANDARDISED diagonal rate beta governs the geometry
#          of the unit-scale base matrix K (off-diagonals N(0,1), diagonal
#          Gamma(alpha, beta)), and how the determinant tilt |K|^delta acts on the
#          lower tail of lambda_min(K).
#   (1) lambda_min crossover: concentrated (large alpha) vs dispersed (small alpha)
#       regime, against the Weyl/Wigner prediction alpha/beta - 2 sqrt(q).
#   (2) beta as a pure 1/beta scale on lambda_min with an (alpha, q)-dependent prefactor.
#   (3) tilt effect: importance sampling p(K) propto |K|^delta * 1{K > 0} shows the
#       tilt LIFTS and CONCENTRATES the lower tail of lambda_min and shrinks the
#       boundary mass Pr[lambda_min < eps]. Effective sample size reported.
# Run:     Rscript R/scripts/standardized_rate_geometry.R
# Output:  R/results/standardized_rate_geometry.rds
# Notes:   base R only; standardised (unit off-diagonal) frame; discrete-variant base.
# See:     notes/standardized-rate-geometry.md, notes/scale-injection-congruence.md

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
set.seed(20260616)

lam_all <- function(M) eigen(M, symmetric = TRUE, only.values = TRUE)$values

## draw one base matrix K = diag(Gamma(alpha,beta)) + Wigner(N(0,1) off-diag)
draw_K <- function(q, alpha, beta) {
  d <- stats::rgamma(q, shape = alpha, rate = beta)
  R <- matrix(0, q, q)
  R[upper.tri(R)] <- stats::rnorm(q * (q - 1) / 2)
  R <- R + t(R)
  diag(R) <- d
  R
}

## =====================================================================
## (1) + (2) lambda_min: concentration crossover and 1/beta scaling
## =====================================================================
summary_lmin <- function(q, alpha, beta, nrep = 3000) {
  lmin <- numeric(nrep); mind <- numeric(nrep)
  Rnorm <- numeric(nrep); sandwich_ok <- 0L; pr_pd <- 0
  for (r in seq_len(nrep)) {
    K <- draw_K(q, alpha, beta)
    ev <- lam_all(K)
    lmin[r] <- min(ev)
    md <- min(diag(K)); mind[r] <- md
    Roff <- K; diag(Roff) <- 0
    rn <- max(abs(lam_all(Roff)))               # operator norm of the off-diagonal part
    Rnorm[r] <- rn
    # deterministic interlacing/Weyl sandwich: md - ||R|| <= lambda_min <= md
    if (lmin[r] <= md + 1e-9 && lmin[r] >= md - rn - 1e-9) sandwich_ok <- sandwich_ok + 1L
    if (min(ev) > 0) pr_pd <- pr_pd + 1
  }
  data.frame(
    q = q, alpha = alpha, beta = beta,
    mean_dia  = alpha / beta,
    cv_dia    = 1 / sqrt(alpha),                 # diagonal CV, independent of beta
    pred_conc = alpha / beta - 2 * sqrt(q),      # concentrated-regime (large alpha)
    pred_disp = mean(mind) - 2 * sqrt(q),        # dispersed-regime (smallest diagonal)
    emp_lmin  = mean(lmin),
    sd_lmin   = stats::sd(lmin),
    drag      = mean(mind - lmin),               # off-diagonal pull below smallest diagonal
    mean_Rnorm = mean(Rnorm),                    # E||R||, ~ 2 sqrt(q) for full slab
    sandwich_frac = sandwich_ok / nrep,          # should be 1 (deterministic theorem)
    pr_PD     = pr_pd / nrep
  )
}

# crossover: hold mean diagonal alpha/beta = 20 at q = 20, vary alpha (concentration)
cfg1 <- data.frame(alpha = c(2, 5, 15, 50, 150))
cfg1$beta <- cfg1$alpha / 20
tab_crossover <- do.call(rbind, Map(function(a, b) summary_lmin(20, a, b),
                                    cfg1$alpha, cfg1$beta))

# 1/beta scaling: fix (q, alpha), vary beta; check emp_lmin * beta is ~constant
cfg2 <- data.frame(beta = c(0.1, 0.2, 0.4, 0.8))
tab_scale <- do.call(rbind, lapply(cfg2$beta, function(b) summary_lmin(20, 3, b)))
tab_scale$lmin_times_beta <- tab_scale$emp_lmin * tab_scale$beta

## =====================================================================
## (3) tilt |K|^delta: importance sampling on the cone-only base
##     weight w = |K|^delta * 1{K > 0}; weighted lower-tail statistics
## =====================================================================
tilt_lower_tail <- function(q, alpha, beta, deltas, eps = 0.05, N = 3e5) {
  # one shared pool of base draws; reuse across deltas (common random numbers)
  lmin <- numeric(N); logdet <- numeric(N); pd <- logical(N)
  for (i in seq_len(N)) {
    K <- draw_K(q, alpha, beta)
    ev <- lam_all(K)
    pd[i] <- min(ev) > 0
    lmin[i] <- min(ev)
    logdet[i] <- if (pd[i]) sum(log(ev)) else NA_real_
  }
  res <- lapply(deltas, function(delta) {
    lw <- ifelse(pd, delta * logdet, -Inf)        # log weight = delta * log|K| on the cone
    lw <- lw - max(lw[is.finite(lw)])
    w  <- exp(lw); w[!is.finite(w)] <- 0
    sw <- sum(w); ess <- sw^2 / sum(w^2)
    wm <- function(x) sum(w * x) / sw
    m  <- wm(lmin)
    s  <- sqrt(wm((lmin - m)^2))
    # weighted lower quantiles of lambda_min via weighted ECDF
    ord <- order(lmin); cw <- cumsum(w[ord]) / sw
    q05 <- lmin[ord][which(cw >= 0.05)[1]]
    q25 <- lmin[ord][which(cw >= 0.25)[1]]
    data.frame(q = q, alpha = alpha, beta = beta, delta = delta,
               ess = ess, w_mean_lmin = m, w_sd_lmin = s,
               w_q05_lmin = q05, w_q25_lmin = q25,
               w_pr_below_eps = wm(as.numeric(lmin < eps)))
  })
  do.call(rbind, res)
}

# choose beta near the boundary so the tilt has room to act (Pr[PD] moderate)
q3 <- 8; alpha3 <- 2; beta3 <- 0.45
tab_tilt <- tilt_lower_tail(q3, alpha3, beta3,
                            deltas = c(0, 0.5, 1, 2, 4), eps = 0.05, N = 3e5)

## =====================================================================
## (4) partial-correlation prior scale: SD(rho_ij) = eta/(alpha-1)
##     rho_ij = -k_ij / sqrt(k_ii k_jj); off-diagonal N(0,1) indep of diagonal.
##     Var(rho) = E[k_ij^2] E[1/k_ii] E[1/k_jj] = (eta/(alpha-1))^2 for alpha > 1.
##     Reported on cone (PD) draws; alpha = 1 has infinite inverse moment (heavy tail).
## =====================================================================
pcor_scale <- function(q, alpha, eta, nrep = 5000) {
  vals <- numeric(0); npd <- 0
  for (r in seq_len(nrep)) {
    K <- draw_K(q, alpha, eta)
    if (min(lam_all(K)) <= 0) next
    npd <- npd + 1
    P <- -K / sqrt(outer(diag(K), diag(K)))
    vals <- c(vals, P[upper.tri(P)])
  }
  data.frame(q = q, alpha = alpha, eta = eta,
             pred_uncond = if (alpha > 1) eta / (alpha - 1) else Inf,
             emp_sd_PD   = stats::sd(vals),
             pr_PD       = npd / nrep)
}
cfg4 <- expand.grid(alpha = c(1.5, 2, 3, 5), eta = c(0.1, 0.5))
tab_pcor <- do.call(rbind, Map(function(a, e) pcor_scale(8, a, e), cfg4$alpha, cfg4$eta))

## ---- report ----
cat("=== (1) lambda_min concentration crossover (mean diagonal fixed = 20, q = 20) ===\n")
print(round(tab_crossover, 3), row.names = FALSE)
cat("\n=== (2) 1/beta scaling of lambda_min (q = 20, alpha = 3) ===\n")
print(round(tab_scale, 3), row.names = FALSE)
cat("\n=== (3) tilt lifts lambda_min + evacuates the boundary (Pr[<eps] down)",
    sprintf("(q = %d, alpha = %g, eta = %g) ===\n", q3, alpha3, beta3))
print(round(tab_tilt, 4), row.names = FALSE)
cat("\n=== (4) partial-correlation SD(rho) ~ eta/(alpha-1) in the diagonally-dominant",
    "regime (q = 8); PD-truncated as eta/(alpha-1) grows ===\n")
print(round(tab_pcor, 4), row.names = FALSE)

## ---- save ----
results <- list(crossover = tab_crossover, scaling = tab_scale, tilt = tab_tilt,
                pcor = tab_pcor,
                meta = list(seed = 20260616, q_tilt = q3, alpha_tilt = alpha3,
                            beta_tilt = beta3))
dir.create("R/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(results, file.path(PROJECT_ROOT, "R", "results", "standardized_rate_geometry.rds"))
cat("\nDone. Results saved to R/results/standardized_rate_geometry.rds\n")
