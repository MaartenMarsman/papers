# R/scripts/cone_interior_diagnostics_expanded.R
#
# Expanded version of the cone-interior diagnostic. Sweeps q in {5, 10, 15, 20,
# 30, 50} crossed with sigma in {0.3, 1, 3}, for both discrete and continuous
# spike-and-slab variants. Records lambda_min and log|K| for every draw so the
# downstream figure scripts can compute importance-weighted statistics and the
# empirical tail-bound for arbitrary delta.
#
# Sample sizes scale with q to give at least ~10K PD draws where the cone-
# interior issue is mild (sigma=0.3) and to push as deep as the brute-force MC
# allows where the issue is severe (sigma=1, sigma=3). At sigma >= 1 and q >=
# 30 the cone-exit rate is below the brute-force resolution; we keep those
# cells in the sweep and report 0 PD draws.

suppressPackageStartupMessages({
  library(Rcpp)
  library(RcppArmadillo)
  library(parallel)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
CPP_SRC      <- file.path(PROJECT_ROOT, "R", "src", "Z_prior_diag.cpp")
OUT_RDS      <- file.path(PROJECT_ROOT, "R", "results", "cone_interior_expanded.rds")

Rcpp::sourceCpp(CPP_SRC)

# ----- single-cell driver -----

run_one_cell <- function(args) {
  q        <- args$q
  variant  <- args$variant
  alpha    <- args$alpha
  beta     <- args$beta
  sigma    <- args$sigma
  v0       <- args$v0
  v1       <- args$v1
  gamma_p  <- args$gamma_prob
  n_draws  <- args$n_draws
  seed     <- args$seed
  delta_grid <- args$delta_grid

  variant_int <- if (variant == "discrete") 0L else 1L
  res <- z_prior_diag_cpp(q, alpha, beta, variant_int,
                          sigma, v0, v1, gamma_p, n_draws, seed)

  lmin      <- res$lmin
  log_det   <- res$log_det
  diag_mean <- res$diag_mean
  pd        <- res$pd == 1L
  pd_rate   <- mean(pd)
  lmin_pd   <- lmin[pd]
  log_det_pd <- log_det[pd]
  diag_pd   <- diag_mean[pd]
  n_pd      <- length(lmin_pd)

  # Marginal stats on PD draws (untilted)
  eps_grid <- c(0, 0.01, 0.025, 0.05, 0.1, 0.2, 0.5)
  pr_above <- vapply(eps_grid, function(e) mean(lmin > e), numeric(1L))

  weighted_q <- function(x, w, probs) {
    if (length(x) < 10L) return(setNames(rep(NA_real_, length(probs)), sprintf("q%g", probs)))
    ord <- order(x)
    x <- x[ord]; w <- w[ord]
    cw <- cumsum(w) / sum(w)
    out <- vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1L))
    setNames(out, sprintf("q%g", probs))
  }

  # Tilt summary: median, IQR, Pr[lmin <= eps] for each delta
  tilt_summary <- lapply(delta_grid, function(d) {
    if (n_pd < 10L) {
      return(list(delta = d, median = NA_real_, q05 = NA_real_, q95 = NA_real_,
                  pr_lmin_le_eps = setNames(rep(NA_real_, length(eps_grid)),
                                             sprintf("eps_%g", eps_grid)),
                  ess = 0))
    }
    log_w <- d * log_det_pd
    log_w <- log_w - max(log_w)
    w <- exp(log_w)
    ess <- sum(w)^2 / sum(w^2)
    qs <- weighted_q(lmin_pd, w, c(0.05, 0.50, 0.95))
    pr_le_w <- vapply(eps_grid, function(e) sum(w * (lmin_pd <= e)) / sum(w), numeric(1L))
    list(delta = d,
         q05    = qs[1L], median = qs[2L], q95 = qs[3L],
         pr_lmin_le_eps = setNames(pr_le_w, sprintf("eps_%g", eps_grid)),
         ess = ess)
  })

  # Down-sample PD-draw lmin / log_det vectors for downstream figures (cap at 20K)
  cap <- 20000L
  if (n_pd > cap) {
    idx <- sample.int(n_pd, cap)
    lmin_pd_keep    <- lmin_pd[idx]
    log_det_pd_keep <- log_det_pd[idx]
  } else {
    lmin_pd_keep    <- lmin_pd
    log_det_pd_keep <- log_det_pd
  }

  list(
    q = q, variant = variant,
    alpha = alpha, beta = beta, sigma = sigma, v0 = v0, v1 = v1,
    gamma_prob = gamma_p, n_draws = n_draws, n_pd = n_pd,
    pd_rate = pd_rate,
    eps_grid = eps_grid, pr_above = pr_above,
    lmin_pd_mean   = if (n_pd > 0L) mean(lmin_pd) else NA_real_,
    lmin_pd_median = if (n_pd > 0L) stats::median(lmin_pd) else NA_real_,
    lmin_pd_q05    = if (n_pd > 10L) stats::quantile(lmin_pd, 0.05, names = FALSE) else NA_real_,
    lmin_pd_q95    = if (n_pd > 10L) stats::quantile(lmin_pd, 0.95, names = FALSE) else NA_real_,
    diag_mean_pd   = if (n_pd > 0L) mean(diag_pd) else NA_real_,
    diag_mean_all  = mean(diag_mean),
    tilt_summary   = tilt_summary,
    # Per-draw vectors (capped) for downstream tail-bound figure
    lmin_pd_samples    = lmin_pd_keep,
    log_det_pd_samples = log_det_pd_keep
  )
}

# ----- sweep design -----

ALPHA       <- 3.0
BETA        <- 3.0
V0          <- 0.1   # only used for continuous; spike SD
V1          <- 1.0   # only used for continuous; slab SD when sigma_eff is 1
DELTA_GRID  <- c(0, 0.5, 1, 2)

# q-dependent edge probability (sparse-asymptotic regime), matching the Z project's setup
p_inc_for_q <- c(`5` = 0.50, `10` = 0.20, `15` = 0.13, `20` = 0.10,
                 `30` = 0.06, `50` = 0.03)

# Sample-size table: n_draws[sigma_label, q_label]
# Chosen so cells with measurable PD rate get >= 5K PD samples; cells likely
# zero (sigma=1, q>=30; sigma=3, q>=10) get 5e6 to give a clean upper bound.
n_draws_by_cell <- list(
  `0.3` = c(`5` = 2e5L,  `10` = 5e5L,  `15` = 1e6L,  `20` = 2e6L,    `30` = 5e6L, `50` = 1e7L),
  `1`   = c(`5` = 5e5L,  `10` = 2e6L,  `15` = 5e6L,  `20` = 1.5e7L,  `30` = 1e7L, `50` = 1e7L),
  `3`   = c(`5` = 1e6L,  `10` = 1e6L,  `15` = 1e6L,  `20` = 1e6L,    `30` = 1e6L, `50` = 1e6L)
)

q_grid     <- c(5L, 10L, 15L, 20L, 30L, 50L)
sigma_grid <- c(0.3, 1.0, 3.0)
variants   <- c("discrete", "continuous")

jobs <- list()
k <- 0L
for (sigma in sigma_grid) {
  for (q in q_grid) {
    for (variant in variants) {
      k <- k + 1L
      n_draws <- n_draws_by_cell[[as.character(sigma)]][as.character(q)]
      jobs[[length(jobs) + 1L]] <- list(
        q          = q,
        variant    = variant,
        alpha      = ALPHA,
        beta       = BETA,
        sigma      = sigma,
        v0         = V0,
        v1         = sigma,   # continuous slab SD set to sigma so the two variants are comparable
        gamma_prob = p_inc_for_q[as.character(q)],
        n_draws    = n_draws,
        delta_grid = DELTA_GRID,
        seed       = 42L + 1000L * k
      )
    }
  }
}

cat(sprintf("Running %d cells on 12 cores...\n", length(jobs)))
cat(sprintf("Total draws: %g\n", sum(vapply(jobs, function(j) j$n_draws, numeric(1L)))))

t0 <- Sys.time()
results <- parallel::mclapply(jobs, run_one_cell, mc.cores = 12L)
t1 <- Sys.time()
cat(sprintf("Done in %.1f s.\n", as.numeric(t1 - t0, units = "secs")))

# Summary tables
tab <- data.frame(
  variant = vapply(results, function(r) r$variant, character(1L)),
  q       = vapply(results, function(r) r$q,       integer(1L)),
  sigma   = vapply(results, function(r) r$sigma,   numeric(1L)),
  n_draws = vapply(results, function(r) r$n_draws, integer(1L)),
  n_pd    = vapply(results, function(r) r$n_pd,    integer(1L)),
  pd_rate = vapply(results, function(r) r$pd_rate, numeric(1L)),
  lmin_pd_median = vapply(results, function(r) r$lmin_pd_median, numeric(1L)),
  diag_mean_pd   = vapply(results, function(r) r$diag_mean_pd,   numeric(1L)),
  stringsAsFactors = FALSE
)
tab <- tab[order(tab$variant, tab$sigma, tab$q), ]
rownames(tab) <- NULL

cat("\n=== PD rate table ===\n")
print(tab[, c("variant", "sigma", "q", "n_draws", "n_pd", "pd_rate", "lmin_pd_median", "diag_mean_pd")],
      row.names = FALSE, digits = 4)

saveRDS(list(jobs = jobs, results = results, summary = tab,
             alpha = ALPHA, beta = BETA, delta_grid = DELTA_GRID,
             p_inc_for_q = p_inc_for_q),
        OUT_RDS)
cat(sprintf("\nSaved expanded diagnostics to %s\n", OUT_RDS))
