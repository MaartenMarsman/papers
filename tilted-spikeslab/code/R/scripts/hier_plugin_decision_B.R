# Section 5.4, component B: the plug-in chain against the exact-ratio chain
# at the posterior level (q = 4, the V6 data set, eta = 1, delta = 1,
# p_inc = 0.5).
#
# (i) Decision-discrepancy probe: along the exact-ratio chain's trajectory,
#     each candidate decision is also taken with the estimated J and the SAME
#     uniform; the discrepancy rate is the fraction of differing decisions.
# (ii) Posterior error: plug-in chains (spec = 2) at each N against the
#     exact-ratio chain (spec = 1), total variation over the 64 graphs and
#     maximum edge-marginal difference, two seeds each.
# (iii) Effective sample size per second of the edge-count series.
#
# N in {2, 10, 100, 1000, 10000}.
#
# Output: R/results/hier_plugin_decision_B.rds (+ log)

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})
PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
sourceCpp(file.path(PROJECT_ROOT, "R", "src", "spikeslab_sampler.cpp"))

q <- 4L
M <- 6L
n_obs <- 10
S <- readRDS(file.path(PROJECT_ROOT, "R", "results",
                       "verify_posterior_enumerated_V6.rds"))$S

zt <- readRDS(file.path(PROJECT_ROOT, "R", "results", "joint_spec_sparsening_J1_z.rds"))
d <- zt[zt$q == 4 & zt$eta == 1 & abs(zt$delta - 1) < 1e-9, ]
s1 <- d[d$seed == min(d$seed), ]; s2 <- d[d$seed == max(d$seed), ]
s1 <- s1[order(s1$mask), ]; s2 <- s2[order(s2$mask), ]
lz <- log((s1$z_hat + s2$z_hat) / 2)

tv <- function(p, phat) 0.5 * sum(abs(p - phat))
Ns <- c(2L, 10L, 100L, 1000L, 10000L)

ess_ims <- function(x) {
  n <- length(x)
  if (stats::var(x) == 0) return(NA_real_)
  rho <- as.numeric(stats::acf(x, lag.max = min(2000L, n - 2L),
                               plot = FALSE)$acf)
  npair <- floor(length(rho) / 2)
  G <- rho[2 * seq_len(npair) - 1] + rho[2 * seq_len(npair)]
  firstneg <- which(G <= 0)[1]
  if (!is.na(firstneg)) G <- G[seq_len(firstneg - 1)]
  if (length(G) == 0) return(n)
  G <- cummin(G)
  n / max(2 * sum(G) - 1, 1e-8)
}

run_chain <- function(spec, N, seed, n_iter) {
  set.seed(seed)
  t0 <- Sys.time()
  ch <- sampler_chain_cpp(S = S, n_obs = n_obs, q = q, delta = 1, beta = 1,
                          sigma = 1, p_inc = 0.5, spec = spec, log_z = lz,
                          n_iter = n_iter, n_burn = 2000L,
                          refresh_every = 50L, track_masks = TRUE,
                          n_aux = N, stat_thin = 1L)
  ch$secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ch
}

## ---- exact-ratio reference chains (spec 1, no probe) ----
cat("== exact-ratio reference (spec 1) ==\n")
ref <- mclapply(c(9601L, 9602L), function(seed)
  run_chain(1L, 0L, seed, 6e4L), mc.cores = 2L)
p_ref <- (ref[[1]]$mask_visits / sum(ref[[1]]$mask_visits) +
          ref[[2]]$mask_visits / sum(ref[[2]]$mask_visits)) / 2
ref_floor <- tv(ref[[1]]$mask_visits / sum(ref[[1]]$mask_visits),
                ref[[2]]$mask_visits / sum(ref[[2]]$mask_visits))
ess_ref <- mean(vapply(ref, function(ch) ess_ims(ch$rec_m), numeric(1)))
essps_ref <- ess_ref / mean(vapply(ref, function(ch) ch$secs, numeric(1)))
cat(sprintf("between-seed TV floor = %.4f; E[m] = %.3f; ESS(m)/s = %.0f\n",
            ref_floor, mean(c(ref[[1]]$E_m, ref[[2]]$E_m)), essps_ref))

## ---- (i) discrepancy probe + (ii) plug-in chains ----
jobs <- list()
for (N in Ns) {
  jobs[[length(jobs) + 1L]] <- list(kind = "probe", N = N,
                                    n_iter = if (N >= 10000L) 1e4L else 2e4L)
  for (seed in c(9701L, 9702L))
    jobs[[length(jobs) + 1L]] <- list(kind = "plugin", N = N, seed = seed,
                                      n_iter = if (N >= 10000L) 1e4L else 2e4L)
}
run_job <- function(job) {
  if (job$kind == "probe") {
    set.seed(9650L + job$N %% 1000L)
    ch <- sampler_chain_cpp(S = S, n_obs = n_obs, q = q, delta = 1, beta = 1,
                            sigma = 1, p_inc = 0.5, spec = 1L, log_z = lz,
                            n_iter = job$n_iter, n_burn = 2000L,
                            refresh_every = 50L, track_masks = FALSE,
                            n_aux = job$N)
    list(kind = "probe", N = job$N, disc_rate = ch$disc_rate)
  } else {
    ch <- run_chain(2L, job$N, job$seed, job$n_iter)
    list(kind = "plugin", N = job$N, seed = job$seed,
         phat = ch$mask_visits / sum(ch$mask_visits), E_m = ch$E_m,
         ess_ps = ess_ims(ch$rec_m) / ch$secs)
  }
}
out <- mclapply(jobs, run_job, mc.cores = 12L, mc.preschedule = FALSE)

rows <- list()
for (N in Ns) {
  pr <- Filter(function(o) o$kind == "probe" && o$N == N, out)[[1]]
  pl <- Filter(function(o) o$kind == "plugin" && o$N == N, out)
  pooled <- (pl[[1]]$phat + pl[[2]]$phat) / 2
  rows[[length(rows) + 1L]] <- data.frame(
    N = N, disc_rate = pr$disc_rate,
    tv_vs_exact = tv(p_ref, pooled),
    tv_between_seeds = tv(pl[[1]]$phat, pl[[2]]$phat),
    E_m = mean(c(pl[[1]]$E_m, pl[[2]]$E_m)),
    ess_per_sec = mean(c(pl[[1]]$ess_ps, pl[[2]]$ess_ps)))
}
res <- do.call(rbind, rows)
attr(res, "reference") <- list(ref_floor = ref_floor,
                               E_m_ref = mean(c(ref[[1]]$E_m, ref[[2]]$E_m)),
                               ess_per_sec_ref = essps_ref)
cat("\n== plug-in chain vs exact-ratio chain (posterior, V6 data) ==\n")
print(res, row.names = FALSE, digits = 3)
cat(sprintf("reference: TV floor %.4f, E[m] %.3f, ESS(m)/s %.0f\n",
            ref_floor, attr(res, "reference")$E_m_ref, essps_ref))
saveRDS(res, file.path(PROJECT_ROOT, "R", "results", "hier_plugin_decision_B.rds"))
cat("Saved R/results/hier_plugin_decision_B.rds\n")
