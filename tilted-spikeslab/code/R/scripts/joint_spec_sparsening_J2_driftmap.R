# Drift map for J2: maximum within-sweep Sigma drift (against the exact
# inverse recomputed at the end of every sweep) by cell class, to identify
# which J2 cells are numerically reliable.
#
# Within a sweep the toggle acceptances read entries of the maintained Sigma;
# large within-sweep drift therefore invalidates the kernel in that sweep even
# though the per-sweep refresh restores Sigma afterwards. Cells whose drift is
# not small must be flagged, not quoted as estimates.
#
# Output: R/results/joint_spec_sparsening_J2_driftmap.rds (+ log)

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})
PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
sourceCpp(file.path(PROJECT_ROOT, "R", "src", "spikeslab_sampler.cpp"))

cells <- list()
for (q in c(20L, 50L)) for (eta in c(0.5, 1, 2)) for (delta in c(0, 1, 2))
  for (p in c(0.1, 0.25, 0.5)) {
    cells[[length(cells) + 1L]] <- list(q = q, eta = eta, delta = delta,
                                        p_inc = p)
  }
cat(sprintf("running %d chains on 12 cores\n", length(cells)))

run_cell <- function(cell) {
  set.seed(4401L)
  ch <- tryCatch(
    sampler_chain_cpp(S = matrix(0, cell$q, cell$q), n_obs = 0, q = cell$q,
                      delta = cell$delta, beta = cell$eta, sigma = 1,
                      p_inc = cell$p_inc, spec = 0L, log_z = numeric(0),
                      n_iter = 4000L, n_burn = 500L, refresh_every = 1L,
                      track_masks = FALSE),
    error = function(e) NULL)
  data.frame(q = cell$q, eta = cell$eta, delta = cell$delta,
             p_inc = cell$p_inc,
             drift = if (is.null(ch)) NA_real_ else ch$max_sigma_drift,
             failed = is.null(ch))
}

res <- do.call(rbind, mclapply(cells, run_cell, mc.cores = 12L,
                               mc.preschedule = FALSE))
saveRDS(res, file.path(PROJECT_ROOT, "R", "results",
                       "joint_spec_sparsening_J2_driftmap.rds"))

res <- res[order(res$delta, res$q, res$p_inc, res$eta), ]
cat("\n== max within-sweep Sigma drift by cell ==\n")
print(transform(res, drift = sprintf("%.2e", drift)), row.names = FALSE)

cat("\n== summary by delta ==\n")
for (d in c(0, 1, 2)) {
  s <- res[res$delta == d & !res$failed, ]
  cat(sprintf("delta = %d: max drift %.2e (n fail = %d)\n",
              d, max(s$drift), sum(res$failed[res$delta == d])))
}
