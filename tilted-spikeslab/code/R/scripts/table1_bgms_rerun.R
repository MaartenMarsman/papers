# Table 1 (Section 4.3) prior-sampling study rerun in bgms 0.2.0.1 with 25
# chains per cell (decisions 2026-11-14; see notes/table1-bgms-rerun.md for
# the audit, the parametrization check, and the analysis rules).
#
# Design: q in {5, 10, 20, 50} x eta in {0.5, 1, 2} x delta in
# {0, 1, 2, 0.5 log q} x p_inc in {0.1, 0.25, 0.5}; q = 100 at p_inc = 0.1 only.
# 3500 retained iterations after 500 warmup, joint specification, Bernoulli
# graph weights, slab sd 1, beta = eta. bgms mapping (verified in
# bgms_param_check_q5.R): normal_prior(scale = 0.5), exponential_prior(eta = eta).
#
# Recording per chain: mean edge count, batch-means MCSE, completion flag and
# error message, run time, and the edge-count trajectory. Failed chains are
# counted and kept. Every chain is written to its own file as soon as it
# finishes, so an interrupted run resumes where it stopped.
#
# Run:    Rscript R/scripts/table1_bgms_rerun.R            (all q blocks)
#         Rscript R/scripts/table1_bgms_rerun.R 50 100     (selected q blocks)
# Output: R/results/table1_bgms/chains/<q>_<eta>_<delta>_<p>_<chain>.rds,
#         R/results/table1_bgms/cells.rds (per-chain summaries, all cells),
#         R/results/table1_bgms/summary.rds (per-cell estimates), + run.log.

suppressPackageStartupMessages(library(parallel))
PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
.libPaths(c(file.path(PROJECT_ROOT, "R", "lib", "bgms_0.2.0.1"), .libPaths()))
suppressPackageStartupMessages(library(bgms))
stopifnot(packageVersion("bgms") == "0.2.0.1")

## ---- parameters ----
n_chains <- 25L
n_samples <- 3500L
n_warmup <- 500L
n_cores <- 12L
out_dir <- file.path(PROJECT_ROOT, "R", "results", "table1_bgms")
chain_dir <- file.path(out_dir, "chains")
dir.create(chain_dir, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
qs <- if (length(args) > 0) as.integer(args) else c(5L, 10L, 20L, 50L, 100L)

make_cells <- function(q) {
  p_set <- if (q == 100L) 0.1 else c(0.1, 0.25, 0.5)
  cells <- list()
  for (eta in c(0.5, 1, 2)) for (delta in c(0, 1, 2, 0.5 * log(q))) for (p in p_set)
    cells[[length(cells) + 1L]] <- list(q = q, eta = eta, delta = delta, p_inc = p)
  cells
}
cell_tag <- function(cell) sprintf("%d_%g_%.4f_%g", cell$q, cell$eta, cell$delta, cell$p_inc)
chain_file <- function(cell, chain_id) {
  file.path(chain_dir, sprintf("%s_%02d.rds", cell_tag(cell), chain_id))
}

batch_mcse <- function(x, n_batch = 50L) {
  b <- floor(length(x) / n_batch)
  means <- vapply(seq_len(n_batch), function(i) mean(x[((i - 1L) * b + 1L):(i * b)]),
                  numeric(1))
  stats::sd(means) / sqrt(n_batch)
}

## one chain; the seed is distinct across the whole study
run_chain <- function(cell, chain_id) {
  f <- chain_file(cell, chain_id)
  if (file.exists(f)) return(invisible(f))
  seed <- as.integer(cell$q * 100000L + round(cell$eta * 10) * 1000L +
                       round(cell$delta * 100) + chain_id)
  t0 <- Sys.time()
  res <- tryCatch(
    sample_ggm_prior(p = cell$q, n_samples = n_samples, n_warmup = n_warmup,
                     interaction_prior = normal_prior(scale = 0.5),
                     precision_scale_prior = exponential_prior(eta = cell$eta),
                     delta = cell$delta, spec = "joint",
                     edge_prior = bernoulli_prior(cell$p_inc),
                     update_method = "gibbs", seed = seed, verbose = FALSE),
    error = function(e) conditionMessage(e))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.character(res)) {
    out <- list(summary = data.frame(q = cell$q, eta = cell$eta, delta = cell$delta,
                                     p_inc = cell$p_inc, chain = chain_id, seed = seed,
                                     E_m = NA_real_, mcse = NA_real_, completed = FALSE,
                                     secs = secs, error = res),
                m = rep(NA_integer_, n_samples))
  } else {
    m <- as.integer(rowSums(res$edge_indicators))
    out <- list(summary = data.frame(q = cell$q, eta = cell$eta, delta = cell$delta,
                                     p_inc = cell$p_inc, chain = chain_id, seed = seed,
                                     E_m = mean(m), mcse = batch_mcse(m), completed = TRUE,
                                     secs = secs, error = NA_character_),
                m = m)
  }
  saveRDS(out, f)
  invisible(f)
}

## ---- run, one q block at a time, all chains of the block in one pool ----
t_all <- Sys.time()
cat(sprintf("[%s] bgms %s; q blocks %s; %d chains per cell; %d cores\n",
            format(Sys.time()), as.character(packageVersion("bgms")),
            paste(qs, collapse = ","), n_chains, n_cores))
for (q in qs) {
  cells <- make_cells(q)
  jobs <- list()
  for (cell in cells) for (i in seq_len(n_chains)) jobs[[length(jobs) + 1L]] <- list(cell = cell, id = i)
  t0 <- Sys.time()
  cat(sprintf("[%s] q = %d: %d cells, %d chains\n", format(Sys.time()), q, length(cells), length(jobs)))
  res <- mclapply(jobs, function(j) tryCatch(run_chain(j$cell, j$id), error = function(e) conditionMessage(e)),
                  mc.cores = n_cores, mc.preschedule = FALSE)
  bad <- vapply(res, function(r) !is.character(r) || !file.exists(r), logical(1))
  if (any(bad)) cat(sprintf("  %d chains did not produce a file (worker crash?)\n", sum(bad)))
  cat(sprintf("[%s] q = %d done in %.1f h\n", format(Sys.time()), q,
              as.numeric(difftime(Sys.time(), t0, units = "hours"))))
}

## ---- assemble everything on disk ----
files <- list.files(chain_dir, pattern = "\\.rds$", full.names = TRUE)
chains <- do.call(rbind, lapply(files, function(f) readRDS(f)$summary))
chains <- chains[order(chains$q, chains$p_inc, chains$eta, chains$delta, chains$chain), ]
rownames(chains) <- NULL
saveRDS(chains, file.path(out_dir, "cells.rds"))

key <- interaction(chains$q, chains$eta, round(chains$delta, 4), chains$p_inc, drop = TRUE)
summ <- do.call(rbind, lapply(split(seq_len(nrow(chains)), key), function(idx) {
  d <- chains[idx, ]
  done <- d$completed
  M <- d$q[1] * (d$q[1] - 1) / 2
  E_w <- M * d$p_inc[1]
  E_m <- if (any(done)) mean(d$E_m[done]) else NA_real_
  mcse_cell <- if (sum(done) >= 2) stats::sd(d$E_m[done]) / sqrt(sum(done)) else NA_real_
  rhat <- NA_real_
  if (sum(done) >= 2 && requireNamespace("posterior", quietly = TRUE)) {
    cell <- list(q = d$q[1], eta = d$eta[1], delta = d$delta[1], p_inc = d$p_inc[1])
    m_mat <- do.call(cbind, lapply(d$chain[done], function(i) readRDS(chain_file(cell, i))$m))
    rhat <- posterior::rhat(m_mat)
  }
  data.frame(q = d$q[1], eta = d$eta[1], delta = d$delta[1], p_inc = d$p_inc[1],
             E_w = E_w, E_m = E_m, pct = 100 * E_m / E_w, mcse_pct = 100 * mcse_cell / E_w,
             n_chains = nrow(d), n_failed = sum(!done), rhat = rhat,
             mean_chain_secs = mean(d$secs))
}))
rownames(summ) <- NULL
summ <- summ[order(summ$q, summ$p_inc, summ$eta, summ$delta), ]
saveRDS(summ, file.path(out_dir, "summary.rds"))

cat("\n== per-cell estimates (percent of specified count, MCSE in points, failures, R-hat) ==\n")
print(transform(summ, pct = round(pct, 2), mcse_pct = round(mcse_pct, 3), rhat = round(rhat, 3),
                mean_chain_secs = round(mean_chain_secs)),
      row.names = FALSE)
if (any(summ$n_failed > 0)) {
  cat("\n== failure messages ==\n")
  print(table(chains$error[!chains$completed]))
}
cat(sprintf("\n[%s] all done in %.1f h\n", format(Sys.time()),
            as.numeric(difftime(Sys.time(), t_all, units = "hours"))))
