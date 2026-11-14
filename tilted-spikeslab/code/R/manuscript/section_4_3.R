# section_4_3.R
#
# Reproduces the prior-sampling study of Section 4.3 of the manuscript (table
# tab:joint-sparsening): the percentage of the specified expected edge count
# M p_inc, M = q(q-1)/2, that the joint prior distribution retains,
# 100 E[m] / (M p_inc), across coupling eta, exponent delta, inclusion
# probability p_inc, and dimension q.
#
# Design
#   Cells: q in {5, 10, 20, 50} x eta in {0.5, 1, 2} x delta in {0, 1, 2,
#   0.5 log q} x p_inc in {0.1, 0.25, 0.5}, and q = 100 at p_inc = 0.1 only
#   (156 cells). sigma = 1 and beta = eta throughout; Bernoulli graph weights.
#   Chains: 25 independent chains per cell of bgms::sample_ggm_prior(spec =
#   "joint", update_method = "gibbs"), the row-block Gibbs sweep and conjugate
#   edge move of the manuscript run without data, 500 warmup and 3500 retained
#   iterations, started from a graph drawn from the Bernoulli weights. Chain
#   seed: q * 1e5 + round(10 eta) * 1e3 + round(100 delta) + chain.
#   Priors: slab N(0, 1) on included off-diagonal entries and Exp(eta) on the
#   diagonal. bgms places its priors on K/2, so these are normal_prior(scale =
#   0.5) and exponential_prior(eta = eta) (the standardized rate is unchanged;
#   verified against enumeration at q = 5).
#   Estimates: per chain, the mean edge count over the retained iterations;
#   per cell, the mean over completed chains. Monte Carlo standard error:
#   standard deviation of the chain means divided by the square root of the
#   number of completed chains. Convergence: split R-hat of the edge count over
#   the chains of a cell. Failed chains are counted and reported with their
#   error messages; nothing is dropped.
#   Every chain is written to its own file as soon as it finishes, so an
#   interrupted run resumes where it stopped.
#
# Requirements
#   R >= 4.3 with the packages bgms (version 0.2.0.1: GitHub
#   Bayesian-Graphical-Modelling-Lab/bgms, branch develop, commit 2da4bf30,
#   remotes::install_github("Bayesian-Graphical-Modelling-Lab/bgms@2da4bf30"))
#   and posterior (CRAN). A project-local library R/lib/bgms_0.2.0.1 is used
#   when present. parallel (base R) for the chains.
#
# Run (from the project root)
#   Rscript R/manuscript/section_4_3.R              runs the missing chains on
#                                                   12 cores (about 10 hours in
#                                                   total, 7.5 of them for q = 100)
#                                                   and builds the table
#   Rscript R/manuscript/section_4_3.R table-only   builds the table from the
#                                                   chains already on disk; cells
#                                                   without their 25 chains are
#                                                   left blank and listed
# Output
#   tables/section_4_3.tex           the table of Section 4.3
#   R/results/section_4_3.rds        per-chain table, per-cell summary, design,
#                                    session information
#   R/results/table1_bgms/chains/    one file per chain (edge-count trajectory,
#                                    mean, batch-means MCSE, timing, error)

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
local_lib <- file.path(PROJECT_ROOT, "R", "lib", "bgms_0.2.0.1")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages({ library(bgms); library(parallel) })
if (packageVersion("bgms") != "0.2.0.1") {
  message("bgms ", packageVersion("bgms"), " found; the manuscript used 0.2.0.1 (commit 2da4bf30)")
}
stopifnot(requireNamespace("posterior", quietly = TRUE))
table_only <- "table-only" %in% commandArgs(trailingOnly = TRUE)

## ---- design ----
n_chains <- 25L
n_samples <- 3500L
n_warmup <- 500L
n_cores <- 12L
etas <- c(0.5, 1, 2)
p_set <- c(0.1, 0.25, 0.5)
qs <- c(5L, 10L, 20L, 50L, 100L)
chain_dir <- file.path(PROJECT_ROOT, "R", "results", "table1_bgms", "chains")
dir.create(chain_dir, showWarnings = FALSE, recursive = TRUE)

make_cells <- function(q) {
  ps <- if (q == 100L) 0.1 else p_set
  cells <- list()
  for (eta in etas) for (delta in c(0, 1, 2, 0.5 * log(q))) for (p in ps)
    cells[[length(cells) + 1L]] <- list(q = q, eta = eta, delta = delta, p_inc = p)
  cells
}
cell_tag <- function(cell) sprintf("%d_%g_%.4f_%g", cell$q, cell$eta, cell$delta, cell$p_inc)
chain_file <- function(cell, chain_id) file.path(chain_dir, sprintf("%s_%02d.rds", cell_tag(cell), chain_id))
chain_seed <- function(cell, chain_id) {
  as.integer(cell$q * 100000L + round(cell$eta * 10) * 1000L + round(cell$delta * 100) + chain_id)
}
batch_mcse <- function(x, n_batch = 50L) {
  b <- floor(length(x) / n_batch)
  means <- vapply(seq_len(n_batch), function(i) mean(x[((i - 1L) * b + 1L):(i * b)]), numeric(1))
  stats::sd(means) / sqrt(n_batch)
}

## ---- one chain, written to its own file ----
run_chain <- function(cell, chain_id) {
  f <- chain_file(cell, chain_id)
  if (file.exists(f)) return(invisible(f))
  seed <- chain_seed(cell, chain_id)
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
  base <- data.frame(q = cell$q, eta = cell$eta, delta = cell$delta, p_inc = cell$p_inc,
                     chain = chain_id, seed = seed)
  if (is.character(res)) {
    out <- list(summary = cbind(base, E_m = NA_real_, mcse = NA_real_, completed = FALSE,
                                secs = secs, error = res),
                m = rep(NA_integer_, n_samples))
  } else {
    m <- as.integer(rowSums(res$edge_indicators))
    out <- list(summary = cbind(base, E_m = mean(m), mcse = batch_mcse(m), completed = TRUE,
                                secs = secs, error = NA_character_),
                m = m)
  }
  saveRDS(out, f)
  invisible(f)
}

## ---- run the missing chains, one q block at a time ----
t_all <- Sys.time()
if (!table_only) {
  cat(sprintf("[%s] Section 4.3: bgms %s; %d chains per cell; %d + %d iterations; %d cores\n",
              format(Sys.time()), as.character(packageVersion("bgms")), n_chains, n_warmup, n_samples, n_cores))
  for (q in qs) {
    jobs <- list()
    for (cell in make_cells(q)) for (i in seq_len(n_chains)) {
      if (!file.exists(chain_file(cell, i))) jobs[[length(jobs) + 1L]] <- list(cell = cell, id = i)
    }
    if (!length(jobs)) next
    t0 <- Sys.time()
    cat(sprintf("[%s] q = %d: %d chains to run\n", format(Sys.time()), q, length(jobs)))
    res <- mclapply(jobs, function(j) tryCatch(run_chain(j$cell, j$id), error = function(e) conditionMessage(e)),
                    mc.cores = n_cores, mc.preschedule = FALSE)
    bad <- vapply(res, function(r) !is.character(r) || !file.exists(r), logical(1))
    if (any(bad)) cat(sprintf("  %d chains did not produce a file\n", sum(bad)))
    cat(sprintf("[%s] q = %d done in %.1f h\n", format(Sys.time()), q,
                as.numeric(difftime(Sys.time(), t0, units = "hours"))))
  }
}

## ---- assemble what is on disk ----
cells <- do.call(c, lapply(qs, make_cells))
chains <- list(); rhats <- numeric(length(cells))
for (k in seq_along(cells)) {
  cell <- cells[[k]]
  objs <- lapply(seq_len(n_chains), function(i) {
    f <- chain_file(cell, i)
    if (file.exists(f)) readRDS(f) else NULL
  })
  objs <- Filter(Negate(is.null), objs)
  if (length(objs)) chains[[length(chains) + 1L]] <- do.call(rbind, lapply(objs, `[[`, "summary"))
  done <- Filter(function(o) isTRUE(o$summary$completed), objs)
  rhats[k] <- if (length(done) >= 2) posterior::rhat(do.call(cbind, lapply(done, `[[`, "m"))) else NA_real_
}
chains <- do.call(rbind, chains); rownames(chains) <- NULL
summ <- do.call(rbind, lapply(seq_along(cells), function(k) {
  cell <- cells[[k]]
  d <- chains[chains$q == cell$q & chains$eta == cell$eta & abs(chains$delta - cell$delta) < 1e-8 &
                chains$p_inc == cell$p_inc, ]
  done <- d$completed
  M <- cell$q * (cell$q - 1) / 2; E_w <- M * cell$p_inc
  E_m <- if (any(done)) mean(d$E_m[done]) else NA_real_
  se <- if (sum(done) >= 2) stats::sd(d$E_m[done]) / sqrt(sum(done)) else NA_real_
  data.frame(q = cell$q, eta = cell$eta, delta = cell$delta, p_inc = cell$p_inc, E_w = E_w, E_m = E_m,
             pct = 100 * E_m / E_w, mcse_pct = 100 * se / E_w, n_run = nrow(d), n_completed = sum(done),
             n_failed = sum(!done), rhat = rhats[k], mean_chain_secs = if (nrow(d)) mean(d$secs) else NA_real_)
}))
summ <- summ[order(summ$q, summ$p_inc, summ$eta, summ$delta), ]; rownames(summ) <- NULL
summ$complete <- summ$n_run == n_chains

cat("\n== per-cell estimates: percent of the specified count, MCSE in points, chains run/completed, R-hat ==\n")
print(transform(summ, pct = round(pct, 2), mcse_pct = round(mcse_pct, 3), rhat = round(rhat, 4),
                mean_chain_secs = round(mean_chain_secs)), row.names = FALSE)
if (any(!summ$complete)) {
  cat("\n== cells without their 25 chains (left blank in the table) ==\n")
  print(summ[!summ$complete, c("q", "eta", "delta", "p_inc", "n_run", "pct")], row.names = FALSE)
}
if (any(summ$n_failed > 0)) {
  cat("\n== failed chains ==\n")
  print(chains[!chains$completed, c("q", "eta", "delta", "p_inc", "chain", "seed", "error")], row.names = FALSE)
}
full <- summ[summ$complete, ]
cat(sprintf("\nComplete cells: %d of %d; failed chains: %d; max R-hat %.4f; max MCSE %.3f points (q = 5), %.3f (q = 10), %.3f (q >= 20)\n",
            nrow(full), nrow(summ), sum(summ$n_failed), max(full$rhat, na.rm = TRUE),
            max(full$mcse_pct[full$q == 5]), max(full$mcse_pct[full$q == 10]), max(full$mcse_pct[full$q >= 20])))

## ---- table: rows eta x delta, column groups p_inc with q columns ----
fmt <- function(q, eta, dl, p) {
  if (q == 100L && p != 0.1) return("")
  r <- summ[summ$q == q & summ$eta == eta & abs(summ$delta - dl) < 1e-8 & summ$p_inc == p, ]
  if (!nrow(r) || !r$complete) return("--")
  sprintf("%.1f", r$pct)
}
delta_label <- c("0", "1", "2", "$\\tfrac12\\log q$")
delta_of <- function(j, q) c(0, 1, 2, 0.5 * log(q))[j]
q_cols <- list(`0.1` = qs, `0.25` = qs[qs < 100], `0.5` = qs[qs < 100])
header_groups <- paste(mapply(function(p, cols) sprintf("\\multicolumn{%d}{c}{$p_{\\mathrm{inc}}=%.2f$}", length(cols), p),
                              as.numeric(names(q_cols)), q_cols), collapse = " & ")
ends <- cumsum(lengths(q_cols)) + 2L; starts <- ends - lengths(q_cols) + 1L
rules <- paste(sprintf("\\cmidrule(lr){%d-%d}", starts, ends), collapse = "")
rows <- character(0)
for (eta in etas) {
  for (j in 1:4) {
    vals <- unlist(lapply(names(q_cols), function(p) vapply(q_cols[[p]], function(q) fmt(q, eta, delta_of(j, q), as.numeric(p)), character(1))))
    rows <- c(rows, sprintf("%s & %s & %s \\\\", if (j == 1) format(eta) else "", delta_label[j], paste(vals, collapse = " & ")))
  }
  if (eta != etas[length(etas)]) rows <- c(rows, "\\addlinespace")
}
blank_note <- if (any(!summ$complete)) sprintf("%% Blank entries (--): %d cells whose 25 chains had not all finished when the table was built.", sum(!summ$complete)) else "% All 156 cells complete."
tex <- c(
  sprintf("%% Produced by R/manuscript/section_4_3.R (bgms %s, %d chains x %d retained iterations per cell).",
          as.character(packageVersion("bgms")), n_chains, n_samples),
  sprintf("%% Entries: 100 E[m] / (M p_inc), one decimal. Monte Carlo standard errors at most %.2f points at q = 5, %.2f at q = 10, %.2f at q >= 20; max split R-hat %.3f; failed chains: %d.",
          max(full$mcse_pct[full$q == 5]), max(full$mcse_pct[full$q == 10]), max(full$mcse_pct[full$q >= 20]),
          max(full$rhat, na.rm = TRUE), sum(summ$n_failed)),
  blank_note,
  sprintf("\\begin{tabular}{ll%s}", strrep("r", sum(lengths(q_cols)))),
  "\\toprule",
  paste0(" & & ", header_groups, " \\\\"),
  rules,
  paste0("$\\eta$ & $\\delta$ & ", paste(unlist(lapply(q_cols, function(cols) sprintf("$%d$", cols))), collapse = " & "), " \\\\"),
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}")
writeLines(tex, file.path(PROJECT_ROOT, "tables", "section_4_3.tex"))
saveRDS(list(summary = summ, chains = chains,
             design = list(qs = qs, etas = etas, p_inc = p_set, deltas = "0, 1, 2, 0.5 log q", sigma = 1,
                           n_chains = n_chains, n_warmup = n_warmup, n_samples = n_samples,
                           chain_seed_rule = "q * 1e5 + round(10 eta) * 1e3 + round(100 delta) + chain",
                           bgms_mapping = "normal_prior(scale = 0.5), exponential_prior(eta = eta), bernoulli_prior(p_inc), spec = joint, update_method = gibbs"),
             bgms_version = as.character(packageVersion("bgms")), table_only = table_only,
             session = sessionInfo(), date = Sys.time()),
        file.path(PROJECT_ROOT, "R", "results", "section_4_3.rds"))
cat(sprintf("\n[%s] done in %.1f h; wrote tables/section_4_3.tex and R/results/section_4_3.rds\n",
            format(Sys.time()), as.numeric(difftime(Sys.time(), t_all, units = "hours"))))
