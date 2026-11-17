# section_4_1.R
#
# Reproduces the fixed-graph calibration study of Section 4.1 of the manuscript
# (figure fig:delta-scaling): the tilt exponent at which the graph-averaged
# median smallest eigenvalue of K reaches one half, at unit slab scale with
# exponential diagonal rate beta = eta, across dimension, coupling, and two
# Erdos-Renyi graph regimes.
#
# Design
#   Dimensions q in {10, 20, 50}; coupling eta in {0.5, 1, 2} (sigma = 1,
#   beta = eta); graph regimes: fixed expected degree dbar in {2, 4, 8} with
#   edge probability dbar/(q-1), and fixed edge probability rho in {0.05, 0.10,
#   0.20}. Exponent grid delta in {0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24}.
#   Graphs: for each (q, regime), independent Erdos-Renyi graphs (no edge
#   forced in), seed 410000 + 100 q + regime index; 25 graphs per cell (the
#   first 25 of a seeded sequence of 200, so a pilot graph is the same graph in
#   the full study; the count was fixed after the pilot to keep the run within
#   one day on 12 cores, the q = 50 cells costing about 0.8 h per graph). The same graphs are used
#   at every eta and every delta, so differences across exponents are paired
#   within graphs.
#   Chains: for each graph and delta, four independent chains of
#   bgms::sample_ggm_prior(spec = "conditional"), a no-U-turn sampler on a
#   Cholesky parametrization of K that enforces the zero pattern of the graph,
#   1000 warmup and 2000 retained iterations each (8000 retained draws per
#   graph and delta). Chain seed: see chain_seed() below; distinct across the
#   whole study. Priors: slab N(0, 1) on included off-diagonal entries and
#   Exp(eta) on the diagonal; bgms places its priors on K/2, so these are
#   normal_prior(scale = 0.5) and exponential_prior(eta = eta) (the same
#   mapping as in section_3_6.R, verified there against an importance-sampling
#   reference); the output is on the K scale.
#   Per draw: the smallest eigenvalue of K (saved for every chain), and over the
#   included edges the squared entries and partial correlations.
#   Target: for each graph and delta, the median of the 8000 pooled draws of
#   the smallest eigenvalue; for each cell (q, regime level, eta) and delta,
#   the mean of these medians over the 200 graphs; the required exponent
#   delta* is where the mean curve reaches 0.5, by linear interpolation
#   between grid points (the previous study's rule). A monotone spline through
#   the grid gives a second value, reported as an interpolation check. Cells
#   whose mean curve does not reach the target on the grid are reported as
#   unbracketed, not extrapolated.
#   Uncertainty: bootstrap over graphs (2000 resamples of the 200 graphs with
#   replacement, each graph carrying its full curve over delta, so the pairing
#   across exponents is preserved); the mean curve and delta* are recomputed
#   for every resample and the 2.5% and 97.5% quantiles are reported. This
#   interval covers graph-to-graph variation within the regime and the Monte
#   Carlo error of each graph's median (which enters through the resampled
#   per-graph estimates); it does not cover a different graph regime or a
#   different target.
#   Monte Carlo precision of the medians: per graph and delta, the standard
#   deviation of the four chain medians divided by two (the Monte Carlo
#   standard error of the mean of the chain medians); reported relative to the
#   spread of the medians across graphs. Convergence: split R-hat of the
#   smallest eigenvalue over the four chains; bulk effective sample size.
#   bgms exports no divergence or tree-depth information for this sampler, so
#   those NUTS diagnostics are unavailable. Failed chains are recorded with
#   their error messages; nothing is dropped.
#
# Requirements
#   R >= 4.3 with bgms 0.2.0.1 (GitHub Bayesian-Graphical-Modelling-Lab/bgms,
#   branch develop, commit 2da4bf30; remotes::install_github(
#   "Bayesian-Graphical-Modelling-Lab/bgms@2da4bf30")), posterior, ggplot2,
#   patchwork (figure, through R/scripts/fig_style.R), parallel (base R). A
#   project-local library R/lib/bgms_0.2.0.1 is used when present.
#
# Run (from the project root)
#   Rscript R/manuscript/section_4_1.R pilot     the first 8 graphs of eight
#                                                representative cells (see
#                                                pilot_cells below) over the full
#                                                grid; these checkpoints are part
#                                                of the full study and are reused
#   Rscript R/manuscript/section_4_1.R           full study; resumes from the
#                                                per-graph checkpoints; then
#                                                summary and figure
#   Rscript R/manuscript/section_4_1.R summary   summary and figure from the
#                                                checkpoints on disk only
#   Environment: SECTION_4_1_CORES (default 12) sets the number of workers.
# Output
#   R/results/section_4_1/graphs/<q>_<axis><level>_eta<eta>_g<graph>.rds
#       (the pilot writes the same files)
#       per graph: adjacency matrix, seeds, smallest-eigenvalue draws
#       (2000 x 4 chains x 12 exponents), per-chain summaries and errors
#   R/results/section_4_1.rds        per-graph and per-cell summaries, delta*
#                                    with bootstrap intervals, diagnostics,
#                                    design, session information
#   figures/section_4_1.pdf          the calibration figure of Section 4.1
#   R/results/section_4_1.log        console output when run with
#                                    "> R/results/section_4_1.log 2>&1"

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
local_lib <- file.path(PROJECT_ROOT, "R", "lib", "bgms_0.2.0.1")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages({ library(bgms); library(parallel) })
if (packageVersion("bgms") != "0.2.0.1") {
  message("bgms ", packageVersion("bgms"), " found; the manuscript used 0.2.0.1 (commit 2da4bf30)")
}
stopifnot(requireNamespace("posterior", quietly = TRUE))
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "run"
stopifnot(mode %in% c("run", "pilot", "summary"))
n_cores <- as.integer(Sys.getenv("SECTION_4_1_CORES", unset = "12"))

## ---- design ----
sigma <- 1
qs <- c(10L, 20L, 50L)
regimes <- rbind(data.frame(axis = "degree", level = c(2, 4, 8)),
                 data.frame(axis = "density", level = c(0.05, 0.10, 0.20)))
etas <- c(0.5, 1, 2)
deltas <- c(0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24)
n_graphs_q <- c(`10` = 25L, `20` = 25L, `50` = 25L)   # graphs per cell (user decision 2026-11-16: one day of compute on 12 cores)
n_chains <- 4L
n_warmup <- 1000L
n_samples <- 2000L
target <- 0.5
n_boot <- 2000L
batch_graphs <- 12L     # graphs per checkpoint batch
## pilot: representative cells, including the demanding dense q = 50 ones; 8 graphs each
pilot_cells <- data.frame(q = c(50, 50, 50, 50, 20, 20, 10, 10),
                          axis = c("degree", "density", "degree", "density", "degree", "density", "degree", "density"),
                          level = c(2, 0.20, 8, 0.05, 2, 0.20, 2, 0.20), eta = c(1, 1, 2, 0.5, 1, 2, 1, 2))
n_graphs_pilot <- 8L
out_dir <- file.path(PROJECT_ROOT, "R", "results", "section_4_1")
ckpt_dir <- file.path(out_dir, "graphs")
dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)

edge_prob <- function(q, axis, level) if (axis == "degree") min(level / (q - 1), 1) else level
graph_seed <- function(q, r) 410000L + 100L * q + r           # r: regime row index (1..6)
chain_seed <- function(qi, ri, ei, di, g, c) {
  as.integer((((qi - 1L) * 6L + (ri - 1L)) * 3L + (ei - 1L)) * 10000000L + di * 1000000L + g * 10L + c)
}
cell_tag <- function(q, axis, level, eta) sprintf("%d_%s%g_eta%g", q, axis, level, eta)
graph_file <- function(q, axis, level, eta, g) file.path(ckpt_dir, sprintf("%s_g%03d.rds", cell_tag(q, axis, level, eta), g))

draw_graphs <- function(q, p, n, seed) {
  set.seed(seed)
  lapply(seq_len(n), function(i) {
    E <- matrix(0L, q, q)
    E[upper.tri(E)] <- as.integer(stats::runif(q * (q - 1) / 2) < p)
    E <- E + t(E); diag(E) <- 1L
    E
  })
}
## indices of the included edges in bgms's row-major upper-triangle column order
edge_index <- function(E) {
  q <- nrow(E)
  ii <- rep(seq_len(q - 1L), times = (q - 1L):1L)
  jj <- unlist(lapply(seq_len(q - 1L), function(i) (i + 1L):q))
  inc <- E[cbind(ii, jj)] == 1L
  list(col = which(inc), i = ii[inc], j = jj[inc], ut = cbind(ii, jj))
}

## ---- one chain: smallest eigenvalue per draw, edge summaries ----
run_chain <- function(q, E, eta, delta, seed, idx) {
  t0 <- Sys.time()
  res <- tryCatch(
    sample_ggm_prior(p = q, n_samples = n_samples, n_warmup = n_warmup,
                     interaction_prior = normal_prior(scale = sigma / 2),
                     precision_scale_prior = exponential_prior(eta = eta * sigma),
                     delta = delta, spec = "conditional", edge_indicators = E,
                     seed = seed, verbose = FALSE),
    error = function(e) conditionMessage(e))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.character(res)) return(list(ok = FALSE, error = res, secs = secs))
  D <- res$K_diag; O <- res$K_offdiag
  lmin <- vapply(seq_len(nrow(D)), function(t) {
    K <- diag(D[t, ]); K[idx$ut] <- O[t, ]; K <- K + t(K) - diag(D[t, ])
    min(eigen(K, symmetric = TRUE, only.values = TRUE)$values)
  }, numeric(1))
  if (length(idx$col)) {
    ke <- O[, idx$col, drop = FALSE]
    rho <- -ke / sqrt(D[, idx$i, drop = FALSE] * D[, idx$j, drop = FALSE])
    var_ratio <- mean(ke^2) / sigma^2; rho_sd <- stats::sd(rho)
  } else { var_ratio <- NA_real_; rho_sd <- NA_real_ }
  list(ok = TRUE, error = NA_character_, secs = secs, lmin = lmin, var_ratio = var_ratio, rho_sd = rho_sd)
}

## ---- one job: all chains of one graph at one exponent ----
run_job <- function(job) {
  idx <- edge_index(job$E)
  chains <- vector("list", n_chains); lmin <- matrix(NA_real_, n_samples, n_chains)
  for (c in seq_len(n_chains)) {
    seed <- chain_seed(job$qi, job$ri, job$ei, job$di, job$g, c)
    r <- run_chain(job$q, job$E, job$eta, job$delta, seed, idx)
    if (r$ok) lmin[, c] <- r$lmin
    chains[[c]] <- data.frame(delta = job$delta, chain = c, seed = seed, completed = r$ok, error = r$error,
                              secs = r$secs, median_lmin = if (r$ok) stats::median(r$lmin) else NA_real_,
                              ess_lmin = if (r$ok) posterior::ess_bulk(r$lmin) else NA_real_,
                              var_ratio = if (r$ok) r$var_ratio else NA_real_,
                              rho_sd = if (r$ok) r$rho_sd else NA_real_)
  }
  chains <- do.call(rbind, chains)
  done <- chains$completed
  list(delta = job$delta, lmin = lmin, chains = chains,
       rhat_lmin = if (sum(done) >= 2) posterior::rhat(lmin[, done, drop = FALSE]) else NA_real_,
       median_pooled = if (any(done)) stats::median(lmin[, done]) else NA_real_)
}

## ---- run: cells, batches of graphs, jobs = (graph, delta) ----
t_all <- Sys.time()
cells <- expand.grid(qi = seq_along(qs), ri = seq_len(nrow(regimes)), ei = seq_along(etas))
cells$n_graphs <- n_graphs_q[as.character(qs[cells$qi])]
if (mode == "pilot") {
  keep <- mapply(function(qi, ri, ei) any(pilot_cells$q == qs[qi] & pilot_cells$axis == regimes$axis[ri] &
                                            pilot_cells$level == regimes$level[ri] & pilot_cells$eta == etas[ei]),
                 cells$qi, cells$ri, cells$ei)
  cells <- cells[keep, ]; cells$n_graphs <- n_graphs_pilot
}
if (mode != "summary") {
  cat(sprintf("[%s] Section 4.1 (%s): bgms %s; %d cells; graphs per cell %s; %d exponents x %d chains; %d + %d iterations; %d cores\n",
              format(Sys.time()), mode, as.character(packageVersion("bgms")), nrow(cells),
              paste(unique(cells$n_graphs), collapse = "/"), length(deltas), n_chains, n_warmup, n_samples, n_cores))
  for (k in seq_len(nrow(cells))) {
    qi <- cells$qi[k]; ri <- cells$ri[k]; ei <- cells$ei[k]; n_graphs <- cells$n_graphs[k]
    q <- qs[qi]; axis <- regimes$axis[ri]; level <- regimes$level[ri]; eta <- etas[ei]
    graphs <- draw_graphs(q, edge_prob(q, axis, level), 200L, seed = graph_seed(q, ri))[seq_len(n_graphs)]
    todo <- which(!file.exists(vapply(seq_len(n_graphs), function(g) graph_file(q, axis, level, eta, g), character(1))))
    if (!length(todo)) next
    t0 <- Sys.time()
    cat(sprintf("[%s] cell %s: %d graphs to run\n", format(Sys.time()), cell_tag(q, axis, level, eta), length(todo)))
    for (b in split(todo, ceiling(seq_along(todo) / batch_graphs))) {
      jobs <- list()
      for (g in b) for (di in seq_along(deltas))
        jobs[[length(jobs) + 1L]] <- list(qi = qi, ri = ri, ei = ei, di = di, g = g, q = q, eta = eta,
                                          delta = deltas[di], E = graphs[[g]])
      res <- mclapply(jobs, function(j) tryCatch(run_job(j), error = function(e) conditionMessage(e)),
                      mc.cores = n_cores, mc.preschedule = FALSE)
      for (g in b) {
        sel <- which(vapply(jobs, function(j) j$g == g, logical(1)))
        rj <- res[sel]
        crashed <- vapply(rj, is.character, logical(1))
        if (any(crashed)) {
          cat(sprintf("  graph %d: %d exponent jobs crashed (%s); not checkpointed\n", g, sum(crashed), rj[crashed][[1]]))
          next
        }
        lmin <- array(NA_real_, c(n_samples, n_chains, length(deltas)), dimnames = list(NULL, NULL, deltas))
        for (m in seq_along(rj)) lmin[, , m] <- rj[[m]]$lmin
        chains <- do.call(rbind, lapply(rj, `[[`, "chains")); chains$graph <- g
        saveRDS(list(q = q, axis = axis, level = level, eta = eta, graph = g, E = graphs[[g]],
                     n_edges = (sum(graphs[[g]]) - q) / 2, lmin = lmin, chains = chains,
                     per_delta = data.frame(delta = deltas, rhat_lmin = vapply(rj, `[[`, numeric(1), "rhat_lmin"),
                                            median_pooled = vapply(rj, `[[`, numeric(1), "median_pooled"))),
                graph_file(q, axis, level, eta, g))
      }
    }
    cat(sprintf("[%s] cell %s done in %.2f h\n", format(Sys.time()), cell_tag(q, axis, level, eta),
                as.numeric(difftime(Sys.time(), t0, units = "hours"))))
  }
}

## ---- summary: per graph, per cell, inversion, bootstrap ----
files <- list.files(ckpt_dir, pattern = "\\.rds$", full.names = TRUE)
cat(sprintf("\n[%s] %d graph checkpoints on disk\n", format(Sys.time()), length(files)))
per_graph <- list(); chains_all <- list()
for (f in files) {
  o <- readRDS(f)
  ch <- o$chains
  md <- tapply(ch$median_lmin, ch$delta, function(x) stats::sd(x) / sqrt(sum(!is.na(x))))
  vr <- tapply(ch$var_ratio, ch$delta, mean, na.rm = TRUE); rs <- tapply(ch$rho_sd, ch$delta, mean, na.rm = TRUE)
  ess <- tapply(ch$ess_lmin, ch$delta, min); nf <- tapply(!ch$completed, ch$delta, sum)
  per_graph[[length(per_graph) + 1L]] <- data.frame(
    q = o$q, axis = o$axis, level = o$level, eta = o$eta, graph = o$graph, n_edges = o$n_edges,
    delta = o$per_delta$delta, median_lmin = o$per_delta$median_pooled, mcse_median = as.numeric(md[as.character(o$per_delta$delta)]),
    rhat_lmin = o$per_delta$rhat_lmin, min_ess_lmin = as.numeric(ess[as.character(o$per_delta$delta)]),
    n_failed = as.integer(nf[as.character(o$per_delta$delta)]),
    var_ratio = as.numeric(vr[as.character(o$per_delta$delta)]), rho_sd = as.numeric(rs[as.character(o$per_delta$delta)]))
  ch$q <- o$q; ch$axis <- o$axis; ch$level <- o$level; ch$eta <- o$eta
  chains_all[[length(chains_all) + 1L]] <- ch
}
per_graph <- do.call(rbind, per_graph); rownames(per_graph) <- NULL
chains_all <- do.call(rbind, chains_all); rownames(chains_all) <- NULL

interp_linear <- function(delta, y, tg) {
  o <- order(delta); delta <- delta[o]; y <- y[o]
  if (any(is.na(y)) || max(y) < tg) return(NA_real_)
  if (y[1] >= tg) return(delta[1])
  i <- which(y >= tg)[1]
  delta[i - 1] + (tg - y[i - 1]) / (y[i] - y[i - 1]) * (delta[i] - delta[i - 1])
}
interp_spline <- function(delta, y, tg) {
  o <- order(delta); delta <- delta[o]; y <- y[o]
  if (any(is.na(y)) || max(y) < tg || y[1] >= tg || any(diff(y) <= 0)) return(NA_real_)
  f <- stats::splinefun(delta, y, method = "hyman")
  stats::uniroot(function(d) f(d) - tg, range(delta))$root
}

cell_key <- interaction(per_graph$q, per_graph$axis, per_graph$level, per_graph$eta, drop = TRUE)
curves <- list(); inversion <- list()
set.seed(4101L)
for (d in split(per_graph, cell_key)) {
  G <- length(unique(d$graph))
  M <- matrix(NA_real_, G, length(deltas), dimnames = list(sort(unique(d$graph)), deltas))
  for (r in seq_len(nrow(d))) M[as.character(d$graph[r]), as.character(d$delta[r])] <- d$median_lmin[r]
  mean_curve <- colMeans(M, na.rm = TRUE)
  se_curve <- apply(M, 2, stats::sd, na.rm = TRUE) / sqrt(colSums(!is.na(M)))
  mcse_med <- tapply(d$mcse_median, d$delta, mean, na.rm = TRUE)
  sd_graph <- tapply(d$median_lmin, d$delta, stats::sd, na.rm = TRUE)
  curves[[length(curves) + 1L]] <- data.frame(
    q = d$q[1], axis = d$axis[1], level = d$level[1], eta = d$eta[1], delta = deltas, n_graphs = G,
    mean_median_lmin = as.numeric(mean_curve), se_mean = as.numeric(se_curve),
    mean_mcse_median = as.numeric(mcse_med[as.character(deltas)]), sd_median_across_graphs = as.numeric(sd_graph[as.character(deltas)]),
    max_rhat = as.numeric(tapply(d$rhat_lmin, d$delta, max, na.rm = TRUE)[as.character(deltas)]),
    min_ess = as.numeric(tapply(d$min_ess_lmin, d$delta, min, na.rm = TRUE)[as.character(deltas)]),
    n_failed = as.integer(tapply(d$n_failed, d$delta, sum)[as.character(deltas)]),
    var_ratio = as.numeric(tapply(d$var_ratio, d$delta, mean, na.rm = TRUE)[as.character(deltas)]),
    rho_sd = as.numeric(tapply(d$rho_sd, d$delta, mean, na.rm = TRUE)[as.character(deltas)]),
    mean_degree = 2 * mean(tapply(d$n_edges, d$graph, `[`, 1)) / d$q[1])
  ds_lin <- interp_linear(deltas, mean_curve, target); ds_spl <- interp_spline(deltas, mean_curve, target)
  boot <- if (G >= 2) vapply(seq_len(n_boot), function(b) {
    interp_linear(deltas, colMeans(M[sample.int(G, G, replace = TRUE), , drop = FALSE], na.rm = TRUE), target)
  }, numeric(1)) else NA_real_
  per_graph_ds <- apply(M, 1, function(y) interp_linear(deltas, y, target))
  bracket <- if (is.na(ds_lin)) c(NA, NA) else { i <- which(mean_curve >= target)[1]; c(deltas[i - 1], deltas[i]) }
  inversion[[length(inversion) + 1L]] <- data.frame(
    q = d$q[1], axis = d$axis[1], level = d$level[1], eta = d$eta[1], n_graphs = G,
    delta_star = ds_lin, delta_star_spline = ds_spl, bracket_lo = bracket[1], bracket_hi = bracket[2],
    boot_lo = if (all(is.na(boot))) NA_real_ else stats::quantile(boot, 0.025, na.rm = TRUE, names = FALSE),
    boot_hi = if (all(is.na(boot))) NA_real_ else stats::quantile(boot, 0.975, na.rm = TRUE, names = FALSE),
    boot_na = sum(is.na(boot)), graph_ds_q10 = stats::quantile(per_graph_ds, 0.10, na.rm = TRUE, names = FALSE),
    graph_ds_q90 = stats::quantile(per_graph_ds, 0.90, na.rm = TRUE, names = FALSE), graph_ds_na = sum(is.na(per_graph_ds)),
    unbracketed = is.na(ds_lin))
}
curves <- do.call(rbind, curves); rownames(curves) <- NULL
inversion <- do.call(rbind, inversion); rownames(inversion) <- NULL
inversion$interp_diff <- inversion$delta_star_spline - inversion$delta_star
inversion <- inversion[order(inversion$axis, inversion$level, inversion$eta, inversion$q), ]

cat("\n== required exponent delta* (target: graph-averaged median smallest eigenvalue 0.5) ==\n")
print(transform(inversion, delta_star = round(delta_star, 3), delta_star_spline = round(delta_star_spline, 3),
                boot_lo = round(boot_lo, 3), boot_hi = round(boot_hi, 3), interp_diff = round(interp_diff, 3),
                graph_ds_q10 = round(graph_ds_q10, 2), graph_ds_q90 = round(graph_ds_q90, 2)),
      row.names = FALSE)
near <- curves[curves$mean_median_lmin > 0.3 & curves$mean_median_lmin < 0.8, ]
cat(sprintf("\nChains: %d run, %d failed. Smallest eigenvalue: max split R-hat %.4f, min bulk ESS %.0f (per chain of %d draws).\n",
            nrow(chains_all), sum(!chains_all$completed), max(per_graph$rhat_lmin, na.rm = TRUE),
            min(per_graph$min_ess_lmin, na.rm = TRUE), n_samples))
cat(sprintf("Near the target (mean median in (0.3, 0.8), %d cell-exponent points): mean MCSE of a graph's median %.4f, sd of medians across graphs %.4f (median over points), SE of the mean curve at most %.4f.\n",
            nrow(near), mean(near$mean_mcse_median), stats::median(near$sd_median_across_graphs), max(near$se_mean)))
cat(sprintf("Unbracketed cells: %d; max |spline - linear| at delta*: %.3f; bootstrap resamples without a bracket: %d.\n",
            sum(inversion$unbracketed), max(abs(inversion$interp_diff), na.rm = TRUE), sum(inversion$boot_na)))
if (any(!chains_all$completed)) {
  cat("\n== failed chains ==\n")
  print(chains_all[!chains_all$completed, c("q", "axis", "level", "eta", "graph", "delta", "chain", "seed", "error")], row.names = FALSE)
}
cat(sprintf("Mean chain time: q = 10 %.2f s, q = 20 %.2f s, q = 50 %.2f s.\n",
            mean(chains_all$secs[chains_all$q == 10]), mean(chains_all$secs[chains_all$q == 20]), mean(chains_all$secs[chains_all$q == 50])))

## ---- figure ----
if (mode != "pilot" && nrow(inversion) > 1) {
  source(file.path(PROJECT_ROOT, "R", "scripts", "fig_style.R"))
  inv <- inversion
  inv$level_f <- factor(inv$level)
  inv$eta_f <- factor(inv$eta, levels = etas, labels = c("eta == 0.5", "eta == 1", "eta == 2"))
  ref <- expand.grid(q = seq(10, 50, by = 1), eta_f = levels(inv$eta_f), stringsAsFactors = FALSE)
  ref$delta_ref <- 0.5 * log(ref$q)
  PAL_LEVEL <- c(OKABE_ITO[["blue"]], OKABE_ITO[["orange"]], OKABE_ITO[["vermillion"]])
  panel_regime <- function(dat, axis_name, legend_title, level_labels) {
    d <- subset(dat, axis == axis_name); d$level_f <- droplevels(d$level_f)
    ggplot(d, aes(x = q, y = delta_star, colour = level_f, group = level_f)) +
      geom_line(data = ref, aes(x = q, y = delta_ref), inherit.aes = FALSE, linetype = "22",
                colour = OKABE_ITO[["grey"]], linewidth = LW_REF) +
      geom_errorbar(aes(ymin = boot_lo, ymax = boot_hi), width = 1.5, linewidth = LW_GUIDE, na.rm = TRUE) +
      geom_line(linewidth = LW_MAIN, na.rm = TRUE) +
      geom_point(size = PT_MAIN, na.rm = TRUE) +
      facet_wrap(~eta_f, nrow = 1, labeller = label_parsed) +
      scale_colour_manual(values = stats::setNames(PAL_LEVEL, levels(d$level_f)), labels = level_labels, name = legend_title) +
      scale_x_continuous(breaks = c(10, 20, 50)) +
      labs(x = LAB_Q, y = expression("required tilt"~delta^"*")) +
      theme_fig()
  }
  pA <- panel_regime(inv, "degree", expression("expected degree"~bar(d)), c("2", "4", "8"))
  pB <- panel_regime(inv, "density", expression("edge density"~rho), c("0.05", "0.10", "0.20"))
  fig <- (pA / pB) + plot_annotation(tag_levels = "A")
  save_fig(file.path(PROJECT_ROOT, "figures", "section_4_1.pdf"), fig, height = fig_height(n_rows = 2L, n_legend_rows = 2L))
  cat("wrote figures/section_4_1.pdf\n")
}

saveRDS(list(inversion = inversion, curves = curves, per_graph = per_graph, chains = chains_all,
             design = list(qs = qs, regimes = regimes, etas = etas, deltas = deltas, sigma = sigma, n_graphs_q = n_graphs_q,
                           pilot_cells = pilot_cells, n_graphs_pilot = n_graphs_pilot,
                           n_chains = n_chains, n_warmup = n_warmup, n_samples = n_samples, target = target, n_boot = n_boot,
                           graph_seed_rule = "410000 + 100 q + regime index (degree 2, 4, 8, density 0.05, 0.10, 0.20)",
                           chain_seed_rule = "(((qi-1) 6 + (ri-1)) 3 + (ei-1)) 1e7 + delta index 1e6 + graph 10 + chain",
                           bgms_mapping = "normal_prior(scale = sigma/2), exponential_prior(eta = eta sigma), spec = conditional, NUTS",
                           nuts_diagnostics = "not exported by bgms 0.2.0.1 for this sampler (no divergence or tree-depth output)"),
             mode = mode, bgms_version = as.character(packageVersion("bgms")), session = sessionInfo(), date = Sys.time()),
        file.path(PROJECT_ROOT, "R", "results", "section_4_1.rds"))
cat(sprintf("\n[%s] done in %.2f h\n", format(Sys.time()), as.numeric(difftime(Sys.time(), t_all, units = "hours"))))
