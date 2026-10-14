# Timing study for Section sec:timing (Computational scaling) and
# tables/timing_scaling.tex.
#
# Joint chain (data present, n = 2q, sparse regime p_inc = 4/(q - 1),
# delta = 0.5 log q, eta = 1): wall-clock per within-model sweep, per
# candidate toggle (rejected path), and the additional cost of an accepted
# toggle (the rank-two Sigma update), from the chain's internal block timers.
#
# Hierarchical between-move: per-candidate cost of the coupled estimator
# (2N auxiliary Cholesky completions) as a function of N and q, measured on
# standalone calls at an Erdos-Renyi graph of expected degree 4.
#
# 2026-10-13: replicate runs added (N_REP independent chains / timing batches
# per cell, table reports the mean, the log and rds also the min and max), and
# the measurement environment is recorded in the rds and printed. The archived
# 2026-07-28 single-run version is in R/results/superseded_2026-10-13/.
#
# Output: R/results/timing_scaling.rds (+ log), tables/timing_scaling.tex

suppressPackageStartupMessages(library(Rcpp))
PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
sourceCpp(file.path(PROJECT_ROOT, "R", "src", "spikeslab_sampler.cpp"))

N_REP <- 5L
env <- list(
  date = format(Sys.time(), "%Y-%m-%d %H:%M"),
  cpu = tryCatch(system("sysctl -n machdep.cpu.brand_string", intern = TRUE), error = function(e) NA_character_),
  os = paste(Sys.info()[["sysname"]], Sys.info()[["release"]]),
  r_version = R.version.string,
  blas = extSoftVersion()[["BLAS"]],
  lapack = La_library(),
  compiler = paste(system("R CMD config CXX", intern = TRUE), system("R CMD config CXXFLAGS", intern = TRUE)),
  rcpp = as.character(packageVersion("Rcpp")),
  rcpparmadillo = as.character(packageVersion("RcppArmadillo")),
  threads = "single R process; reference BLAS; no OpenMP in the sampler")
cat("== environment ==\n"); for (k in names(env)) cat(sprintf("  %-14s %s\n", k, env[[k]]))

qs <- c(20L, 50L, 100L)
Ns <- c(100L, 1000L, 10000L)

## ---- joint chain block timing ----
rows <- list()
for (rep in seq_len(N_REP)) for (q in qs) {
  set.seed(8801L + rep)
  n <- 2L * q
  Y <- matrix(rnorm(n * q), n, q)
  S <- crossprod(Y)
  delta <- 0.5 * log(q)
  n_iter <- 200L
  ch <- sampler_chain_cpp(S = S, n_obs = n, q = q, delta = delta, beta = 1,
                          sigma = 1, p_inc = 4 / (q - 1), spec = 0L,
                          log_z = numeric(0), n_iter = n_iter, n_burn = 50L,
                          refresh_every = 50L, track_masks = FALSE,
                          n_aux = 0L, stat_thin = 0L, time_blocks = TRUE)
  tm <- ch$timing
  M <- q * (q - 1) / 2
  rows[[length(rows) + 1L]] <- data.frame(
    rep = rep, q = q, n_att = unname(tm["n_att"]), n_acc = unname(tm["n_acc"]),
    sweep_ms = 1e3 * tm["sweep_s"] / n_iter,
    node_ms = 1e3 * tm["sweep_s"] / (n_iter * q),
    cand_us = 1e6 * (tm["toggle_s"] - tm["accept_s"]) / tm["n_att"],
    accept_extra_us = 1e6 * tm["accept_s"] / max(tm["n_acc"], 1),
    acc_rate = tm["n_acc"] / tm["n_att"],
    E_m = ch$E_m)
}
chain_rep <- do.call(rbind, rows)
rownames(chain_rep) <- NULL
metrics <- c("sweep_ms", "node_ms", "cand_us", "accept_extra_us", "acc_rate", "E_m")
chain_tab <- aggregate(chain_rep[, metrics], by = list(q = chain_rep$q), FUN = mean)
chain_min <- aggregate(chain_rep[, metrics], by = list(q = chain_rep$q), FUN = min)
chain_max <- aggregate(chain_rep[, metrics], by = list(q = chain_rep$q), FUN = max)
cat(sprintf("== joint chain block timing, mean over %d replicate chains (delta = 0.5 log q, p_inc = 4/(q-1), n = 2q, %d iterations) ==\n", N_REP, 200L))
print(chain_tab, row.names = FALSE, digits = 3)
cat("-- min over replicates --\n"); print(chain_min, row.names = FALSE, digits = 3)
cat("-- max over replicates --\n"); print(chain_max, row.names = FALSE, digits = 3)
cat("-- candidate evaluations per chain (n_att) and accepted toggles (n_acc), replicate 1 --\n")
print(chain_rep[chain_rep$rep == 1, c("q", "n_att", "n_acc")], row.names = FALSE)

## ---- hierarchical estimator per-candidate cost ----
er_graph <- function(q, seed) {
  set.seed(seed)
  G <- matrix(0L, q, q)
  p <- 4 / (q - 1)
  for (a in 1:(q - 1)) for (b in (a + 1):q)
    if (runif(1) < p) { G[a, b] <- 1L; G[b, a] <- 1L }
  G
}
est_rows <- list()
for (trep in seq_len(N_REP)) for (q in qs) {
  G <- er_graph(q, 8901L)
  absent <- which(upper.tri(matrix(0, q, q)) & G == 0, arr.ind = TRUE)
  pick <- absent[1 + (seq_len(5) * 37) %% nrow(absent), , drop = FALSE]
  for (N in Ns) {
    reps <- max(2L, as.integer(2e6 / N / q))
    set.seed(8902L + trep)
    t0 <- Sys.time()
    for (r in seq_len(reps)) {
      e <- pick[1 + (r - 1) %% nrow(pick), ]
      invisible(akm_ratio_cpp(G, e[1], e[2], 0.5 * log(q), 1, 1, N))
    }
    per_call <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / reps
    est_rows[[length(est_rows) + 1L]] <- data.frame(
      trep = trep, q = q, N = N, reps = reps, cand_ms = 1e3 * per_call)
  }
}
est_rep <- do.call(rbind, est_rows)
est_tab <- aggregate(est_rep["cand_ms"], by = list(q = est_rep$q, N = est_rep$N), FUN = mean)
est_tab$reps <- aggregate(est_rep["reps"], by = list(q = est_rep$q, N = est_rep$N), FUN = max)$reps
est_tab$cand_ms_min <- aggregate(est_rep["cand_ms"], by = list(q = est_rep$q, N = est_rep$N), FUN = min)$cand_ms
est_tab$cand_ms_max <- aggregate(est_rep["cand_ms"], by = list(q = est_rep$q, N = est_rep$N), FUN = max)$cand_ms
est_tab <- est_tab[order(est_tab$q, est_tab$N), ]
cat(sprintf("\n== hierarchical per-candidate cost (coupled estimator, ms per candidate; mean, min, max over %d timing batches of reps calls) ==\n", N_REP))
print(est_tab, row.names = FALSE, digits = 3)

saveRDS(list(chain = chain_tab, chain_min = chain_min, chain_max = chain_max, chain_replicates = chain_rep,
             estimator = est_tab, estimator_replicates = est_rep, env = env, n_rep = N_REP),
        file.path(PROJECT_ROOT, "R", "results", "timing_scaling.rds"))

## ---- table ----
fmt <- function(x, d) formatC(x, format = "f", digits = d)
lines <- c("\\begin{tabular}{rrrrrrrr}", "\\toprule",
  " & \\multicolumn{4}{c}{joint chain} & \\multicolumn{3}{c}{hierarchical, per candidate (ms)} \\\\",
  "\\cmidrule(lr){2-5}\\cmidrule(lr){6-8}",
  "$q$ & sweep (ms) & node (ms) & candidate ($\\mu$s) & accept extra ($\\mu$s) & $N = 10^2$ & $N = 10^3$ & $N = 10^4$ \\\\",
  "\\midrule")
for (q in qs) {
  ct <- chain_tab[chain_tab$q == q, ]
  et <- est_tab[est_tab$q == q, ]
  lines <- c(lines, sprintf("%d & %s & %s & %s & %s & %s & %s & %s \\\\",
    q, fmt(ct$sweep_ms, 2), fmt(ct$node_ms, 3), fmt(ct$cand_us, 2),
    fmt(ct$accept_extra_us, 1),
    fmt(et$cand_ms[et$N == 100], 2), fmt(et$cand_ms[et$N == 1000], 1),
    fmt(et$cand_ms[et$N == 10000], 0)))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(PROJECT_ROOT, "tables", "timing_scaling.tex"))
cat("\nSaved tables/timing_scaling.tex and R/results/timing_scaling.rds\n")
