# Check of the bgms 0.2.0.1 fixed-graph prior sampler (sample_ggm_prior with
# spec = "conditional": NUTS on a free-element Cholesky parametrization) before
# the Section 3.6 marginal study is moved to it (2026-11-14).
#
# Target at a fixed graph Gamma, in our notation (sigma = beta = 1, eta = 1):
#   p(K | Gamma) propto prod_l exp(-k_ll) prod_{(i,j) in Gamma} N(k_ij; 0, 1)
#                       |K|^delta 1{K PD}, k_ij = 0 off the graph.
# bgms mapping (same as for the joint spec, see bgms_param_check_q5.R):
# normal_prior(scale = 0.5), exponential_prior(eta = 1); output on the K scale.
#
# Exact reference at q = 3: importance sampling from the specified entry-wise
# prior with weight |K|^delta 1{K PD} (4e6 proposals). Two graphs: the full
# graph and the chain 1-2-3 (edge (1,3) absent), which exercises the
# zero-constraint path of the Cholesky parametrization. delta in {0, 1, 2}.
# Compared: E[k_ll], E[k_12], E[k_12^2], Var[k_12], with Monte Carlo errors.
#
# Also records run time and effective sample sizes of one chain at q = 10 and
# q = 20 on a Bernoulli graph with edge (1,2) forced in, to size the study.
#
# Run:    R_PROFILE_USER=/dev/null Rscript R/scripts/section36_bgms_conditional_check.R
# Output: R/results/section36_bgms_conditional_check.rds (+ .log)

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
# bgms 0.2.0.1 (GitHub Bayesian-Graphical-Modelling-Lab/bgms, branch develop,
# commit 2da4bf30): remotes::install_github("Bayesian-Graphical-Modelling-Lab/bgms@2da4bf30").
# A project-local library R/lib/bgms_0.2.0.1 is used when present.
local_lib <- file.path(PROJECT_ROOT, "R", "lib", "bgms_0.2.0.1")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages(library(bgms))
if (packageVersion("bgms") != "0.2.0.1") {
  message("bgms ", packageVersion("bgms"), " found; the archived results used 0.2.0.1 (commit 2da4bf30)")
}
has_posterior <- requireNamespace("posterior", quietly = TRUE)

set.seed(36001L)
deltas <- c(0, 1, 2)
n_is <- 4000000L
n_warmup <- 1000L
n_samples <- 50000L

batch_mcse <- function(x, n_batch = 50L) {
  b <- floor(length(x) / n_batch)
  means <- vapply(seq_len(n_batch), function(i) mean(x[((i - 1L) * b + 1L):(i * b)]),
                  numeric(1))
  stats::sd(means) / sqrt(n_batch)
}
w_mean <- function(f, w) {
  est <- sum(w * f) / sum(w)
  se <- sqrt(sum((w * (f - est))^2)) / sum(w)
  c(est = est, se = se)
}

## ---- exact reference by importance sampling, q = 3 ----
graphs <- list(full = matrix(1L, 3, 3), chain = matrix(c(1L, 1L, 0L, 1L, 1L, 1L, 0L, 1L, 1L), 3, 3))
k11 <- stats::rexp(n_is); k22 <- stats::rexp(n_is); k33 <- stats::rexp(n_is)
k12 <- stats::rnorm(n_is); k23 <- stats::rnorm(n_is); k13_full <- stats::rnorm(n_is)
ref <- list()
for (g in names(graphs)) {
  k13 <- if (g == "full") k13_full else 0
  det3 <- k11 * (k22 * k33 - k23^2) - k12 * (k12 * k33 - k23 * k13) + k13 * (k12 * k23 - k22 * k13)
  pd <- (k11 * k22 - k12^2 > 0) & (det3 > 0)
  for (d in deltas) {
    w <- ifelse(pd, det3^d, 0)
    ref[[length(ref) + 1L]] <- data.frame(
      graph = g, delta = d, acc = mean(pd), ess_w = sum(w)^2 / sum(w^2),
      k11 = w_mean(k11, w)[1], k11_se = w_mean(k11, w)[2],
      k22 = w_mean(k22, w)[1], k22_se = w_mean(k22, w)[2],
      k33 = w_mean(k33, w)[1], k33_se = w_mean(k33, w)[2],
      k12 = w_mean(k12, w)[1], k12_se = w_mean(k12, w)[2],
      k12sq = w_mean(k12^2, w)[1], k12sq_se = w_mean(k12^2, w)[2])
  }
}
ref <- do.call(rbind, ref); rownames(ref) <- NULL
ref$var_k12 <- ref$k12sq - ref$k12^2
rm(k11, k22, k33, k12, k23, k13_full)

## ---- bgms conditional chains, q = 3 ----
run_cond <- function(q, E, delta, seed, n_samples, n_warmup) {
  t0 <- Sys.time()
  res <- sample_ggm_prior(p = q, n_samples = n_samples, n_warmup = n_warmup,
                          interaction_prior = normal_prior(scale = 0.5),
                          precision_scale_prior = exponential_prior(eta = 1),
                          delta = delta, spec = "conditional", edge_indicators = E,
                          seed = seed, verbose = FALSE)
  res$secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  res
}
bg <- list()
for (g in names(graphs)) for (d in deltas) {
  res <- run_cond(3L, graphs[[g]], d, seed = 1000L + 10L * d + (g == "chain"), n_samples, n_warmup)
  K12 <- res$K_offdiag[, 1]
  bg[[length(bg) + 1L]] <- data.frame(
    graph = g, delta = d, secs = res$secs,
    k11 = mean(res$K_diag[, 1]), k11_se = batch_mcse(res$K_diag[, 1]),
    k22 = mean(res$K_diag[, 2]), k22_se = batch_mcse(res$K_diag[, 2]),
    k33 = mean(res$K_diag[, 3]), k33_se = batch_mcse(res$K_diag[, 3]),
    k12 = mean(K12), k12_se = batch_mcse(K12),
    k12sq = mean(K12^2), k12sq_se = batch_mcse(K12^2),
    k13_zero = if (g == "chain") all(res$K_offdiag[, 2] == 0) else NA,
    ess_k12 = if (has_posterior) posterior::ess_basic(K12) else NA_real_,
    ess_k11 = if (has_posterior) posterior::ess_basic(res$K_diag[, 1]) else NA_real_)
}
bg <- do.call(rbind, bg); rownames(bg) <- NULL
bg$var_k12 <- bg$k12sq - bg$k12^2

cmp <- data.frame(graph = ref$graph, delta = ref$delta,
                  z_k11 = (bg$k11 - ref$k11) / sqrt(bg$k11_se^2 + ref$k11_se^2),
                  z_k22 = (bg$k22 - ref$k22) / sqrt(bg$k22_se^2 + ref$k22_se^2),
                  z_k33 = (bg$k33 - ref$k33) / sqrt(bg$k33_se^2 + ref$k33_se^2),
                  z_k12 = (bg$k12 - ref$k12) / sqrt(bg$k12_se^2 + ref$k12_se^2),
                  z_k12sq = (bg$k12sq - ref$k12sq) / sqrt(bg$k12sq_se^2 + ref$k12sq_se^2),
                  var_ref = ref$var_k12, var_bgms = bg$var_k12)

cat("== importance-sampling reference (q = 3) ==\n")
print(transform(ref, acc = round(acc, 3), ess_w = round(ess_w)), digits = 4, row.names = FALSE)
cat("\n== bgms conditional chains (q = 3, ", n_samples, " draws) ==\n", sep = "")
print(bg, digits = 4, row.names = FALSE)
cat("\n== standardized differences bgms - reference ==\n")
print(cmp, digits = 3, row.names = FALSE)

## ---- timing and mixing at q = 10 and q = 20 ----
make_graph <- function(q, p_inc, seed) {
  set.seed(seed)
  E <- matrix(0L, q, q)
  E[upper.tri(E)] <- as.integer(stats::runif(q * (q - 1) / 2) < p_inc)
  E <- E + t(E); diag(E) <- 1L; E[1, 2] <- E[2, 1] <- 1L
  E
}
tim <- list()
for (q in c(10L, 20L)) for (d in c(0, 2)) {
  E <- make_graph(q, c(`10` = 0.20, `20` = 0.10)[as.character(q)], seed = 500L + q)
  res <- run_cond(q, E, d, seed = 2000L + q + d, n_samples = 2000L, n_warmup = 1000L)
  K12 <- res$K_offdiag[, 1]
  tim[[length(tim) + 1L]] <- data.frame(
    q = q, delta = d, n_edges = (sum(E) - q) / 2, secs = res$secs,
    ess_k12 = if (has_posterior) posterior::ess_basic(K12) else NA_real_,
    ess_k12sq = if (has_posterior) posterior::ess_basic(K12^2) else NA_real_,
    ess_k11 = if (has_posterior) posterior::ess_basic(res$K_diag[, 1]) else NA_real_,
    mean_diag = mean(res$K_diag), var_k12 = stats::var(K12))
}
tim <- do.call(rbind, tim); rownames(tim) <- NULL
cat("\n== one chain, 1000 warmup + 2000 draws ==\n")
print(tim, digits = 4, row.names = FALSE)

saveRDS(list(reference = ref, bgms = bg, comparison = cmp, timing = tim,
             n_is = n_is, n_samples = n_samples, n_warmup = n_warmup,
             bgms_version = as.character(packageVersion("bgms"))),
        file.path(PROJECT_ROOT, "R", "results", "section36_bgms_conditional_check.rds"))
cat("\nSaved R/results/section36_bgms_conditional_check.rds\n")
