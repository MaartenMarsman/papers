# Second-sampler cross-check for the Section 3.6 study (2026-11-14): the
# exported NUTS fixed-graph sampler (sample_ggm_prior, spec = "conditional")
# against bgms's internal row-block Gibbs sweep run on the same fixed graph
# (bgms:::sample_ggm with edge_selection = FALSE, sampler_type = "gibbs",
# n = 0; the edge move is gated off, so the row draws keep the initial graph).
# Same graphs, same priors (sigma = beta = 1), delta in {0, 1, 2}; compared
# per graph on E[k_12], E[k_12^2] and the mean diagonal entry with
# batch-means Monte Carlo errors. Also records the largest absolute value of an
# excluded entry (NUTS enforces the zeros through a null-space construction,
# so they are zero to floating-point precision rather than exactly).
#
# Run:    R_PROFILE_USER=/dev/null Rscript R/scripts/section36_gibbs_crosscheck.R
# Output: R/results/section36_gibbs_crosscheck.rds (+ .log)

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

qs <- c(5L, 10L, 20L)
p_inc_for_q <- c(`5` = 0.50, `10` = 0.20, `20` = 0.10)
deltas <- c(0, 1, 2)
n_graphs <- 4L
n_warmup <- 1000L
n_samples <- 20000L
sigma <- 1; beta <- 1

batch_mcse <- function(x, n_batch = 50L) {
  b <- floor(length(x) / n_batch)
  means <- vapply(seq_len(n_batch), function(i) mean(x[((i - 1L) * b + 1L):(i * b)]), numeric(1))
  stats::sd(means) / sqrt(n_batch)
}
draw_graphs <- function(q, p_inc, n, seed) {
  set.seed(seed)
  lapply(seq_len(n), function(i) {
    E <- matrix(0L, q, q)
    E[upper.tri(E)] <- as.integer(stats::runif(q * (q - 1) / 2) < p_inc)
    E <- E + t(E); diag(E) <- 1L; E[1, 2] <- E[2, 1] <- 1L
    E
  })
}
## index of K(i, j) in bgms's (i <= j) row-wise upper-triangle layout (1-based)
ut_index <- function(p, i, j) {
  e <- 0L
  for (a in seq_len(p)) for (b in a:p) { e <- e + 1L; if (a == i && b == j) return(e) }
}
run_nuts <- function(q, E, d, seed) {
  res <- sample_ggm_prior(p = q, n_samples = n_samples, n_warmup = n_warmup,
                          interaction_prior = normal_prior(scale = sigma / 2),
                          precision_scale_prior = exponential_prior(eta = beta * sigma),
                          delta = d, spec = "conditional", edge_indicators = E,
                          seed = seed, verbose = FALSE)
  list(k12 = res$K_offdiag[, 1], dbar = rowMeans(res$K_diag),
       excl_max = max(abs(res$K_offdiag[, t(E)[lower.tri(E)] == 0L])))  # columns are row-major
}
run_gibbs <- function(q, E, d, seed) {
  ip <- bgms:::unpack_interaction_prior(normal_prior(scale = sigma / 2))
  sp <- bgms:::unpack_scale_prior(exponential_prior(eta = beta * sigma))
  rate <- bgms:::resolve_scale_rate(sp$scale_rate, sp$scale_eta, ip$pairwise_scale)
  inputFromR <- list(n = 0L, suf_stat = matrix(0, q, q), pairwise_scale = ip$pairwise_scale,
                     interaction_prior_type = ip$interaction_prior_type,
                     scale_prior_type = sp$scale_prior_type, scale_shape = sp$scale_shape,
                     scale_rate = rate)
  out <- bgms:::sample_ggm(inputFromR = inputFromR, prior_inclusion_prob = matrix(0.5, q, q),
                           initial_edge_indicators = E, no_iter = n_samples, no_warmup = n_warmup,
                           no_chains = 1L, edge_selection = FALSE, sampler_type = "gibbs",
                           seed = seed, no_threads = 1L, progress_type = 0L,
                           edge_prior = "Bernoulli", delta = d)
  if (isTRUE(out[[1]]$error)) stop(out[[1]]$error_msg)
  S <- out[[1]]$samples  # p(p+1)/2 x n_iter
  diag_idx <- vapply(seq_len(q), function(i) ut_index(q, i, i), integer(1))
  excl <- which(E == 0L & upper.tri(E), arr.ind = TRUE)
  excl_idx <- if (nrow(excl)) apply(excl, 1, function(ij) ut_index(q, ij[1], ij[2])) else integer(0)
  list(k12 = S[ut_index(q, 1L, 2L), ], dbar = colMeans(S[diag_idx, , drop = FALSE]),
       excl_max = if (length(excl_idx)) max(abs(S[excl_idx, ])) else NA_real_, rate = rate)
}

rows <- list()
for (q in qs) {
  graphs <- draw_graphs(q, p_inc_for_q[as.character(q)], n_graphs, seed = 36300L + q)
  for (g in seq_len(n_graphs)) for (d in deltas) {
    E <- graphs[[g]]
    t0 <- Sys.time(); a <- run_nuts(q, E, d, seed = 7000L + 100L * g + d); ta <- as.numeric(Sys.time() - t0, units = "secs")
    t0 <- Sys.time(); b <- run_gibbs(q, E, d, seed = 8000L + 100L * g + d); tb <- as.numeric(Sys.time() - t0, units = "secs")
    z <- function(x, y) (mean(x) - mean(y)) / sqrt(batch_mcse(x)^2 + batch_mcse(y)^2)
    rows[[length(rows) + 1L]] <- data.frame(
      q = q, graph = g, delta = d, n_edges = (sum(E) - q) / 2,
      k12_nuts = mean(a$k12), k12_gibbs = mean(b$k12), z_k12 = z(a$k12, b$k12),
      k12sq_nuts = mean(a$k12^2), k12sq_gibbs = mean(b$k12^2), z_k12sq = z(a$k12^2, b$k12^2),
      dbar_nuts = mean(a$dbar), dbar_gibbs = mean(b$dbar), z_dbar = z(a$dbar, b$dbar),
      excl_max_nuts = a$excl_max, excl_max_gibbs = b$excl_max,
      secs_nuts = ta, secs_gibbs = tb, gibbs_rate = b$rate)
  }
}
rows <- do.call(rbind, rows); rownames(rows) <- NULL
cat("== NUTS vs internal Gibbs on the same fixed graphs (", n_samples, " draws each) ==\n", sep = "")
print(transform(rows, k12_nuts = round(k12_nuts, 3), k12_gibbs = round(k12_gibbs, 3), z_k12 = round(z_k12, 2),
                k12sq_nuts = round(k12sq_nuts, 3), k12sq_gibbs = round(k12sq_gibbs, 3), z_k12sq = round(z_k12sq, 2),
                dbar_nuts = round(dbar_nuts, 3), dbar_gibbs = round(dbar_gibbs, 3), z_dbar = round(z_dbar, 2),
                secs_nuts = round(secs_nuts, 1), secs_gibbs = round(secs_gibbs, 1)), row.names = FALSE)
cat(sprintf("\nmax |z|: k12 %.2f, k12^2 %.2f, diag %.2f over %d comparisons; max |excluded entry|: NUTS %.1e, Gibbs %.1e\n",
            max(abs(rows$z_k12)), max(abs(rows$z_k12sq)), max(abs(rows$z_dbar)), nrow(rows),
            max(rows$excl_max_nuts), max(rows$excl_max_gibbs)))
saveRDS(list(rows = rows, n_samples = n_samples, n_warmup = n_warmup,
             bgms_version = as.character(packageVersion("bgms"))),
        file.path(PROJECT_ROOT, "R", "results", "section36_gibbs_crosscheck.rds"))
