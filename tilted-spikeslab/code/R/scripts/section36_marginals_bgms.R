# Section 3.6 marginal study in bgms 0.2.0.1 (2026-11-14): realized diagonal
# mean and realized variance of the included reference entry k_12 under the
# fixed-graph (conditional) prior, averaged with equal weights over a sample
# of reference graphs. Replaces the rejection-sampling study
# (marginal_offdiag_sweep.R, archived results R/results/marginal_offdiag.rds).
#
# Design: q in {5, 10, 15, 20}; graphs drawn edge-wise with inclusion
# probability p_inc = 0.50, 0.20, 0.13, 0.10 and edge (1,2) forced in; the same
# graphs are used at every delta in {0, 1, 2}; sigma = beta = 1
# (normal_prior(scale = 0.5), exponential_prior(eta = 1); mapping verified in
# section36_bgms_conditional_check.R). K | graph is sampled with
# sample_ggm_prior(spec = "conditional") (NUTS on a Cholesky parametrization;
# bgms returns no divergence or tree-depth information, so the within-chain
# diagnostics here are ESS and, from the chains, split R-hat).
#
# Per (graph, delta): n_chains chains of n_warmup + n_samples iterations with
# distinct seeds. Per chain: means of k_12, k_12^2 and the mean diagonal entry,
# batch means (n_batch), ESS, run time, completion flag and error message.
# Pre-specified analysis (notes/section36-bgms-marginals.md): the per-graph
# value is the mean over its completed chains (all chains have equal length,
# so graphs keep equal weight); the cell estimate is the mean over graphs; the
# mixture variance is mean(k_12^2) - mean(k_12)^2 across graphs; standard
# errors use the between-graph spread of the per-graph values (delta method
# for the variance), which carries both chain and graph variation. Failed
# chains and R-hat above 1.01 are counted and reported, never dropped silently.
#
# Run:    R_PROFILE_USER=/dev/null Rscript R/scripts/section36_marginals_bgms.R
# Output: R/results/section36_marginals_bgms.rds (+ .log): graphs, per-chain
#         table, per-graph table, per-cell summary, design, bgms version.

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
stopifnot(requireNamespace("posterior", quietly = TRUE))
OUT_RDS <- file.path(PROJECT_ROOT, "R", "results", "section36_marginals_bgms.rds")

## ---- design (fixed from the pilot, section36_marginals_pilot.R) ----
qs <- c(5L, 10L, 15L, 20L)
p_inc_for_q <- c(`5` = 0.50, `10` = 0.20, `15` = 0.13, `20` = 0.10)
deltas <- c(0, 1, 2)
sigma <- 1; beta <- 1
n_graphs <- 200L
n_chains <- 2L
n_warmup <- 1000L
n_samples <- 2000L
n_batch <- 50L
graph_seed <- function(q) 36200L + q
chain_seed <- function(q, d, g, c) as.integer(q * 1000000L + round(d) * 100000L + g * 10L + c)

batch_means <- function(x, n_batch) {
  b <- floor(length(x) / n_batch)
  vapply(seq_len(n_batch), function(i) mean(x[((i - 1L) * b + 1L):(i * b)]), numeric(1))
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
run_chain <- function(q, E, delta, seed) {
  t0 <- Sys.time()
  res <- tryCatch(
    sample_ggm_prior(p = q, n_samples = n_samples, n_warmup = n_warmup,
                     interaction_prior = normal_prior(scale = sigma / 2),
                     precision_scale_prior = exponential_prior(eta = beta * sigma),
                     delta = delta, spec = "conditional", edge_indicators = E,
                     seed = seed, verbose = FALSE),
    error = function(e) conditionMessage(e))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.character(res)) return(list(ok = FALSE, error = res, secs = secs))
  list(ok = TRUE, error = NA_character_, secs = secs,
       k12 = res$K_offdiag[, 1], dbar = rowMeans(res$K_diag))
}

t_all <- Sys.time()
cat(sprintf("[%s] bgms %s; %d graphs x %d chains per (q, delta); %d + %d iterations\n",
            format(Sys.time()), as.character(packageVersion("bgms")), n_graphs, n_chains,
            n_warmup, n_samples))
graphs_all <- list(); chains <- list(); batches <- list(); rhats <- list()
for (q in qs) {
  graphs <- draw_graphs(q, p_inc_for_q[as.character(q)], n_graphs, graph_seed(q))
  graphs_all[[as.character(q)]] <- graphs
  t0 <- Sys.time()
  for (g in seq_len(n_graphs)) for (d in deltas) {
    out <- vector("list", n_chains)
    for (c in seq_len(n_chains)) {
      seed <- chain_seed(q, d, g, c)
      r <- run_chain(q, graphs[[g]], d, seed)
      out[[c]] <- r
      chains[[length(chains) + 1L]] <- data.frame(
        q = q, delta = d, graph = g, chain = c, seed = seed,
        n_edges = (sum(graphs[[g]]) - q) / 2, completed = r$ok, error = r$error, secs = r$secs,
        m_k12 = if (r$ok) mean(r$k12) else NA_real_,
        m_k12sq = if (r$ok) mean(r$k12^2) else NA_real_,
        m_dbar = if (r$ok) mean(r$dbar) else NA_real_,
        ess_k12 = if (r$ok) posterior::ess_basic(r$k12) else NA_real_,
        ess_k12sq = if (r$ok) posterior::ess_basic(r$k12^2) else NA_real_,
        ess_dbar = if (r$ok) posterior::ess_basic(r$dbar) else NA_real_)
      if (r$ok) batches[[length(batches) + 1L]] <- list(
        q = q, delta = d, graph = g, chain = c,
        k12 = batch_means(r$k12, n_batch), k12sq = batch_means(r$k12^2, n_batch),
        dbar = batch_means(r$dbar, n_batch))
    }
    ok <- vapply(out, function(r) r$ok, logical(1))
    if (sum(ok) >= 2L) {
      out <- out[ok]
      rhats[[length(rhats) + 1L]] <- data.frame(
        q = q, delta = d, graph = g,
        rhat_k12 = posterior::rhat(do.call(cbind, lapply(out, `[[`, "k12"))),
        rhat_k12sq = posterior::rhat(do.call(cbind, lapply(out, function(r) r$k12^2))),
        rhat_dbar = posterior::rhat(do.call(cbind, lapply(out, `[[`, "dbar"))))
    }
  }
  cat(sprintf("[%s] q = %d done in %.1f min\n", format(Sys.time()), q,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
chains <- do.call(rbind, chains); rownames(chains) <- NULL
rhats <- do.call(rbind, rhats); rownames(rhats) <- NULL

## ---- per-graph values (mean over completed chains) ----
gkey <- interaction(chains$q, chains$delta, chains$graph, drop = TRUE)
per_graph <- do.call(rbind, lapply(split(chains, gkey), function(d) {
  ok <- d$completed
  data.frame(q = d$q[1], delta = d$delta[1], graph = d$graph[1], n_edges = d$n_edges[1],
             n_completed = sum(ok),
             m_k12 = if (any(ok)) mean(d$m_k12[ok]) else NA_real_,
             m_k12sq = if (any(ok)) mean(d$m_k12sq[ok]) else NA_real_,
             m_dbar = if (any(ok)) mean(d$m_dbar[ok]) else NA_real_)
}))
rownames(per_graph) <- NULL

## ---- per-cell summary ----
ckey <- interaction(per_graph$q, per_graph$delta, drop = TRUE)
summ <- do.call(rbind, lapply(split(per_graph, ckey), function(d) {
  q <- d$q[1]; dl <- d$delta[1]
  ok <- d$n_completed > 0
  G <- sum(ok); m1 <- d$m_k12[ok]; m2 <- d$m_k12sq[ok]; db <- d$m_dbar[ok]
  var_mix <- mean(m2) - mean(m1)^2
  grad <- c(1, -2 * mean(m1))
  se_var <- sqrt(drop(t(grad) %*% stats::cov(cbind(m2, m1)) %*% grad) / G)
  ch <- chains[chains$q == q & chains$delta == dl, ]
  rk <- rhats[rhats$q == q & rhats$delta == dl, ]
  data.frame(q = q, delta = dl, p_inc = p_inc_for_q[as.character(q)],
             n_graphs = nrow(d), n_graphs_used = G, n_failed_chains = sum(!ch$completed),
             diag_ratio = mean(db) * beta, diag_ratio_se = stats::sd(db) / sqrt(G) * beta,
             diag_bound = dl + 1,
             var_ratio = var_mix / sigma^2, var_ratio_se = se_var / sigma^2,
             mean_k12 = mean(m1), mean_k12_se = stats::sd(m1) / sqrt(G),
             min_ess_k12sq = min(ch$ess_k12sq, na.rm = TRUE), min_ess_dbar = min(ch$ess_dbar, na.rm = TRUE),
             max_rhat = if (nrow(rk)) max(rk$rhat_k12, rk$rhat_k12sq, rk$rhat_dbar) else NA_real_,
             n_rhat_above_1.01 = if (nrow(rk)) sum(pmax(rk$rhat_k12, rk$rhat_k12sq, rk$rhat_dbar) > 1.01) else NA_integer_,
             mean_chain_secs = mean(ch$secs))
}))
rownames(summ) <- NULL
summ <- summ[order(summ$q, summ$delta), ]

cat("\n== per-cell summary ==\n")
print(transform(summ, diag_ratio = round(diag_ratio, 3), diag_ratio_se = round(diag_ratio_se, 4),
                var_ratio = round(var_ratio, 3), var_ratio_se = round(var_ratio_se, 4),
                mean_k12 = round(mean_k12, 4), mean_k12_se = round(mean_k12_se, 4),
                min_ess_k12sq = round(min_ess_k12sq), min_ess_dbar = round(min_ess_dbar),
                max_rhat = round(max_rhat, 3), mean_chain_secs = round(mean_chain_secs, 2)),
      row.names = FALSE)
if (any(!chains$completed)) {
  cat("\n== failed chains ==\n")
  print(chains[!chains$completed, c("q", "delta", "graph", "chain", "seed", "error")], row.names = FALSE)
}
cat(sprintf("\n[%s] all done in %.1f h\n", format(Sys.time()),
            as.numeric(difftime(Sys.time(), t_all, units = "hours"))))

saveRDS(list(summary = summ, per_graph = per_graph, chains = chains, rhats = rhats,
             batches = batches, graphs = graphs_all,
             design = list(qs = qs, p_inc = p_inc_for_q, deltas = deltas, sigma = sigma, beta = beta,
                           n_graphs = n_graphs, n_chains = n_chains, n_warmup = n_warmup,
                           n_samples = n_samples, n_batch = n_batch,
                           graph_seed = vapply(qs, graph_seed, integer(1)),
                           chain_seed_rule = "q*1e6 + delta*1e5 + graph*10 + chain"),
             bgms_version = as.character(packageVersion("bgms")),
             r_version = R.version.string, date = format(Sys.Date())),
        OUT_RDS)
cat("Saved", OUT_RDS, "\n")
