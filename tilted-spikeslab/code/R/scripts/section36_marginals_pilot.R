# Pilot for the Section 3.6 marginal study in bgms 0.2.0.1 (2026-11-14).
#
# Design under test: q in {5, 10, 15, 20}; graphs drawn edge-wise with
# p_inc = 0.50, 0.20, 0.13, 0.10 and edge (1,2) forced in; the same graphs
# across delta in {0, 1, 2}; K | graph sampled with sample_ggm_prior(spec =
# "conditional") (NUTS), sigma = beta = 1 (normal_prior(scale = 0.5),
# exponential_prior(eta = 1)); 1000 warmup + 2000 retained iterations.
#
# Pilot: 40 graphs per q, two chains (different seeds) per graph and delta.
# Recorded per chain: means of k_12, k_12^2, mean diagonal entry; batch-means
# MCSE; ESS; run time; completion. Per (graph, delta): split R-hat over the two
# chains. Used to choose the number of graphs and to check warmup and chain
# length for the main study (section36_marginals_bgms.R).
#
# Run:    R_PROFILE_USER=/dev/null Rscript R/scripts/section36_marginals_pilot.R
# Output: R/results/section36_marginals_pilot.rds (+ .log)

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

qs <- c(5L, 10L, 15L, 20L)
p_inc_for_q <- c(`5` = 0.50, `10` = 0.20, `15` = 0.13, `20` = 0.10)
deltas <- c(0, 1, 2)
n_graphs <- 40L
n_chains <- 2L
n_warmup <- 1000L
n_samples <- 2000L

batch_mcse <- function(x, n_batch = 50L) {
  b <- floor(length(x) / n_batch)
  means <- vapply(seq_len(n_batch), function(i) mean(x[((i - 1L) * b + 1L):(i * b)]),
                  numeric(1))
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
run_chain <- function(q, E, delta, seed) {
  t0 <- Sys.time()
  res <- tryCatch(
    sample_ggm_prior(p = q, n_samples = n_samples, n_warmup = n_warmup,
                     interaction_prior = normal_prior(scale = 0.5),
                     precision_scale_prior = exponential_prior(eta = 1),
                     delta = delta, spec = "conditional", edge_indicators = E,
                     seed = seed, verbose = FALSE),
    error = function(e) conditionMessage(e))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.character(res)) return(list(ok = FALSE, error = res, secs = secs))
  list(ok = TRUE, error = NA_character_, secs = secs,
       k12 = res$K_offdiag[, 1], dbar = rowMeans(res$K_diag))
}

t_all <- Sys.time()
chains <- list(); rhats <- list(); graphs_all <- list()
for (q in qs) {
  graphs <- draw_graphs(q, p_inc_for_q[as.character(q)], n_graphs, seed = 36100L + q)
  graphs_all[[as.character(q)]] <- graphs
  t0 <- Sys.time()
  for (g in seq_len(n_graphs)) for (d in deltas) {
    out <- vector("list", n_chains)
    for (c in seq_len(n_chains)) {
      seed <- as.integer(q * 100000L + round(d) * 10000L + g * 10L + c)
      r <- run_chain(q, graphs[[g]], d, seed)
      out[[c]] <- r
      chains[[length(chains) + 1L]] <- data.frame(
        q = q, delta = d, graph = g, chain = c, seed = seed,
        n_edges = (sum(graphs[[g]]) - q) / 2, completed = r$ok, error = r$error, secs = r$secs,
        m_k12 = if (r$ok) mean(r$k12) else NA_real_,
        m_k12sq = if (r$ok) mean(r$k12^2) else NA_real_,
        m_dbar = if (r$ok) mean(r$dbar) else NA_real_,
        mcse_k12 = if (r$ok) batch_mcse(r$k12) else NA_real_,
        mcse_k12sq = if (r$ok) batch_mcse(r$k12^2) else NA_real_,
        mcse_dbar = if (r$ok) batch_mcse(r$dbar) else NA_real_,
        ess_k12 = if (r$ok) posterior::ess_basic(r$k12) else NA_real_,
        ess_k12sq = if (r$ok) posterior::ess_basic(r$k12^2) else NA_real_,
        ess_dbar = if (r$ok) posterior::ess_basic(r$dbar) else NA_real_)
    }
    if (all(vapply(out, function(r) r$ok, logical(1)))) {
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

## ---- summary per (q, delta): between-graph spread of chain means, MC error ----
key <- interaction(chains$q, chains$delta, drop = TRUE)
summ <- do.call(rbind, lapply(split(chains, key), function(d) {
  d1 <- d[d$chain == 1L & d$completed, ]
  rk <- rhats[rhats$q == d$q[1] & rhats$delta == d$delta[1], ]
  data.frame(q = d$q[1], delta = d$delta[1], n_failed = sum(!d$completed),
             secs = mean(d$secs),
             sd_graph_k12sq = stats::sd(d1$m_k12sq), med_mcse_k12sq = stats::median(d1$mcse_k12sq),
             sd_graph_dbar = stats::sd(d1$m_dbar), med_mcse_dbar = stats::median(d1$mcse_dbar),
             min_ess_k12sq = min(d1$ess_k12sq), min_ess_dbar = min(d1$ess_dbar),
             max_rhat = max(rk$rhat_k12, rk$rhat_k12sq, rk$rhat_dbar),
             var_k12 = mean(d1$m_k12sq) - mean(d1$m_k12)^2, mean_diag = mean(d1$m_dbar))
}))
rownames(summ) <- NULL
summ <- summ[order(summ$q, summ$delta), ]
cat("\n== pilot summary (chain 1 of each graph; sd over graphs of chain means vs median within-chain MCSE) ==\n")
print(summ, digits = 3, row.names = FALSE)
if (any(!chains$completed)) { cat("\n== failures ==\n"); print(chains[!chains$completed, c("q", "delta", "graph", "chain", "error")]) }
cat(sprintf("\n[%s] pilot done in %.1f min\n", format(Sys.time()),
            as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
saveRDS(list(chains = chains, rhats = rhats, summary = summ, graphs = graphs_all,
             design = list(n_graphs = n_graphs, n_chains = n_chains, n_warmup = n_warmup,
                           n_samples = n_samples, p_inc = p_inc_for_q, deltas = deltas),
             bgms_version = as.character(packageVersion("bgms"))),
        file.path(PROJECT_ROOT, "R", "results", "section36_marginals_pilot.rds"))
