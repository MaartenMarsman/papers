# Shared machinery for the known-truth calibration study of indicator and
# mixed diagnostics (see the simulation appendix of the manuscript).
#
# Three groups of helpers:
#   (1) truth generation: networks, threshold calibration, datasets;
#   (2) per-edge observables from a fitted bgms object;
#   (3) checkpointed fit bookkeeping so the grid is resumable.
#
# The two-state machinery (split_rhat, chain_pip_stats,
# indicator_agreement_stats) is reused from R/rhat-diagnostics-helpers.R
# rather than reimplemented.

source("R/rhat-diagnostics-helpers.R")

# --- Paths ---------------------------------------------------------------------

# Overridable so the analysis can be exercised end to end against fabricated
# checkpoints in a scratch directory without touching the real run.
study_dir <- Sys.getenv("THRESHOLD_STUDY_DIR", "results/threshold-study")

study_paths <- function() {
  paths <- list(
    root = study_dir,
    datasets = file.path(study_dir, "datasets"),
    fits = file.path(study_dir, "fits"),
    logs = file.path(study_dir, "logs")
  )
  for (p in paths) {
    if (!dir.exists(p)) dir.create(p, recursive = TRUE)
  }
  paths
}

# --- Inclusion Bayes factors: scale hygiene -------------------------------------

# Stored inclusion Bayes factors are read on the natural-log scale whichever
# scale bgms::extract_inclusion_bf() returns by default: older builds return
# the log, newer builds return the Bayes factor itself with log = TRUE
# recovering the log. This helper is explicit about the scale rather than
# trusting a default, and works on either build. Checkpoints produced with
# earlier builds are put on one scale by R/26-recompute-inclusion-bf.R.
#
# The log scale is also the better one on its own terms: the Bayes factor
# scale saturates at double precision (a log above ~709.78 nats reports +Inf,
# a log of -Inf reports 0), while the accumulators are exact on the log scale
# everywhere.
extract_log_inclusion_bf <- function(fit) {
  if ("log" %in% names(formals(bgms::extract_inclusion_bf))) {
    log_bf <- bgms::extract_inclusion_bf(fit, log = TRUE)
    bayes_factor <- bgms::extract_inclusion_bf(fit, log = FALSE)
    # Assert the two scales are consistent where the BF scale has not
    # saturated, rather than assuming which one came back.
    comparable <- is.finite(log_bf) & is.finite(bayes_factor) &
      abs(log_bf) < 700 & bayes_factor > 0
    if (any(comparable) &&
        !isTRUE(all.equal(exp(log_bf[comparable]), bayes_factor[comparable],
                          tolerance = 1e-6))) {
      stop("extract_inclusion_bf(): log and Bayes factor scales disagree",
           call. = FALSE)
    }
    return(log_bf)
  }
  # Older build without a log argument: the single return value is the
  # natural log.
  bgms::extract_inclusion_bf(fit)
}

# --- Verdicts ------------------------------------------------------------------

bf_hi <- 10
bf_lo <- 1 / 10

# Three-way verdict from the log10 inclusion Bayes factor, at the primary
# thresholds of Table 1. Secondary thresholds (3, 30) are applied in the
# analysis script from the same log10_bf column.
verdict_from_log10_bf <- function(log10_bf, upper = log10(bf_hi),
                                  lower = log10(bf_lo)) {
  out <- rep("undecided", length(log10_bf))
  out[!is.na(log10_bf) & log10_bf >= upper] <- "presence"
  out[!is.na(log10_bf) & log10_bf <= lower] <- "absence"
  out[is.na(log10_bf)] <- NA_character_
  out
}

# --- (1) Truth generation ------------------------------------------------------

# Edge list in the order bgms names edges: row-major over the upper triangle,
# (1,2), (1,3), ..., (1,p), (2,3), ... `upper.tri()` enumerates COLUMN-major
# and would silently misalign every truth vector against the fitted
# diagnostics (the same trap R/01-worked-example.R flags).
edge_pairs <- function(num_variables) {
  t(utils::combn(num_variables, 2))
}

# Read a symmetric matrix out in bgms edge order.
edge_vector <- function(m) {
  pairs <- edge_pairs(nrow(m))
  m[cbind(pairs[, 1], pairs[, 2])]
}

# Erdos-Renyi graph at the target density with weights drawn from a
# calibration pool. `weight_pool` is a vector of absolute weights (the
# empirical conditional interaction strengths of the calibrating fit) and
# `p_negative` the share of negative signs in that fit.
simulate_network <- function(num_variables, density, weight_pool, p_negative,
                             seed) {
  set.seed(seed)
  adjacency <- matrix(0, num_variables, num_variables)
  upper_index <- which(upper.tri(adjacency))
  num_edges <- round(density * length(upper_index))
  present <- sample(upper_index, num_edges)
  values <- rep(0, length(upper_index))
  drawn <- sample(weight_pool, num_edges, replace = TRUE)
  signs <- sample(c(1, -1), num_edges, replace = TRUE,
                  prob = c(1 - p_negative, p_negative))
  values[match(present, upper_index)] <- signs * drawn
  adjacency[upper_index] <- values
  adjacency <- adjacency + t(adjacency)
  true_weight <- edge_vector(adjacency)
  list(pairwise = adjacency,
       true_weight = true_weight,
       true_edge = true_weight != 0,
       density_realized = num_edges / length(upper_index))
}

# Iteratively adjust category thresholds so that the simulated marginal
# category proportions match those of the calibrating data set.
#
# Needed because the tutorial-scale interaction field is strong: pairing the
# calibrating fit's thresholds with a denser or stronger graph than the one
# they were estimated under drives most mass into the top category (for
# example, density .4 with RMS weight .215 puts 62% of responses in the top
# Wenchuan category against an observed 10%). Calibrating the thresholds
# per network keeps density a property of the graph rather than of the
# marginal skew.
#
# Update on the multinomial-logit scale relative to the base category, damped
# (0.15) because the field is self-reinforcing and an undamped step
# oscillates without converging.
calibrate_thresholds <- function(thresholds, pairwise, target_props,
                                 num_categories, num_calibration = 4000L,
                                 gibbs_iter = 400L, max_iter = 40L,
                                 damping = 0.15, average_last = 8L,
                                 seed = 1L, verbose = TRUE) {
  num_variables <- nrow(thresholds)
  history <- vector("list", max_iter)
  deviation <- numeric(max_iter)
  for (it in seq_len(max_iter)) {
    x <- bgms::simulate_mrf(num_states = num_calibration,
                            num_variables = num_variables,
                            num_categories = num_categories,
                            pairwise = pairwise, main = thresholds,
                            variable_type = "ordinal", iter = gibbs_iter,
                            seed = seed * 1000L + it)
    props <- t(apply(x, 2, function(v) {
      prop.table(table(factor(v, levels = 0:num_categories)))
    }))
    props <- pmax(props, 1e-4)
    step <- (log(target_props[, -1, drop = FALSE]) -
               log(props[, -1, drop = FALSE])) -
      (log(target_props[, 1]) - log(props[, 1]))
    thresholds <- thresholds + damping * step
    history[[it]] <- thresholds
    deviation[it] <- max(abs(props - target_props))
    if (verbose) {
      cat(sprintf("    calibration %2d/%2d  max|dprop| = %.4f\n",
                  it, max_iter, deviation[it]))
    }
  }
  # Average the tail of the trajectory: single iterates carry the Monte Carlo
  # noise of one calibration sample.
  tail_index <- seq.int(max_iter - average_last + 1L, max_iter)
  averaged <- Reduce(`+`, history[tail_index]) / length(tail_index)
  list(thresholds = averaged,
       deviation_final = mean(deviation[tail_index]),
       deviation_trace = deviation)
}

# Partial-correlation network for the GGM arm. Psychological GGMs carry
# partial correlations in roughly the 0.1-0.3 range; at higher densities a
# matrix with those entries need not be positive definite, so the off-diagonal
# block is shrunk until the smallest eigenvalue clears `min_eigen`. The
# shrink factor is recorded and reported with the calibration.
simulate_ggm_network <- function(num_variables, density, pcor_range,
                                 p_negative, seed, min_eigen = 0.05) {
  set.seed(seed)
  pcor <- matrix(0, num_variables, num_variables)
  upper_index <- which(upper.tri(pcor))
  num_edges <- round(density * length(upper_index))
  present <- sample(upper_index, num_edges)
  values <- rep(0, length(upper_index))
  drawn <- stats::runif(num_edges, pcor_range[1], pcor_range[2])
  signs <- sample(c(1, -1), num_edges, replace = TRUE,
                  prob = c(1 - p_negative, p_negative))
  values[match(present, upper_index)] <- signs * drawn
  pcor[upper_index] <- values
  pcor <- pcor + t(pcor)

  shrink <- 1
  repeat {
    omega <- diag(num_variables) - shrink * pcor
    eigen_min <- min(eigen(omega, symmetric = TRUE, only.values = TRUE)$values)
    if (eigen_min > min_eigen) break
    shrink <- shrink * 0.95
    if (shrink < 0.05) stop("GGM shrinkage failed to reach a PD matrix",
                            call. = FALSE)
  }
  # With a unit diagonal the partial correlation of an edge is exactly the
  # shrunk off-diagonal entry, so the truth vector is read straight off the
  # matrix (in bgms edge order).
  omega <- diag(num_variables) - shrink * pcor
  true_pcor <- edge_vector(shrink * pcor)
  list(precision = omega,
       true_weight = true_pcor,
       true_edge = true_pcor != 0,
       shrink = shrink,
       eigen_min = eigen_min,
       density_realized = num_edges / length(upper_index))
}

# --- (2) Per-edge observables --------------------------------------------------

# Second-order transition counts, pooled over chains: the input for the
# first-order-Markov validity check that the fragile flag's standard error
# assumes. Returns counts n_{xy->z} in the
# order 00->1, 01->1, 10->1, 11->1 together with their denominators.
second_order_counts <- function(x) {
  n <- length(x)
  if (n < 3L) return(c(n_00 = 0, n_001 = 0, n_01 = 0, n_011 = 0,
                       n_10 = 0, n_101 = 0, n_11 = 0, n_111 = 0))
  a <- x[seq_len(n - 2L)]
  b <- x[seq.int(2L, n - 1L)]
  z <- x[seq.int(3L, n)]
  c(n_00 = sum(a == 0L & b == 0L), n_001 = sum(a == 0L & b == 0L & z == 1L),
    n_01 = sum(a == 0L & b == 1L), n_011 = sum(a == 0L & b == 1L & z == 1L),
    n_10 = sum(a == 1L & b == 0L), n_101 = sum(a == 1L & b == 0L & z == 1L),
    n_11 = sum(a == 1L & b == 1L), n_111 = sum(a == 1L & b == 1L & z == 1L))
}

# Pooled Jeffreys-smoothed two-state summary over a chain set: interior
# PIP, its logit standard error, and the smoothed ESS. Boundary-stable, so it
# is defined for saturated (zero-flip) edges where the empirical variance and
# hence the raw ESS collapse.
pooled_two_state <- function(draws) {
  n_iter <- nrow(draws)
  n01 <- n10 <- n0 <- n1 <- 0
  for (ch in seq_len(ncol(draws))) {
    x <- draws[, ch]
    x_from <- x[-n_iter]
    x_to <- x[-1L]
    n0 <- n0 + sum(x_from == 0L)
    n1 <- n1 + sum(x_from == 1L)
    n01 <- n01 + sum(x_from == 0L & x_to == 1L)
    n10 <- n10 + sum(x_from == 1L & x_to == 0L)
  }
  a_hat <- (n01 + 0.5) / (n0 + 1)
  b_hat <- (n10 + 0.5) / (n1 + 1)
  p_hat <- a_hat / (a_hat + b_hat)
  ess <- n_iter * ncol(draws) * (a_hat + b_hat) / (2 - a_hat - b_hat)
  se_logit <- sqrt(1 / (ess * p_hat * (1 - p_hat)))
  c(a = a_hat, b = b_hat, p = p_hat, ess = ess, se_logit = se_logit,
    flips = n01 + n10)
}

# Vectorized two-state chain simulator: k chains of length n_iter,
# started from stationarity.
simulate_indicator_chains <- function(k, n_iter, p, lambda) {
  a <- (1 - lambda) * p
  b <- (1 - lambda) * (1 - p)
  states <- matrix(0L, nrow = k, ncol = n_iter)
  s <- stats::rbinom(k, 1L, p)
  states[, 1L] <- s
  for (t in 2:n_iter) {
    s <- stats::rbinom(k, 1L, a + s * (1 - b - a))
    states[, t] <- s
  }
  states
}

agreement_stat_from_draws <- function(draws) {
  per_chain <- vapply(seq_len(ncol(draws)),
                      function(ch) chain_pip_stats(draws[, ch]),
                      c(p_hat = 0, var_p = 0, flips = 0, ess_smooth = 0))
  indicator_agreement_stats(per_chain)
}

# Null distribution of the agreement statistic for a (p, lambda, T, chains)
# cell, cached in `env`. The naive chi-square reference false-alarms at up to
# 66% in the few-flip regime, so the test is calibrated by parametric
# bootstrap from the fitted pooled two-state chain. Simulating a fresh null
# per edge would cost more than the fits: the cell key rounds p on the logit
# scale and lambda on the log(1 - lambda) scale, so edges in the same regime
# share one null. The null varies smoothly in both, so the rounding is a
# resolution choice, not an approximation of a different quantity.
agreement_null_cached <- function(env, p, lambda, n_iter, n_chains,
                                  n_boot = 99L, logit_step = 0.25,
                                  lambda_step = 0.25) {
  p <- min(max(p, 1e-6), 1 - 1e-6)
  lambda <- min(max(lambda, 0), 1 - 1e-6)
  key_logit <- round(stats::qlogis(p) / logit_step) * logit_step
  key_lambda <- round(log1p(-lambda) / lambda_step) * lambda_step
  key <- sprintf("%.3f|%.3f|%d|%d", key_logit, key_lambda, n_iter, n_chains)
  if (!is.null(env[[key]])) return(env[[key]])
  p_cell <- stats::plogis(key_logit)
  lambda_cell <- 1 - exp(key_lambda)
  states <- simulate_indicator_chains(n_boot * n_chains, n_iter, p_cell,
                                      lambda_cell)
  null_stats <- vapply(seq_len(n_boot), function(bb) {
    rows <- seq.int((bb - 1L) * n_chains + 1L, bb * n_chains)
    agreement_stat_from_draws(t(states[rows, , drop = FALSE]))
  }, numeric(1))
  env[[key]] <- null_stats
  null_stats
}

# One row of observables per edge for one fitted model: the posterior
# summaries, the two-state and agreement diagnostics, the per-chain verdicts,
# and the second-order transition
# counts. Diagnostics come from summary(fit) (so the study grades what a
# reader sees); the unmasked manual split-Rhat on the RB draws is recorded
# alongside so the analysis can measure what the zero-flip mask costs.
edge_observables <- function(fit, cache_env = new.env(parent = emptyenv()),
                             n_boot = 99L) {
  fit_summary <- summary(fit)
  indicator_summary <- fit_summary$indicator
  pairwise_summary <- fit_summary$pairwise
  main_summary <- fit_summary$main

  raw <- S7::prop(fit, "raw_samples")
  n_chains <- as.integer(raw$nchains)
  n_iter <- as.integer(raw$niter)
  edge_names <- raw$parameter_names$pairwise
  num_edges <- length(edge_names)

  log_bf_matrix <- extract_log_inclusion_bf(fit)
  log10_bf <- log_bf_matrix[lower.tri(log_bf_matrix)] / log(10)

  # Per-chain log inclusion Bayes factors from the same accumulator identity
  # bgms uses when it pools (see bgms:::rb_log_odds_from_counts).
  per_chain_log_bf <- vapply(seq_len(n_chains), function(ch) {
    counts <- raw$rb_counts[[ch]]
    num <- counts[, 1] + counts[, 4] - counts[, 2]
    den <- counts[, 3] - counts[, 1] + counts[, 2]
    out <- log(num) - log(den)
    out[num <= 0] <- -Inf
    out[den <= 0] <- Inf
    out
  }, numeric(num_edges))

  # Cross-check the extractor against the per-chain accumulators, which are
  # build-independent because they are computed here from rb_counts.
  prior_pip <- bgms::extract_prior_inclusion_probabilities(fit)
  prior_log_odds <- log(prior_pip) - log1p(-prior_pip)
  prior_vec <- prior_log_odds[lower.tri(prior_log_odds)]
  per_chain_log10_bf <- (per_chain_log_bf - prior_vec) / log(10)

  # The continuous parameter class gets its own decision tree,
  # which needs per-parameter estimates to grade against the reference run,
  # not just the fit-level max Rhat and min ESS.
  #
  # Which block holds the continuous parameters depends on the model: the
  # ordinal and binary arms put category thresholds in `main` and leave
  # `quadratic` empty, while the GGM arm leaves `main` empty and puts the 15
  # precision parameters in `quadratic`. Taking `main` alone would have left
  # the GGM arm with no continuous diagnostics at all (and a fit-level max
  # Rhat of -Inf).
  continuous_block <- function(block, label) {
    if (is.null(block) || !nrow(block)) return(NULL)
    data.frame(block = label, parameter = rownames(block),
               mean = block$mean, mcse = block$mcse, sd = block$sd,
               n_eff = block$n_eff, rhat = block$Rhat,
               stringsAsFactors = FALSE)
  }
  main_table <- rbind(continuous_block(main_summary, "main"),
                      continuous_block(fit_summary$quadratic, "quadratic"))
  if (is.null(main_table)) {
    stop("no continuous parameter block found in the fit summary",
         call. = FALSE)
  }
  fit_main_max_rhat <- max(main_table$rhat, na.rm = TRUE)
  fit_main_min_ess <- min(main_table$n_eff, na.rm = TRUE)

  rows <- vector("list", num_edges)
  per_chain_rows <- vector("list", num_edges)
  for (e in seq_len(num_edges)) {
    indicator_draws <- vapply(seq_len(n_chains),
                              function(ch) raw$indicator[[ch]][, e],
                              numeric(n_iter))
    rb_draws <- vapply(seq_len(n_chains),
                       function(ch) raw$rb_inclusion[[ch]][, e],
                       numeric(n_iter))

    per_chain <- vapply(seq_len(n_chains),
                        function(ch) chain_pip_stats(indicator_draws[, ch]),
                        c(p_hat = 0, var_p = 0, flips = 0, ess_smooth = 0))
    pooled <- pooled_two_state(indicator_draws)
    agreement <- indicator_agreement_stats(per_chain)
    # Reference runs are the truth, not a graded procedure, so their agreement
    # test is never read; n_boot = 0 skips a bootstrap that would cost more
    # than the reference fit itself at 12 chains x 25,000 draws.
    if (n_boot > 0L) {
      null_stats <- agreement_null_cached(cache_env, pooled["p"],
                                          1 - pooled["a"] - pooled["b"],
                                          n_iter, n_chains, n_boot = n_boot)
      agreement_p_boot <- (1 + sum(null_stats >= agreement)) / (n_boot + 1)
    } else {
      agreement_p_boot <- NA_real_
    }
    agreement_p_chisq <- stats::pchisq(agreement, df = n_chains - 1,
                                       lower.tail = FALSE)

    logit_pip <- stats::qlogis(pooled["p"])
    distance <- min(abs(logit_pip - log(bf_hi)), abs(logit_pip - log(bf_lo)))
    fragile <- distance < 2 * pooled["se_logit"]

    second_order <- rowSums(vapply(seq_len(n_chains), function(ch) {
      second_order_counts(indicator_draws[, ch])
    }, numeric(8)))

    chain_verdicts <- verdict_from_log10_bf(per_chain_log10_bf[e, ])

    # Keep the per-chain quantities in long form. Chain counts differ between
    # test (4) and reference (12) fits, so they cannot be columns, and a fit
    # is far too expensive to repeat if the analysis later wants a
    # cross-chain statistic that was not anticipated here.
    per_chain_rows[[e]] <- data.frame(
      edge = edge_names[e], chain = seq_len(n_chains),
      pip = colMeans(indicator_draws),
      pip_smoothed = per_chain["p_hat", ],
      rb_mean = colMeans(rb_draws),
      flips = per_chain["flips", ],
      ess_smoothed = per_chain["ess_smooth", ],
      log10_bf = per_chain_log10_bf[e, ],
      verdict = chain_verdicts,
      stringsAsFactors = FALSE
    )

    rows[[e]] <- data.frame(
      edge = edge_names[e],
      # --- posterior summaries the reader sees ---
      pip = indicator_summary$mean[e],
      pip_mcse = indicator_summary$mcse[e],
      pip_sd = indicator_summary$sd[e],
      rb_n_eff = indicator_summary$n_eff[e],
      n_eff_mixt = indicator_summary$n_eff_mixt[e],
      rhat_rb = indicator_summary$Rhat[e],
      n01 = indicator_summary[["n0->1"]][e],
      n10 = indicator_summary[["n1->0"]][e],
      total_flips = indicator_summary[["n0->1"]][e] +
        indicator_summary[["n1->0"]][e],
      # --- unmasked counterparts (what the mask suppresses) ---
      rhat_rb_unmasked = split_rhat(rb_draws),
      rhat_indicator_raw = split_rhat(indicator_draws),
      # --- mixed (model-averaged) weight block ---
      weight_mean = pairwise_summary$mean[e],
      weight_mcse = pairwise_summary$mcse[e],
      weight_sd = pairwise_summary$sd[e],
      weight_n_eff = pairwise_summary$n_eff[e],
      weight_share_incl = pairwise_summary$share_incl[e],
      weight_rhat = pairwise_summary$Rhat[e],
      # --- two-state / boundary-stable tools ---
      pip_smoothed = unname(pooled["p"]),
      se_logit_smoothed = unname(pooled["se_logit"]),
      ess_smoothed = unname(pooled["ess"]),
      fragile = unname(fragile),
      agreement_stat = unname(agreement),
      agreement_p_boot = unname(agreement_p_boot),
      agreement_p_chisq = unname(agreement_p_chisq),
      # --- per-chain agreement inputs ---
      pip_chain_min = min(per_chain["p_hat", ]),
      pip_chain_max = max(per_chain["p_hat", ]),
      pip_chain_range = diff(range(per_chain["p_hat", ])),
      flips_chain_min = min(per_chain["flips", ]),
      chains_verdict_agree = length(unique(chain_verdicts)) == 1L,
      # --- verdict ---
      log10_bf = log10_bf[e],
      verdict = verdict_from_log10_bf(log10_bf[e]),
      # --- first-order-Markov check inputs ---
      as.list(second_order),
      # --- fit-level continuous-parameter diagnostics ---
      main_max_rhat = fit_main_max_rhat,
      main_min_ess = fit_main_min_ess,
      stringsAsFactors = FALSE
    )
  }
  list(observables = do.call(rbind, rows),
       per_chain = do.call(rbind, per_chain_rows),
       main = main_table)
}

# --- (3) Checkpointed bookkeeping ----------------------------------------------

# Which arms are Gaussian graphical models. Matched by prefix rather than by
# exact name, so that the near-boundary GGM arm "ggmnb" is fitted as a
# continuous model too; an equality test against "ggm" would route it to the
# ordinal model, where continuous data cannot be fitted.
arm_is_continuous <- function(arm) startsWith(arm, "ggm")

fit_checkpoint_path <- function(fit_id) {
  file.path(study_paths()$fits, paste0(fit_id, ".rds"))
}

fit_is_done <- function(fit_id) {
  file.exists(fit_checkpoint_path(fit_id))
}

# Run one bgm fit for a manifest row and checkpoint its per-edge observables.
# Returns the elapsed seconds. Fit objects are never stored: a 12-chain
# reference run holds roughly a gigabyte of draws, and only the derived
# per-edge table is needed downstream.
run_and_checkpoint <- function(row, dataset, n_boot = 99L, cores = NULL,
                               force = FALSE) {
  if (!force && fit_is_done(row$fit_id)) return(NA_real_)
  if (is.null(cores)) cores <- row$chains
  started <- Sys.time()
  fit_args <- list(
    x = dataset$x,
    variable_type = if (arm_is_continuous(row$arm)) "continuous" else "ordinal",
    iter = row$iter,
    warmup = row$warmup,
    edge_selection = TRUE,
    edge_prior = bgms::bernoulli_prior(0.5),
    interaction_prior = bgms::normal_prior(scale = 1),
    chains = row$chains,
    cores = cores,
    seed = row$sampler_seed,
    display_progress = "none",
    verbose = FALSE
  )
  if (arm_is_continuous(row$arm)) {
    fit_args$precision_graph_prior <- "hierarchical"
  }
  fit <- do.call(bgms::bgm, fit_args)
  elapsed_fit <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  # The parametric-bootstrap agreement p-value draws from R's own generator;
  # seed it per fit so a checkpoint does not depend on which worker, or in
  # which order, produced it.
  set.seed(row$sampler_seed)
  extracted <- edge_observables(fit, n_boot = n_boot)
  elapsed_total <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  identifiers <- data.frame(
    fit_id = row$fit_id, stage = row$stage, arm = row$arm, kind = row$kind,
    dataset_id = row$dataset_id, density = row$density, n = row$n,
    data_seed = row$data_seed, setting = row$setting, iter = row$iter,
    warmup = row$warmup, chains = row$chains,
    sampler_seed = row$sampler_seed, stringsAsFactors = FALSE
  )
  observables <- cbind(identifiers, extracted$observables)
  per_chain <- cbind(identifiers[rep(1L, nrow(extracted$per_chain)), ],
                     extracted$per_chain)
  main <- cbind(identifiers[rep(1L, nrow(extracted$main)), ], extracted$main)
  # Guard the join between the generating truth and the fitted diagnostics:
  # both must be in bgms edge order, and a silent misalignment here would
  # corrupt every error rate the study reports.
  pairs <- edge_pairs(ncol(dataset$x))
  expected_names <- paste(colnames(dataset$x)[pairs[, 1]],
                          colnames(dataset$x)[pairs[, 2]], sep = "-")
  if (!identical(observables$edge, expected_names)) {
    stop("edge order mismatch between the fit and the generating network",
         call. = FALSE)
  }
  observables$true_weight <- dataset$true_weight
  observables$true_edge <- dataset$true_edge

  # Write the checkpoint atomically, so an interrupted worker never leaves a
  # truncated file that fit_is_done() would mistake for a finished fit.
  final_path <- fit_checkpoint_path(row$fit_id)
  tmp_path <- paste0(final_path, ".tmp")
  saveRDS(list(observables = observables,
               per_chain = per_chain,
               main = main,
               seconds_fit = elapsed_fit,
               seconds_total = elapsed_total,
               finished = Sys.time()),
          tmp_path)
  file.rename(tmp_path, final_path)
  elapsed_total
}
