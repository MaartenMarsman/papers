# Known-truth calibration study, step 4: truth, rule search, decision trees.
# Design: the simulation appendix of the manuscript.
#
# The quantity every threshold implicitly promises is a CONDITIONAL error
# rate: given that a diagnostic rule passes, how often is the reported answer
# wrong? Nothing in the existing evidence estimates it, because nothing ties
# an observed diagnostic profile to a known truth. This script estimates it.
#
# Truth is the reference run (12 chains, 10,000 warmup + 25,000 sampling) for
# the same data set. Edges whose reference Bayes factor sits within Monte
# Carlo error of a verdict boundary are dropped from the denominator: no
# procedure can be graded there, so their share is reported instead.
#
# What counts as an error differs by parameter class, because what the class
# reports differs:
#   indicator     the three-way verdict, so error = wrong verdict;
#   mixed weight  a value with a stated MCSE, so error = the reference value
#                 lies further away than that MCSE claims;
#   continuous    likewise, for the category thresholds (ordinal, binary) or
#                 the precisions (GGM).
#
# Usage: Rscript R/23-threshold-study-analysis.R [--stage 1]
# Outputs: tables/threshold-study-truth.csv,
#          tables/threshold-study-rules-<class>.csv,
#          tables/threshold-study-stability.csv,
#          tables/threshold-study-markov.csv,
#          tables/threshold-study-bf-sensitivity.csv

source("R/threshold-study-helpers.R")

parse_args <- function(args) {
  out <- list(stage = "1")
  i <- 1L
  while (i <= length(args)) {
    out[[sub("^--", "", args[i])]] <- args[i + 1L]
    i <- i + 2L
  }
  out
}
options <- parse_args(commandArgs(trailingOnly = TRUE))
target_stage <- as.integer(options$stage)
paths <- study_paths()

# Overridable alongside THRESHOLD_STUDY_DIR so a smoke test can write its
# results somewhere harmless.
output_prefix <- Sys.getenv("THRESHOLD_STUDY_OUTPUT", "tables/")
if (!grepl("/$", output_prefix)) output_prefix <- paste0(output_prefix, "/")
dir.create(output_prefix, showWarnings = FALSE, recursive = TRUE)
out_path <- function(name) paste0(output_prefix, name)

# --- Load the stage's checkpoints ----------------------------------------------

manifest <- utils::read.csv(file.path(paths$root, "manifest.csv"),
                            stringsAsFactors = FALSE)
manifest <- manifest[manifest$stage == target_stage, ]

available <- manifest[vapply(manifest$fit_id, fit_is_done, logical(1)), ]
cat(sprintf("stage %d: %d of %d fits checkpointed (%d test, %d reference)\n",
            target_stage, nrow(available), nrow(manifest),
            sum(available$kind == "test"),
            sum(available$kind == "reference")))
if (!any(available$kind == "reference")) {
  stop("no reference runs yet: truth cannot be built", call. = FALSE)
}

read_block <- function(fit_ids, block) {
  do.call(rbind, lapply(fit_ids, function(id) {
    readRDS(fit_checkpoint_path(id))[[block]]
  }))
}

reference_ids <- available$fit_id[available$kind == "reference"]
test_ids <- available$fit_id[available$kind == "test"]

reference <- read_block(reference_ids, "observables")
test <- read_block(test_ids, "observables")
reference_main <- read_block(reference_ids, "main")
test_main <- read_block(test_ids, "main")

# --- Truth, and which edges can be graded at all --------------------------------

# The reference verdict is truth. Its own Monte Carlo uncertainty is taken
# from the boundary-stable smoothed two-state standard error rather than the
# summary's mcse column, because the latter is masked exactly on the settled
# edges (the decisive ones) where a verdict still has to be graded.
reference$se_log10 <- reference$se_logit_smoothed / log(10)
reference$distance_to_boundary <- pmin(abs(reference$log10_bf - 1),
                                       abs(reference$log10_bf + 1))
reference$gradeable <- reference$distance_to_boundary >
  4 * reference$se_log10

truth <- data.frame(
  dataset_id = reference$dataset_id, edge = reference$edge,
  truth_verdict = reference$verdict,
  truth_log10_bf = reference$log10_bf,
  truth_weight = reference$weight_mean,
  truth_pip = reference$pip,
  truth_se_log10 = reference$se_log10,
  gradeable = reference$gradeable,
  true_edge = reference$true_edge,
  true_weight = reference$true_weight,
  stringsAsFactors = FALSE
)
utils::write.csv(truth, out_path("threshold-study-truth.csv"), row.names = FALSE)

cat(sprintf("\nreference edges: %d; gradeable %d (%.1f%%), within 4 SE of a boundary %d (%.1f%%)\n",
            nrow(truth), sum(truth$gradeable), 100 * mean(truth$gradeable),
            sum(!truth$gradeable), 100 * mean(!truth$gradeable)))
cat("reference verdicts:\n")
print(table(truth$truth_verdict))

# Agreement between the reference posterior and the generating graph is a
# property of the data, not of the sampler; it is reported so the reader can
# see the study grades sampling error, not model recovery.
cat("\nreference verdict against the generating edge (data, not sampler):\n")
print(table(verdict = truth$truth_verdict, true_edge = truth$true_edge))

graded <- merge(test, truth, by = c("dataset_id", "edge"))
graded <- graded[graded$gradeable, ]
graded$wrong <- graded$verdict != graded$truth_verdict
cat(sprintf("\ngraded test edge-fits: %d; raw verdict error rate %.3f%%\n",
            nrow(graded), 100 * mean(graded$wrong)))

# --- Rule machinery -------------------------------------------------------------

# A gate is a predicate on the observables. NA means the diagnostic could not
# be computed, which is a failure to certify, not a pass: that is the whole
# point of the settled-edge case.
gate <- function(x) !is.na(x) & x

# Wilson upper bound: with error rates near zero the normal interval is
# useless, and the honest claim "this rule controls the error at 5%" has to
# survive sampling uncertainty in the estimate itself.
wilson_upper <- function(successes, trials, level = 0.95) {
  if (trials == 0) return(NA_real_)
  z <- stats::qnorm(level)
  phat <- successes / trials
  centre <- phat + z^2 / (2 * trials)
  half <- z * sqrt(phat * (1 - phat) / trials + z^2 / (4 * trials^2))
  (centre + half) / (1 + z^2 / trials)
}

evaluate_rules <- function(data, rules, wrong_column = "wrong") {
  wrong <- data[[wrong_column]]
  out <- do.call(rbind, lapply(names(rules), function(name) {
    pass <- rules[[name]](data)
    n_pass <- sum(pass)
    errors <- sum(wrong[pass])
    data.frame(rule = name, n_total = length(pass), n_pass = n_pass,
               pass_rate = n_pass / length(pass), errors = errors,
               error_rate = if (n_pass) errors / n_pass else NA_real_,
               error_upper95 = wilson_upper(errors, n_pass),
               stringsAsFactors = FALSE)
  }))
  out[order(-out$pass_rate), ]
}

frontier <- function(rules_table, target) {
  ok <- rules_table[!is.na(rules_table$error_upper95) &
                      rules_table$error_upper95 <= target, ]
  if (!nrow(ok)) return(ok)
  ok[order(-ok$pass_rate), ][seq_len(min(8L, nrow(ok))), ]
}

# --- Candidate rules: the indicator class ---------------------------------------

indicator_rules <- list()

for (threshold in c(1.005, 1.01, 1.02, 1.05, 1.1)) {
  local({
    t_local <- threshold
    indicator_rules[[sprintf("rhat_rb < %.3f", t_local)]] <<-
      function(d) gate(d$rhat_rb < t_local)
    indicator_rules[[sprintf("rhat_rb (unmasked) < %.3f", t_local)]] <<-
      function(d) gate(d$rhat_rb_unmasked < t_local)
  })
}
for (threshold in c(100, 200, 400, 1000, 2000, 4000)) {
  local({
    t_local <- threshold
    indicator_rules[[sprintf("rb_n_eff >= %d", t_local)]] <<-
      function(d) gate(d$rb_n_eff >= t_local)
    indicator_rules[[sprintf("n_eff_mixt >= %d", t_local)]] <<-
      function(d) gate(d$n_eff_mixt >= t_local)
  })
}
for (threshold in c(1, 10, 50, 100, 500, 2000)) {
  local({
    t_local <- threshold
    indicator_rules[[sprintf("total_flips >= %d", t_local)]] <<-
      function(d) gate(d$total_flips >= t_local)
  })
}
indicator_rules[["not fragile"]] <- function(d) gate(!d$fragile)
indicator_rules[["chains agree on verdict"]] <-
  function(d) gate(d$chains_verdict_agree)
for (alpha in c(0.01, 0.05)) {
  local({
    a_local <- alpha
    indicator_rules[[sprintf("agreement p_boot > %.2f", a_local)]] <<-
      function(d) gate(d$agreement_p_boot > a_local)
  })
}
for (threshold in c(0.02, 0.05, 0.10, 0.20)) {
  local({
    t_local <- threshold
    indicator_rules[[sprintf("pip chain range <= %.2f", t_local)]] <<-
      function(d) gate(d$pip_chain_range <= t_local)
  })
}
# Baseline any diagnostic has to beat: distance of the verdict from the
# boundary, which uses no convergence information at all.
for (threshold in c(0.0, 0.3, 0.5, 1.0)) {
  local({
    t_local <- threshold
    indicator_rules[[sprintf("|log10 BF| - 1 >= %.1f", t_local)]] <<-
      function(d) gate(abs(d$log10_bf) - 1 >= t_local)
  })
}
indicator_rules[["accept everything"]] <- function(d) rep(TRUE, nrow(d))

# Small trees. The manuscript reads the ESS PAIR and lets cross-chain
# agreement arbitrate when the transition ESS is small or missing, so that
# shape is tested directly rather than only its parts.
for (ess in c(200, 400, 1000)) {
  for (threshold in c(1.01, 1.05)) {
    local({
      e_local <- ess
      t_local <- threshold
      indicator_rules[[sprintf(
        "TREE pair: n_eff_mixt >= %d & rhat_rb < %.2f, else chains agree",
        e_local, t_local)]] <<- function(d) {
          computable <- !is.na(d$n_eff_mixt) & !is.na(d$rhat_rb)
          ifelse(computable,
                 gate(d$n_eff_mixt >= e_local) & gate(d$rhat_rb < t_local),
                 gate(d$chains_verdict_agree))
        }
      indicator_rules[[sprintf(
        "TREE pair: n_eff_mixt >= %d & rhat_rb < %.2f, else agreement test",
        e_local, t_local)]] <<- function(d) {
          computable <- !is.na(d$n_eff_mixt) & !is.na(d$rhat_rb)
          ifelse(computable,
                 gate(d$n_eff_mixt >= e_local) & gate(d$rhat_rb < t_local),
                 gate(d$agreement_p_boot > 0.01))
        }
    })
  }
}
# The fragility flag has to be read together with the flip count, because the
# standard error behind it is itself unreliable when flips are scarce.
for (flips in c(50, 100, 500)) {
  local({
    f_local <- flips
    indicator_rules[[sprintf(
      "TREE fragile+flips: not fragile & (flips >= %d or settled)",
      f_local)]] <<- function(d) {
        gate(!d$fragile) & (gate(d$total_flips >= f_local) |
                              gate(d$total_flips == 0))
      }
  })
}
indicator_rules[["TREE settled: chains agree if settled, else not fragile"]] <-
  function(d) {
    ifelse(d$total_flips == 0, gate(d$chains_verdict_agree),
           gate(!d$fragile))
  }

indicator_table <- evaluate_rules(graded, indicator_rules)
utils::write.csv(indicator_table,
                 out_path("threshold-study-rules-indicator.csv"),
                 row.names = FALSE)

cat("\n=== INDICATOR CLASS: frontier at 5% conditional verdict error ===\n")
print(frontier(indicator_table, 0.05), row.names = FALSE, digits = 3)
cat("\n=== INDICATOR CLASS: frontier at 1% ===\n")
print(frontier(indicator_table, 0.01), row.names = FALSE, digits = 3)

# Error catching, as a second criterion next to the conditional error rate.
# When the unconditional verdict error rate is already below a target, every
# rule sits on that frontier and "accept everything" wins it on pass rate, so
# the frontier alone cannot separate the rules. For each rule this computes
# the share of wrong verdicts it rejects (recall), the share of correct
# verdicts it rejects (cost), and the share of its rejections that are wrong
# verdicts (precision).
cat("\n=== INDICATOR CLASS: error catching ===\n")
catch_table <- do.call(rbind, lapply(names(indicator_rules), function(name) {
  pass <- indicator_rules[[name]](graded)
  wrong <- graded$wrong
  data.frame(rule = name,
             wrong_total = sum(wrong),
             caught = sum(!pass & wrong),
             recall = sum(!pass & wrong) / sum(wrong),
             cost = sum(!pass & !wrong) / sum(!wrong),
             precision = if (sum(!pass)) sum(!pass & wrong) / sum(!pass)
                         else NA_real_,
             stringsAsFactors = FALSE)
}))
catch_table <- catch_table[order(-catch_table$recall, catch_table$cost), ]
utils::write.csv(catch_table, out_path("threshold-study-catching.csv"),
                 row.names = FALSE)
print(utils::head(catch_table, 12), row.names = FALSE, digits = 3)

# Where the errors live. Verdict errors concentrate near the boundary by
# construction, so a rate pooled over all edges understates what a reader
# faces on the edges they are actually unsure about.
graded$boundary_distance <- pmin(abs(graded$truth_log10_bf - 1),
                                 abs(graded$truth_log10_bf + 1))
graded$boundary_band <- cut(graded$boundary_distance,
                            breaks = c(0, 0.1, 0.25, 0.5, 1, Inf),
                            labels = c("<0.1", "0.1-0.25", "0.25-0.5",
                                       "0.5-1", ">1"))
band_table <- do.call(rbind, lapply(levels(graded$boundary_band), function(b) {
  block <- graded[graded$boundary_band == b, ]
  data.frame(band = b, edge_fits = nrow(block),
             error_rate = mean(block$wrong),
             errors = sum(block$wrong), stringsAsFactors = FALSE)
}))
cat("\n=== where verdict errors live (distance of the TRUE log10 BF from a boundary) ===\n")
print(band_table, row.names = FALSE, digits = 3)

# Re-run the rule search on the near-boundary edges alone: this is the
# subpopulation a threshold is for.
near_boundary <- graded[graded$boundary_distance <= 0.25, ]
if (nrow(near_boundary) > 50) {
  near_table <- evaluate_rules(near_boundary, indicator_rules)
  utils::write.csv(near_table, out_path("threshold-study-rules-nearboundary.csv"),
                   row.names = FALSE)
  cat(sprintf(
    "\n=== INDICATOR CLASS, near-boundary edges only (n = %d, base error %.2f%%) ===\n",
    nrow(near_boundary), 100 * mean(near_boundary$wrong)))
  print(utils::head(near_table[order(near_table$error_rate,
                                     -near_table$pass_rate), ], 10),
        row.names = FALSE, digits = 3)
}

# --- The mixed-weight and continuous classes ------------------------------------

# Both classes report a number, so "wrong" has to mean "off by enough to
# matter". The criterion is Monte Carlo error exceeding a tenth of a
# posterior standard deviation, the usual precision target and the one the
# manuscript's ESS discussion speaks in (MCSE <= 0.1 SD is ESS >= 100). This
# is not circular with the ESS rules being tested: the question is precisely
# whether the REPORTED effective sample size delivers the accuracy it claims.
#
# A +/- 2 MCSE coverage criterion was tried first and discarded: under a
# well-calibrated MCSE it fails 5% of the time by construction, so it cannot
# separate a good rule from a perfect one at a 5% target. It is kept as a
# secondary column because miscalibration of the MCSE is itself worth seeing.
accuracy_target <- 0.1

weight_graded <- merge(test, truth, by = c("dataset_id", "edge"))
weight_graded$absolute_error <- abs(weight_graded$weight_mean -
                                      weight_graded$truth_weight)
weight_graded$wrong <- weight_graded$absolute_error >
  accuracy_target * weight_graded$weight_sd
weight_graded$mcse_miscovered <- weight_graded$absolute_error >
  2 * weight_graded$weight_mcse

weight_rules <- list()
for (threshold in c(1.005, 1.01, 1.02, 1.05)) {
  local({
    t_local <- threshold
    weight_rules[[sprintf("weight rhat < %.3f", t_local)]] <<-
      function(d) gate(d$weight_rhat < t_local)
  })
}
for (threshold in c(100, 400, 1000, 4000)) {
  local({
    t_local <- threshold
    weight_rules[[sprintf("weight n_eff >= %d", t_local)]] <<-
      function(d) gate(d$weight_n_eff >= t_local)
  })
}
for (threshold in c(0.25, 0.5, 0.75)) {
  local({
    t_local <- threshold
    weight_rules[[sprintf("share_incl <= %.2f", t_local)]] <<-
      function(d) gate(d$weight_share_incl <= t_local)
  })
}
weight_rules[["TREE weight: rhat < 1.01 & n_eff >= 400"]] <-
  function(d) gate(d$weight_rhat < 1.01) & gate(d$weight_n_eff >= 400)
weight_rules[["TREE weight: rhat < 1.01 & n_eff >= 400 & share_incl <= 0.5"]] <-
  function(d) gate(d$weight_rhat < 1.01) & gate(d$weight_n_eff >= 400) &
    gate(d$weight_share_incl <= 0.5)
weight_rules[["accept everything"]] <- function(d) rep(TRUE, nrow(d))

weight_table <- evaluate_rules(weight_graded, weight_rules)
weight_miscoverage <- evaluate_rules(weight_graded, weight_rules,
                                     wrong_column = "mcse_miscovered")
weight_table$mcse_miscoverage <- weight_miscoverage$error_rate[
  match(weight_table$rule, weight_miscoverage$rule)]
utils::write.csv(weight_table, out_path("threshold-study-rules-weight.csv"),
                 row.names = FALSE)
cat("\n=== MIXED-WEIGHT CLASS: frontier at 5% (MC error > 0.1 posterior SD) ===\n")
print(frontier(weight_table, 0.05), row.names = FALSE, digits = 3)

continuous <- merge(
  test_main,
  data.frame(dataset_id = reference_main$dataset_id,
             parameter = reference_main$parameter,
             truth_mean = reference_main$mean, stringsAsFactors = FALSE),
  by = c("dataset_id", "parameter")
)
continuous$absolute_error <- abs(continuous$mean - continuous$truth_mean)
continuous$wrong <- continuous$absolute_error > accuracy_target * continuous$sd
continuous$mcse_miscovered <- continuous$absolute_error > 2 * continuous$mcse

continuous_rules <- list()
for (threshold in c(1.005, 1.01, 1.02, 1.05)) {
  local({
    t_local <- threshold
    continuous_rules[[sprintf("rhat < %.3f", t_local)]] <<-
      function(d) gate(d$rhat < t_local)
  })
}
for (threshold in c(100, 400, 1000, 2000)) {
  local({
    t_local <- threshold
    continuous_rules[[sprintf("n_eff >= %d", t_local)]] <<-
      function(d) gate(d$n_eff >= t_local)
  })
}
for (threshold in c(1.01, 1.02)) {
  for (ess in c(400, 1000)) {
    local({
      t_local <- threshold
      e_local <- ess
      continuous_rules[[sprintf("TREE: rhat < %.2f & n_eff >= %d", t_local,
                                e_local)]] <<-
        function(d) gate(d$rhat < t_local) & gate(d$n_eff >= e_local)
    })
  }
}
continuous_rules[["accept everything"]] <- function(d) rep(TRUE, nrow(d))

continuous_table <- evaluate_rules(continuous, continuous_rules)
continuous_miscoverage <- evaluate_rules(continuous, continuous_rules,
                                         wrong_column = "mcse_miscovered")
continuous_table$mcse_miscoverage <- continuous_miscoverage$error_rate[
  match(continuous_table$rule, continuous_miscoverage$rule)]
utils::write.csv(continuous_table,
                 out_path("threshold-study-rules-continuous.csv"),
                 row.names = FALSE)
cat("\n=== CONTINUOUS CLASS: frontier at 5% (MC error > 0.1 posterior SD) ===\n")
print(frontier(continuous_table, 0.05), row.names = FALSE, digits = 3)
cat("\nby block (main = category thresholds, quadratic = GGM precisions):\n")
print(table(block = continuous$block, wrong = continuous$wrong))

# --- Stability of the chosen rules across density, n, and setting ---------------

# The rule to track is not the frontier winner (which, when the target does
# not bind, is just "accept everything") but the one that catches the errors:
# highest recall, ties broken by the smallest share of correct verdicts
# rejected.
stability <- NULL
best_catching <- catch_table[catch_table$rule != "accept everything", ]
if (nrow(best_catching)) {
  chosen <- best_catching$rule[1]
  cat(sprintf("\n=== STABILITY of \"%s\" ===\n", chosen))
  pass <- indicator_rules[[chosen]](graded)
  cells <- unique(graded[, c("density", "n", "setting")])
  for (k in seq_len(nrow(cells))) {
    in_cell <- graded$density == cells$density[k] & graded$n == cells$n[k] &
      graded$setting == cells$setting[k]
    n_pass <- sum(pass & in_cell)
    errors <- sum(graded$wrong[pass & in_cell])
    stability <- rbind(stability, data.frame(
      rule = chosen, density = cells$density[k], n = cells$n[k],
      setting = cells$setting[k], n_total = sum(in_cell), n_pass = n_pass,
      pass_rate = n_pass / sum(in_cell),
      error_rate = if (n_pass) errors / n_pass else NA_real_,
      error_upper95 = wilson_upper(errors, n_pass), stringsAsFactors = FALSE
    ))
  }
  stability <- stability[order(stability$setting, stability$n,
                               stability$density), ]
  print(stability, row.names = FALSE, digits = 3)
  utils::write.csv(stability, out_path("threshold-study-stability.csv"),
                   row.names = FALSE)
}

# --- Side product: is the indicator chain first-order Markov? -------------------

# The fragile-verdict flag's standard error comes from a two-state Markov
# model, so it assumes the indicator chain is first order. Two-state
# simulations assume it by construction; a real sampler need not obey it.
# Under first-order Markov, P(x_t = 1 | x_{t-1}) does not depend on x_{t-2}.
markov <- test[test$n_10 > 30 & test$n_00 > 30 & test$n_11 > 30 &
                 test$n_01 > 30, ]
markov$p1_after_00 <- markov$n_001 / markov$n_00
markov$p1_after_10 <- markov$n_101 / markov$n_10
markov$p1_after_01 <- markov$n_011 / markov$n_01
markov$p1_after_11 <- markov$n_111 / markov$n_11
markov$lag2_gap_from0 <- markov$p1_after_10 - markov$p1_after_00
markov$lag2_gap_from1 <- markov$p1_after_11 - markov$p1_after_01

# Score test on the pooled 2x2: does the lag-2 state shift the transition
# probability by more than sampling noise?
markov$z_from0 <- markov$lag2_gap_from0 /
  sqrt(markov$p1_after_00 * (1 - markov$p1_after_00) / markov$n_00 +
         markov$p1_after_10 * (1 - markov$p1_after_10) / markov$n_10)
markov$z_from1 <- markov$lag2_gap_from1 /
  sqrt(markov$p1_after_01 * (1 - markov$p1_after_01) / markov$n_01 +
         markov$p1_after_11 * (1 - markov$p1_after_11) / markov$n_11)

markov_summary <- data.frame(
  edge_fits = nrow(markov),
  median_gap_from0 = stats::median(markov$lag2_gap_from0, na.rm = TRUE),
  median_gap_from1 = stats::median(markov$lag2_gap_from1, na.rm = TRUE),
  q95_abs_gap_from0 = stats::quantile(abs(markov$lag2_gap_from0), 0.95,
                                      na.rm = TRUE),
  q95_abs_gap_from1 = stats::quantile(abs(markov$lag2_gap_from1), 0.95,
                                      na.rm = TRUE),
  reject_from0_at_1pct = mean(abs(markov$z_from0) >
                                stats::qnorm(0.995), na.rm = TRUE),
  reject_from1_at_1pct = mean(abs(markov$z_from1) >
                                stats::qnorm(0.995), na.rm = TRUE)
)
utils::write.csv(markov_summary, out_path("threshold-study-markov.csv"),
                 row.names = FALSE)
cat("\n=== FIRST-ORDER MARKOV CHECK on real indicator chains ===\n")
print(markov_summary, row.names = FALSE, digits = 3)
cat("(rejection rates far above 1% mean the two-state standard error behind",
    "the fragile flag rests on an assumption these chains do not satisfy)\n")

# --- Secondary Bayes factor thresholds, descriptive only ------------------------

sensitivity <- do.call(rbind, lapply(list(c(3, 1 / 3), c(30, 1 / 30)),
                                     function(bounds) {
  upper <- log10(bounds[1])
  lower <- log10(bounds[2])
  test_verdict <- verdict_from_log10_bf(graded$log10_bf, upper, lower)
  truth_verdict <- verdict_from_log10_bf(graded$truth_log10_bf, upper, lower)
  data.frame(bf_upper = bounds[1], bf_lower = bounds[2],
             edge_fits = length(test_verdict),
             disagreement = mean(test_verdict != truth_verdict),
             stringsAsFactors = FALSE)
}))
sensitivity <- rbind(
  data.frame(bf_upper = 10, bf_lower = 0.1, edge_fits = nrow(graded),
             disagreement = mean(graded$wrong), stringsAsFactors = FALSE),
  sensitivity
)
utils::write.csv(sensitivity, out_path("threshold-study-bf-sensitivity.csv"),
                 row.names = FALSE)
cat("\n=== Bayes factor threshold sensitivity (descriptive) ===\n")
print(sensitivity, row.names = FALSE, digits = 3)

cat("\nstage", target_stage, "analysis complete\n")
