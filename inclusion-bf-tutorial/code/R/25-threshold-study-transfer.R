# Known-truth calibration study, step 6: the STAGE 2 TRANSFER CHECK.
# Design: the simulation appendix of the manuscript.
#
# Stage 2 validates the stage-1 decision tree on the other variable types. It
# does NOT refit it: the rules and their thresholds are fixed here exactly as
# stage 1 chose them, and the only question is whether they still hold on the
# binary and GGM arms. Refitting per arm would turn a transfer check into
# three separate calibrations and would leave nothing to transfer.
#
# Only if transfer fails does an arm earn its own calibration pass; this
# script reports such failures and does not recalibrate anything itself.
#
# Output: tables/threshold-study-transfer.csv and
#         tables/threshold-study-pair-rule.csv

source("R/threshold-study-helpers.R")

paths <- study_paths()
accuracy_target <- 0.1     # as in R/23: MC error above a tenth of a posterior SD

# --- The stage-1 rules, frozen ---------------------------------------------------

gate <- function(x) !is.na(x) & x

wilson_upper <- function(successes, trials, level = 0.95) {
  if (trials == 0) return(NA_real_)
  z <- stats::qnorm(level)
  phat <- successes / trials
  centre <- phat + z^2 / (2 * trials)
  half <- z * sqrt(phat * (1 - phat) / trials + z^2 / (4 * trials^2))
  (centre + half) / (1 + z^2 / trials)
}

# Three indicator rules are carried side by side, because the two channels
# have complementary blind spots (see R/28-ggm-rule-search.R): a flag built on
# the indicator chain alone need not transfer to the GGM arm.
#
#   two-state   the stage-1 fragile-verdict flag, built on the smoothed
#               two-state summary of the INDICATOR chain: the exploration
#               channel;
#   rb          the same two-standard-error idea with the standard error
#               taken from the continuous Rao-Blackwellized draws (the
#               indicator block's mcse, carried to the log10 Bayes factor
#               scale by the delta method): the precision channel;
#   pair        both channels must pass.
#
# Settled edges have no RB standard error (summary() masks the mcse at zero
# flips), so the RB branch falls back to cross-chain agreement there.
indicator_rule_two_state <- function(d) gate(!d$fragile)
indicator_rule_rb <- function(d) {
  se_rb <- d$pip_mcse / (d$pip * (1 - d$pip) * log(10))
  distance <- pmin(abs(d$log10_bf - 1), abs(d$log10_bf + 1))
  ifelse(is.na(se_rb), gate(d$chains_verdict_agree),
         gate(distance > 2 * se_rb))
}
indicator_rule <- function(d) indicator_rule_two_state(d) & indicator_rule_rb(d)
weight_rule <- function(d) gate(d$weight_n_eff >= 1000)
continuous_rule <- function(d) gate(d$rhat < 1.01)

# --- Load every stage's checkpoints ---------------------------------------------

manifest <- utils::read.csv(file.path(paths$root, "manifest.csv"),
                            stringsAsFactors = FALSE)
manifest <- manifest[vapply(manifest$fit_id, fit_is_done, logical(1)), ]

read_block <- function(fit_ids, block) {
  do.call(rbind, lapply(fit_ids, function(id) {
    readRDS(fit_checkpoint_path(id))[[block]]
  }))
}

reference <- read_block(manifest$fit_id[manifest$kind == "reference"],
                        "observables")
test <- read_block(manifest$fit_id[manifest$kind == "test"], "observables")
reference_main <- read_block(manifest$fit_id[manifest$kind == "reference"],
                             "main")
test_main <- read_block(manifest$fit_id[manifest$kind == "test"], "main")

reference$se_log10 <- reference$se_logit_smoothed / log(10)
reference$distance_to_boundary <- pmin(abs(reference$log10_bf - 1),
                                       abs(reference$log10_bf + 1))
truth <- data.frame(
  dataset_id = reference$dataset_id, edge = reference$edge,
  truth_verdict = reference$verdict, truth_log10_bf = reference$log10_bf,
  truth_weight = reference$weight_mean,
  gradeable = reference$distance_to_boundary > 4 * reference$se_log10,
  stringsAsFactors = FALSE
)

# --- Evaluate, per arm ------------------------------------------------------------

summarize <- function(pass, wrong, arm, class, rule_label) {
  n_pass <- sum(pass)
  errors <- sum(wrong[pass])
  data.frame(
    arm = arm, class = class, rule = rule_label,
    n_total = length(pass), n_pass = n_pass, pass_rate = n_pass / length(pass),
    base_error_rate = mean(wrong),
    errors_after_rule = errors,
    error_rate = if (n_pass) errors / n_pass else NA_real_,
    error_upper95 = wilson_upper(errors, n_pass),
    recall = if (sum(wrong)) sum(!pass & wrong) / sum(wrong) else NA_real_,
    cost = if (sum(!wrong)) sum(!pass & !wrong) / sum(!wrong) else NA_real_,
    stringsAsFactors = FALSE
  )
}

graded <- merge(test, truth, by = c("dataset_id", "edge"))
graded <- graded[graded$gradeable, ]
graded$wrong <- graded$verdict != graded$truth_verdict

weight_graded <- merge(test, truth, by = c("dataset_id", "edge"))
weight_graded$wrong <- abs(weight_graded$weight_mean -
                             weight_graded$truth_weight) >
  accuracy_target * weight_graded$weight_sd

continuous <- merge(
  test_main,
  data.frame(dataset_id = reference_main$dataset_id,
             parameter = reference_main$parameter,
             truth_mean = reference_main$mean, stringsAsFactors = FALSE),
  by = c("dataset_id", "parameter")
)
continuous$wrong <- abs(continuous$mean - continuous$truth_mean) >
  accuracy_target * continuous$sd

transfer <- NULL
for (arm in sort(unique(test$arm))) {
  block <- graded[graded$arm == arm, ]
  if (nrow(block)) {
    transfer <- rbind(transfer,
      summarize(indicator_rule_two_state(block), block$wrong, arm,
                "indicator", "two-state flag"),
      summarize(indicator_rule_rb(block), block$wrong, arm,
                "indicator", "RB 2 SE, else chains agree"),
      summarize(indicator_rule(block), block$wrong, arm,
                "indicator", "PAIR: both channels"))
  }
  block <- weight_graded[weight_graded$arm == arm, ]
  if (nrow(block)) {
    transfer <- rbind(transfer, summarize(weight_rule(block), block$wrong,
                                          arm, "weight",
                                          "derived ESS >= 1000"))
  }
  block <- continuous[continuous$arm == arm, ]
  if (nrow(block)) {
    transfer <- rbind(transfer, summarize(continuous_rule(block), block$wrong,
                                          arm, "continuous",
                                          "split-Rhat < 1.01"))
  }
}

utils::write.csv(transfer, "tables/threshold-study-transfer.csv",
                 row.names = FALSE)

# Pooled pair rule over all arms; R/24-fig-decision-tree.R reads this file for
# the indicator leaf of the decision tree.
pooled_pair <- summarize(indicator_rule(graded), graded$wrong, "all",
                         "indicator", "PAIR: both channels pass")
utils::write.csv(
  data.frame(rule = pooled_pair$rule, n_total = pooled_pair$n_total,
             n_pass = pooled_pair$n_pass, pass_rate = pooled_pair$pass_rate,
             errors = pooled_pair$errors_after_rule,
             error_rate = pooled_pair$error_rate,
             error_upper95 = pooled_pair$error_upper95,
             stringsAsFactors = FALSE),
  "tables/threshold-study-pair-rule.csv", row.names = FALSE)

cat("=== TRANSFER CHECK: stage-1 rules applied unchanged to every arm ===\n")
cat("(ordinal is the arm the rules were calibrated on; binary and GGM are the test)\n\n")
print(transfer, row.names = FALSE, digits = 3)

cat("\n=== gradeable share and reference verdicts by arm ===\n")
reference_arm <- merge(truth, unique(reference[, c("dataset_id", "arm")]),
                       by = "dataset_id")
print(table(arm = reference_arm$arm, gradeable = reference_arm$gradeable))
print(table(arm = reference_arm$arm, verdict = reference_arm$truth_verdict))

cat("\n=== where verdict errors live, by arm ===\n")
graded$boundary_distance <- pmin(abs(graded$truth_log10_bf - 1),
                                 abs(graded$truth_log10_bf + 1))
graded$band <- cut(graded$boundary_distance, c(-Inf, 0.25, Inf),
                   labels = c("within 0.25 of a boundary", "beyond 0.25"))
print(stats::xtabs(wrong ~ arm + band, graded))
cat("\nedge-fits in each cell:\n")
print(table(arm = graded$arm, band = graded$band))

# Transfer has to be judged on the criterion the rules were CHOSEN on, which
# is error catching, not the 5% target. Every arm's unconditional verdict
# error rate is already far below 5%, so a target test would pass everything
# and report transfer where none happened. A rule transfers if it still
# rejects the errors it was picked for; the ordinal arm is the benchmark it is
# measured against.
recall_floor <- 0.8
cat("\n=== verdict: does each rule still catch what it was chosen to catch? ===\n")
judged <- transfer[!is.na(transfer$recall), ]
judged$transfers <- judged$recall >= recall_floor
print(judged[, c("arm", "class", "rule", "base_error_rate", "recall", "cost",
                 "transfers")], row.names = FALSE, digits = 3)

failures <- judged[!judged$transfers, ]
if (nrow(failures)) {
  cat(sprintf("\nTRANSFER FAILS (recall below %.0f%%) for:\n", 100 * recall_floor))
  for (i in seq_len(nrow(failures))) {
    cat(sprintf("  %s / %s: recall %.2f -- the rule misses %d of %d errors\n",
                failures$arm[i], failures$class[i], failures$recall[i],
                round((1 - failures$recall[i]) *
                        failures$base_error_rate[i] * failures$n_total[i]),
                round(failures$base_error_rate[i] * failures$n_total[i])))
  }
  cat("\nOnly an arm whose transfer fails earns its own calibration pass;\n",
      "this script does not recalibrate.\n", sep = "")
} else {
  cat(sprintf("\nEvery rule keeps recall at or above %.0f%% on every arm.\n",
              100 * recall_floor))
}

# The one finding that should transfer regardless of the rules: errors live
# at the boundary.
cat("\n=== does the boundary localization transfer? ===\n")
localization <- do.call(rbind, lapply(sort(unique(graded$arm)), function(a) {
  block <- graded[graded$arm == a, ]
  data.frame(arm = a, errors = sum(block$wrong),
             errors_within_025 = sum(block$wrong &
                                       block$boundary_distance <= 0.25),
             share_within_025 = if (sum(block$wrong))
               mean(block$boundary_distance[block$wrong] <= 0.25) else NA_real_,
             stringsAsFactors = FALSE)
}))
print(localization, row.names = FALSE, digits = 3)
