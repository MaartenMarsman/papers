# GGM calibration pass, step 1: candidate rules on the fits already in hand.
# Design: the simulation appendix of the manuscript (follow-up stage).
#
# The stage-1 fragile-verdict flag need not transfer to the GGM: it rests on
# the two-state summary of the INDICATOR chain, but in a GGM the indicator
# couples to a continuous precision entry, so indicator mixing is not the
# binding constraint and a statistic built on it alone cannot see the error.
# Rather than force it, this evaluates GGM-specific candidates, cheapest
# first:
#
#   1. BOUNDARY BAND. Flag any verdict whose own log10 Bayes factor lies
#      within delta of a threshold. Uses no convergence information at all;
#      the candidate widths bracket the near-boundary band in which the
#      earlier stages locate the verdict errors.
#   2. RB FRAGILITY. The two-MCSE rule the tutorial already teaches, with the
#      standard error taken from the CONTINUOUS Rao-Blackwellized draws (the
#      indicator block's mcse) rather than from the two-state indicator
#      summary, carried to the log10 Bayes factor scale by the delta method.
#   3. Combinations, and the raw precision/ESS columns for comparison.
#
# This step spends no new fits: it runs on the GGM test fits and reference
# runs already checkpointed (with inclusion Bayes factors put on one scale by
# R/26-recompute-inclusion-bf.R). Its purpose is to find which candidates are
# worth spending a targeted grid on.
#
# Output: tables/threshold-study-ggm-rules.csv

source("R/threshold-study-helpers.R")

paths <- study_paths()
target_arm <- "ggm"

manifest <- utils::read.csv(file.path(paths$root, "manifest.csv"),
                            stringsAsFactors = FALSE)
manifest <- manifest[manifest$arm == target_arm &
                       vapply(manifest$fit_id, fit_is_done, logical(1)), ]

read_block <- function(fit_ids, block) {
  do.call(rbind, lapply(fit_ids, function(id) {
    readRDS(fit_checkpoint_path(id))[[block]]
  }))
}

reference <- read_block(manifest$fit_id[manifest$kind == "reference"],
                        "observables")
test <- read_block(manifest$fit_id[manifest$kind == "test"], "observables")

reference$se_log10 <- reference$se_logit_smoothed / log(10)
reference$distance <- pmin(abs(reference$log10_bf - 1),
                           abs(reference$log10_bf + 1))
truth <- data.frame(dataset_id = reference$dataset_id, edge = reference$edge,
                    truth_verdict = reference$verdict,
                    truth_log10_bf = reference$log10_bf,
                    gradeable = reference$distance > 4 * reference$se_log10,
                    stringsAsFactors = FALSE)

graded <- merge(test, truth, by = c("dataset_id", "edge"))
graded <- graded[graded$gradeable, ]
graded$wrong <- graded$verdict != graded$truth_verdict

# Distance of the OBSERVED verdict from a threshold: what a reader can see.
graded$observed_distance <- pmin(abs(graded$log10_bf - 1),
                                 abs(graded$log10_bf + 1))

# Standard error of the log10 Bayes factor from the continuous RB draws.
# summary()'s indicator mcse is the Monte Carlo standard error of the RB
# inclusion probability; on the logit (equivalently log Bayes factor at a
# prior of 1/2) scale the delta method gives mcse / (p (1 - p)).
graded$se_log10_rb <- graded$pip_mcse /
  (graded$pip * (1 - graded$pip) * log(10))

cat(sprintf("GGM: %d graded edge-fits, %d verdict errors (%.3f%%)\n",
            nrow(graded), sum(graded$wrong), 100 * mean(graded$wrong)))
cat(sprintf("RB standard error computable for %d of them (%.1f%%); the rest are settled edges whose mcse is masked\n",
            sum(!is.na(graded$se_log10_rb)),
            100 * mean(!is.na(graded$se_log10_rb))))

gate <- function(x) !is.na(x) & x

rules <- list()
for (delta in c(0.10, 0.15, 0.25, 0.40)) {
  local({
    d_local <- delta
    rules[[sprintf("BAND: |log10 BF| outside %.2f of a threshold", d_local)]] <<-
      function(d) gate(d$observed_distance > d_local)
  })
}
for (k in c(1, 2, 3)) {
  local({
    k_local <- k
    rules[[sprintf("RB FRAGILITY: threshold beyond %d RB standard errors", k_local)]] <<-
      function(d) gate(d$observed_distance > k_local * d$se_log10_rb)
  })
}
# The masked (settled) edges have no RB standard error; route them the way the
# stage-1 tree does, on cross-chain agreement, instead of failing them all.
for (k in c(2, 3)) {
  local({
    k_local <- k
    rules[[sprintf("TREE RB: %d RB SE if computable, else chains agree", k_local)]] <<-
      function(d) {
        ifelse(is.na(d$se_log10_rb), gate(d$chains_verdict_agree),
               gate(d$observed_distance > k_local * d$se_log10_rb))
      }
  })
}
rules[["BAND 0.25 and chains agree"]] <- function(d) {
  gate(d$observed_distance > 0.25) & gate(d$chains_verdict_agree)
}
rules[["TREE RB 2 SE, else chains agree, and BAND 0.15"]] <- function(d) {
  branch <- ifelse(is.na(d$se_log10_rb), gate(d$chains_verdict_agree),
                   gate(d$observed_distance > 2 * d$se_log10_rb))
  branch & gate(d$observed_distance > 0.15)
}
rules[["stage-1 flag: not fragile (two-state)"]] <- function(d) gate(!d$fragile)
rules[["chains agree on verdict"]] <- function(d) gate(d$chains_verdict_agree)
for (threshold in c(1000, 4000)) {
  local({
    t_local <- threshold
    rules[[sprintf("RB n_eff >= %d", t_local)]] <<-
      function(d) gate(d$rb_n_eff >= t_local)
    rules[[sprintf("weight n_eff >= %d", t_local)]] <<-
      function(d) gate(d$weight_n_eff >= t_local)
  })
}
rules[["accept everything"]] <- function(d) rep(TRUE, nrow(d))

wilson_upper <- function(successes, trials, level = 0.95) {
  if (trials == 0) return(NA_real_)
  z <- stats::qnorm(level)
  phat <- successes / trials
  centre <- phat + z^2 / (2 * trials)
  half <- z * sqrt(phat * (1 - phat) / trials + z^2 / (4 * trials^2))
  (centre + half) / (1 + z^2 / trials)
}

evaluate <- function(data, label) {
  do.call(rbind, lapply(names(rules), function(name) {
    pass <- rules[[name]](data)
    wrong <- data$wrong
    data.frame(subset = label, rule = name, n_total = length(pass),
               n_pass = sum(pass), pass_rate = mean(pass),
               errors_after = sum(wrong[pass]),
               error_rate = if (sum(pass)) sum(wrong[pass]) / sum(pass)
                            else NA_real_,
               error_upper95 = wilson_upper(sum(wrong[pass]), sum(pass)),
               recall = if (sum(wrong)) sum(!pass & wrong) / sum(wrong)
                        else NA_real_,
               cost = sum(!pass & !wrong) / sum(!wrong),
               stringsAsFactors = FALSE)
  }))
}

all_rules <- evaluate(graded, "all gradeable")
near <- graded[graded$observed_distance <= 0.5, ]
near_rules <- if (nrow(near) > 50) evaluate(near, "near boundary") else NULL

results <- rbind(all_rules, near_rules)
utils::write.csv(results, "tables/threshold-study-ggm-rules.csv",
                 row.names = FALSE)

cat("\n=== GGM candidates, all gradeable edge-fits (ranked by recall, then cost) ===\n")
ranked <- all_rules[order(-all_rules$recall, all_rules$cost), ]
print(utils::head(ranked[, c("rule", "pass_rate", "recall", "cost",
                             "errors_after", "error_rate")], 14),
      row.names = FALSE, digits = 3)

if (!is.null(near_rules)) {
  cat(sprintf("\n=== near-boundary subset (observed distance <= 0.5, n = %d, base error %.2f%%) ===\n",
              nrow(near), 100 * mean(near$wrong)))
  ranked_near <- near_rules[order(-near_rules$recall, near_rules$cost), ]
  print(utils::head(ranked_near[, c("rule", "pass_rate", "recall", "cost",
                                    "errors_after")], 10),
        row.names = FALSE, digits = 3)
}
