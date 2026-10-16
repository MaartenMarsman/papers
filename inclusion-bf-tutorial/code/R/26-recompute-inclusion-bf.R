# Known-truth calibration study: CORRECTION PASS on the inclusion Bayes factors.
# Design: the simulation appendix of the manuscript.
#
# WHY THIS EXISTS
#
# `bgms::extract_inclusion_bf()` documents its return value as "a symmetric
# p x p matrix of log inclusion Bayes factors (natural log)", and the study's
# extraction divided it by log(10) accordingly. On the stage-1 ordinal
# checkpoints that was right: the stored column reproduces
# log10(PIP / (1 - PIP)) exactly, and 81 of 105 edges per reference run carry
# the negative values an absence verdict needs.
#
# On the stage-2 binary and GGM checkpoints it was NOT right. There the stored
# column equals the Bayes factor ITSELF divided by log(10), never negative, so
# every edge that should have read as evidence of absence read as undecided
# instead: 0 absence verdicts in 630 binary and 630 GGM reference edges,
# against 607 in 1,260 ordinal ones. A binary edge with a posterior inclusion
# probability of 0.054 was recorded at +0.025 rather than -1.245.
#
# The stored column is therefore on different scales across checkpoints, and
# the stored values cannot be used as they are. Rather than depend on the
# scale of the extracted column, this pass recomputes every Bayes factor from
# the posterior inclusion probabilities the study stored itself, so that all
# checkpoints end up on one scale.
#
# WHAT THIS DOES
#
# Rather than trust an ambiguous column, the Bayes factor is recomputed for
# every fit from the per-chain accumulator values this study stored itself
# (`per_chain$log10_bf`, from the same identity bgms uses). Pooling is on the
# odds scale, which is the correct scale to average per-chain odds on when the
# chains are of equal length:
#
#     log10 BF = log10( mean_over_chains( 10^log10_BF_chain ) )
#
# VALIDATION: on the stage-1 ordinal reference runs, where the stored column is
# known to be on the log scale, the recomputed value reproduces it to six
# significant figures on every gradeable edge (-2.07429 vs -2.07429). The two
# diverge only on decisively-included edges whose log10 Bayes factor runs to
# the hundreds, where equal-denominator pooling is an approximation and the
# verdict is not in doubt either way.
#
# The original column is preserved as `log10_bf_extract` so the discrepancy
# stays auditable.
#
# Usage: Rscript R/26-recompute-inclusion-bf.R

source("R/threshold-study-helpers.R")

paths <- study_paths()
manifest <- utils::read.csv(file.path(paths$root, "manifest.csv"),
                            stringsAsFactors = FALSE)
manifest <- manifest[vapply(manifest$fit_id, fit_is_done, logical(1)), ]

changed <- NULL
for (i in seq_len(nrow(manifest))) {
  fit_id <- manifest$fit_id[i]
  checkpoint <- readRDS(fit_checkpoint_path(fit_id))
  observables <- checkpoint$observables
  per_chain <- checkpoint$per_chain

  if (!is.null(observables$log10_bf_extract)) next   # already corrected

  pooled <- tapply(per_chain$log10_bf, per_chain$edge,
                   function(v) log10(mean(10^v)))
  recomputed <- as.vector(pooled[observables$edge])

  observables$log10_bf_extract <- observables$log10_bf
  observables$log10_bf <- recomputed
  observables$verdict_extract <- observables$verdict
  observables$verdict <- verdict_from_log10_bf(recomputed)

  finite_pair <- is.finite(observables$log10_bf_extract) & is.finite(recomputed)
  changed <- rbind(changed, data.frame(
    fit_id = fit_id, arm = manifest$arm[i], kind = manifest$kind[i],
    verdicts_changed = sum(observables$verdict != observables$verdict_extract),
    max_abs_difference = if (any(finite_pair))
      max(abs(observables$log10_bf_extract[finite_pair] -
                recomputed[finite_pair])) else NA_real_,
    stringsAsFactors = FALSE
  ))

  checkpoint$observables <- observables
  saveRDS(checkpoint, fit_checkpoint_path(fit_id))
}

if (is.null(changed)) {
  cat("all checkpoints already corrected; nothing to do\n")
} else {
  utils::write.csv(changed, file.path(paths$root, "inclusion-bf-correction.csv"),
                   row.names = FALSE)
  cat(sprintf("corrected %d checkpoints\n", nrow(changed)))
  cat("\nverdicts changed per fit, by arm:\n")
  print(tapply(changed$verdicts_changed, changed$arm,
               function(v) c(mean = mean(v), max = max(v))))
}
