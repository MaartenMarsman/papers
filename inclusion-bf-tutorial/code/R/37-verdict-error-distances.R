# 37-verdict-error-distances.R
# How far from a verdict boundary did the simulation study's errors sit,
# and did any error ever cross between presence and absence?
#
# For every erroneous verdict in the graded test fits, the script computes
# the distance of the true log Bayes factor (from the reference run on the
# same data set) to the nearest evidence threshold, on the log10, natural-log,
# and Bayes-factor scales, and reports the observed maximum, median, and 90th
# percentile together with the number of presence-absence swaps.
#
# Inputs:  tables/threshold-study/fits/*.rds (per-fit checkpoints; test
#          fits and reference runs; the two "-truncated" pilot reference
#          checkpoints are excluded from truth)
# Outputs: tables/verdict-error-distances.csv (one row per verdict error)
#          and a printed summary quoted in the simulation appendix.
#
# Gradeability is recomputed here from the reference runs: an edge is graded
# when its reference log10 Bayes factor sits more than 4 smoothed-logit SEs
# from both verdict boundaries.

# Iteration count of the reference analyses; see 20-threshold-study-design.R.
# Reference checkpoints with any other count are truncated pilot runs.
reference_iter <- 25000

fit_dir <- file.path(Sys.getenv("THRESHOLD_STUDY_DIR", "results/threshold-study"), "fits")
files <- list.files(fit_dir, full.names = TRUE)
if (length(files) == 0) {
  stop("No simulation-study checkpoints found in ", fit_dir,
       "; see README, section on the simulation study.")
}

all_obs <- lapply(files, function(f) readRDS(f)$observables)
kind1 <- vapply(all_obs, function(o) o$kind[1], character(1))
is_truncated_ref <- vapply(
  all_obs,
  function(o) o$kind[1] == "reference" && o$iter[1] != reference_iter,
  logical(1)
)

test_cols <- c("dataset_id", "arm", "setting", "edge",
               "log10_bf", "verdict", "fragile")
test <- do.call(rbind, lapply(all_obs[kind1 == "test"],
                              function(o) o[, test_cols]))

ref_cols <- c("dataset_id", "edge", "log10_bf", "verdict",
              "se_logit_smoothed")
refs <- do.call(rbind, lapply(all_obs[kind1 == "reference" & !is_truncated_ref],
                              function(o) o[, ref_cols]))

# Gradeability: the reference Bayes factor must sit more than 4 SEs from
# both verdict boundaries (log10 BF = +1 and -1); SE carried from the
# smoothed logit scale to log10 odds.
ref_se_log10 <- refs$se_logit_smoothed / log(10)
ref_dist <- pmin(abs(refs$log10_bf - 1), abs(refs$log10_bf + 1))
refs$gradeable <- is.finite(refs$log10_bf) & is.finite(ref_se_log10) &
  ref_dist > 4 * ref_se_log10

names(refs)[names(refs) == "log10_bf"] <- "truth_log10_bf"
names(refs)[names(refs) == "verdict"] <- "truth_verdict"

graded <- merge(
  test,
  refs[refs$gradeable, c("dataset_id", "edge", "truth_verdict", "truth_log10_bf")],
  by = c("dataset_id", "edge")
)

errors <- graded[graded$verdict != graded$truth_verdict, ]
errors$dist_log10 <- pmin(abs(errors$truth_log10_bf - 1),
                          abs(errors$truth_log10_bf + 1))
errors$dist_ln <- errors$dist_log10 * log(10)
errors$factor_from_threshold <- 10^errors$dist_log10
errors <- errors[order(-errors$dist_log10), ]

utils::write.csv(errors, "tables/verdict-error-distances.csv",
                 row.names = FALSE)

crossed <- sum(
  (errors$truth_verdict == "presence" & errors$verdict == "absence") |
  (errors$truth_verdict == "absence" & errors$verdict == "presence")
)

cat("graded edge-fits:", nrow(graded), "\n")
cat("verdict errors:", nrow(errors), "\n")
print(table(errors$arm))
cat("max distance:", round(max(errors$dist_log10), 3), "log10 =",
    round(max(errors$dist_ln), 3), "natural log = factor",
    round(max(errors$factor_from_threshold), 2), "\n")
cat("median factor from threshold:",
    round(10^stats::median(errors$dist_log10), 2), "\n")
cat("90th percentile factor:",
    round(10^stats::quantile(errors$dist_log10, 0.9), 2), "\n")
cat("presence<->absence swaps:", crossed, "\n")
print(table(truth = errors$truth_verdict, reported = errors$verdict))
