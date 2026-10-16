# 38-consolidated-classification.R
# Folds every check of the tutorial into one per-edge classification
# (manuscript Section: The Robustness of the Verdicts, consolidated table).
#
# Motivation: read separately, the checks fragment the
# analysis into overlapping soft categories (fragile, scale mover, scale
# noisy, structure mover); read together they measure one condition,
# distance to an evidence threshold. This script computes the combined
# classification and the distance summaries quoted in the manuscript.
#
# Inputs:  tables/edge-supplement.csv          (verdict, fragile, log_bf)
#          tables/prior-sensitivity-edges.csv  (scale-check classes)
#          tables/edge-prior-check.csv         (Beta-Bernoulli verdicts)
# Output:  tables/consolidated-robustness.csv  (one row per edge)
#
# A consolidated robustness classification of this kind may later be
# provided by the bgms package itself, in which case this script becomes a
# check against that output.

sup <- read.csv("tables/edge-supplement.csv")
ps <- read.csv("tables/prior-sensitivity-edges.csv")
bb <- read.csv("tables/edge-prior-check.csv")
bb$edge <- paste(bb$from, bb$to, sep = "-")

m <- merge(sup[, c("edge", "verdict", "fragile", "log_bf")],
           ps[, c("edge", "mover", "insufficient")], by = "edge")
m <- merge(m, bb[, c("edge", "verdict_bb")], by = "edge")
stopifnot(nrow(m) == 136)

m$scale_moved <- m$mover == "moved-beyond-wobble" & !m$insufficient
m$scale_noisy <- m$insufficient
m$structure_moved <- m$verdict != m$verdict_bb
m$flagged <- m$fragile | m$scale_moved | m$scale_noisy | m$structure_moved

m$final <- ifelse(m$verdict == "undecided", "undecided",
           ifelse(m$flagged, "decisive, movable", paste(m$verdict, "settled")))

utils::write.csv(m, "tables/consolidated-robustness.csv", row.names = FALSE)

n_flag <- sum(m$flagged)
n_multi <- sum(m$fragile + m$scale_moved + m$scale_noisy + m$structure_moved > 1)
dec <- m[m$verdict != "undecided", ]
dist <- abs(abs(dec$log_bf) - log(10))
med_soft <- stats::median(dist[dec$flagged])
med_robust <- stats::median(dist[!dec$flagged])

cat("final classification:\n")
print(table(m$final))
cat("flagged edges:", n_flag, " with more than one flag:", n_multi, "\n")
cat("decisive flagged: median distance", round(med_soft, 2), "log BF = factor",
    round(exp(med_soft), 2), "\n")
cat("decisive unflagged: median distance", round(med_robust, 2), "log BF = factor",
    round(exp(med_robust), 2), "\n")
