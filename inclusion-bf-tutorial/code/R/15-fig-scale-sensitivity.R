# 15-fig-scale-sensitivity.R
# Sensitivity figure for the manuscript's prior_sensitivity_check() passage:
# the package's own plot method for the check object, so the manuscript shows
# exactly what a reader of the tutorial will see on their own screen.
# Continuous log inclusion-Bayes-factor curves over 0.4x-2.5x of the chosen
# slab scale (bgms anchored reweighting engine, 0.2.0.0), with anchor dots,
# the undecided band, and named movers. bgms reports the curve on the natural
# log scale, so the axis needs no conversion here.
#
# Reads tables/prior-sensitivity-check.rds (saved check object); rerun
# R/16-run-sensitivity-check.R if the rds is absent.

set.seed(123)

suppressMessages(library(bgms))

ps <- readRDS("tables/prior-sensitivity-check.rds")

grDevices::pdf("figures/fig-scale-sensitivity.pdf", width = 8, height = 4.8)
plot(ps)
grDevices::dev.off()

edges <- ps$edges
full_range <- edges$stability_lower <= min(ps$multipliers) &
  edges$stability_upper >= max(ps$multipliers)
cat("verdicts holding across the whole range:",
    sum(full_range, na.rm = TRUE), "of", nrow(edges), "\n")
cat("genuine movers:",
    sum(edges$mover == "moved-beyond-wobble" & !edges$insufficient),
    "of", nrow(edges), "\n")
