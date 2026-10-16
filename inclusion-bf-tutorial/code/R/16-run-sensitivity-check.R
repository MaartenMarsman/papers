# 16-run-sensitivity-check.R
# Fits the worked-example model (settings identical to R/01-worked-example.R)
# and runs the anchored prior_sensitivity_check() from bgms (CRAN bgms
# >= 0.2.0.0): the
# continuous inclusion-Bayes-factor curve over 0.4x-2.5x of the chosen slab
# scale, stitched by importance reweighting between warm-started anchor
# refits. The 1x anchor is the original fit, so chosen-scale verdicts equal
# the reported analysis exactly.
#
# Outputs: tables/prior-sensitivity-check.rds  (full check object)
#          tables/prior-sensitivity-edges.csv  (per-edge table)

set.seed(123)

data("Wenchuan", package = "bgms")
wenchuan <- Wenchuan[stats::complete.cases(Wenchuan), ]

t0 <- Sys.time()
fit <- bgms::bgm(
  x = wenchuan,
  variable_type = "ordinal",
  iter = 1e4,
  warmup = 5e3,
  edge_selection = TRUE,
  edge_prior = bgms::bernoulli_prior(0.5),
  interaction_prior = bgms::normal_prior(scale = 1),
  chains = 4,
  seed = 123
)
fit_seconds <- as.numeric(Sys.time() - t0, units = "secs")
cat("fit:", round(fit_seconds), "s\n")

t0 <- Sys.time()
check <- bgms::prior_sensitivity_check(fit, seed = 123)
check_seconds <- as.numeric(Sys.time() - t0, units = "secs")
cat("check:", round(check_seconds), "s =",
    round(check_seconds / fit_seconds, 2), "x the fit\n")

saveRDS(check, "tables/prior-sensitivity-check.rds")
utils::write.csv(check$edges, "tables/prior-sensitivity-edges.csv",
                 row.names = FALSE)

print(check)

cat("\nmover counts:\n")
print(table(check$edges$mover, useNA = "ifany"))
cat("\nnot certifiable (insufficient):", sum(check$edges$insufficient), "\n")
full_range <- check$edges$stability_lower <= min(check$multipliers) &
  check$edges$stability_upper >= max(check$multipliers)
cat("verdicts holding across the whole range:",
    sum(full_range, na.rm = TRUE), "of", nrow(check$edges), "\n")
cat("data-preferred scale:", signif(check$preferred_scale$s_hat, 3), "\n")
