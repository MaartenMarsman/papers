# R/scripts/Z_prior_inspect.R
# Print the full brute-force Z_prior table from R/results/Z_prior_brute_q234.rds.
# Used to inspect non-empty graphs after the main sweep finishes.

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
brute <- readRDS(file.path(PROJECT_ROOT, "R", "results", "Z_prior_brute_q234.rds"))

show_cols <- c("graph", "variant", "delta", "epsilon", "z_hat", "z_se", "log_z",
               "n_accept", "rel_se")

for (var in c("discrete", "continuous")) {
  cat(sprintf("\n=== %s ===\n", var))
  sub <- brute[brute$variant == var, ]
  sub <- sub[order(sub$graph, sub$delta, sub$epsilon), ]
  print(sub[, show_cols], row.names = FALSE, digits = 5)
}

cat("\n\nLog-Z by (graph, variant) at (delta=0.5, epsilon=0.05):\n")
sub2 <- brute[brute$delta == 0.5 & brute$epsilon == 0.05,
              c("graph", "variant", "log_z", "z_se")]
print(sub2[order(sub2$graph, sub2$variant), ], row.names = FALSE, digits = 5)

cat("\n\nAcceptance rates by graph and variant (epsilon=0.1, delta=0):\n")
sub3 <- brute[brute$delta == 0 & brute$epsilon == 0.1,
              c("graph", "variant", "n_accept", "n_draws")]
sub3$acc_rate <- sub3$n_accept / sub3$n_draws
print(sub3[order(sub3$graph, sub3$variant), ], row.names = FALSE, digits = 5)
