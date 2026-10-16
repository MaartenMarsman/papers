# 20-fig-centrality.R
# Posterior strength centrality for the worked example: the centrality is
# evaluated on each sampled network, so uncertainty about which edges exist
# propagates into the centralities automatically. The manuscript's
# interpretation section shows bgms::extract_centrality(); this script keeps
# the computation explicit so the figure does not depend on the extractor's
# summary format (the two agree).
#
# Reads tables/fit-wenchuan.rds (saved by R/01-worked-example.R).
# Outputs: figures/fig-centrality.pdf, tables/centrality-summary.csv

set.seed(123)

fit <- readRDS("tables/fit-wenchuan.rds")

# --- The four lines from the manuscript ----------------------------------------

samples <- bgms::extract_pairwise_interactions(fit)
pairs <- do.call(rbind, strsplit(colnames(samples), "-", fixed = TRUE))
node_names <- sort(unique(c(pairs)))
# An edge excluded in a draw enters that draw's row sum as a zero weight,
# which is how structure uncertainty propagates into the centralities.
strength <- sapply(node_names, function(node)
  rowSums(abs(samples[, pairs[, 1] == node | pairs[, 2] == node])))

# Each row of `strength` is one posterior draw of all 17 centralities.
# Probability of being the most central node, per node:
p_top <- table(factor(colnames(strength)[max.col(strength)],
                      levels = node_names)) / nrow(strength)

# --- Summary table -------------------------------------------------------------

summary_df <- data.frame(
  node = node_names,
  mean = colMeans(strength),
  lower = apply(strength, 2, stats::quantile, 0.025),
  upper = apply(strength, 2, stats::quantile, 0.975),
  p_most_central = as.numeric(p_top)
)
summary_df <- summary_df[order(summary_df$mean), ]
utils::write.csv(summary_df, "tables/centrality-summary.csv", row.names = FALSE)

cat("posterior strength centrality (top of the ranking last):\n")
print(round(summary_df[, c("mean", "lower", "upper", "p_most_central")], 3))

# --- Figure: posterior means with 95% credible intervals -----------------------

summary_df$node <- factor(summary_df$node, levels = summary_df$node)

fig_centrality <- ggplot2::ggplot(summary_df, ggplot2::aes(x = mean, y = node)) +
  ggplot2::geom_segment(
    ggplot2::aes(x = lower, xend = upper, yend = node),
    color = "#12655F", linewidth = 1.6, alpha = 0.45
  ) +
  ggplot2::geom_point(color = "#12655F", size = 2.2) +
  ggplot2::labs(
    x = "Strength centrality (posterior mean and 95% credible interval)",
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

ggplot2::ggsave(
  filename = "figures/fig-centrality.pdf",
  plot = fig_centrality,
  width = 6.5,
  height = 4.5
)
