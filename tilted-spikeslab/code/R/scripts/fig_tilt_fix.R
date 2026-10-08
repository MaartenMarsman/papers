# R/scripts/fig_tilt_fix.R
#
# Figure 2 for the manuscript: empirical demonstration that the determinant
# tilt alone is sufficient to push the bulk of the prior off the cone
# boundary. For each q in {5, 10, 15, 20} and delta in {0, 0.5, 1, 2},
# plot the (weighted) median and 5-95% interval of lambda_min(K) under the
# PD-conditional prior, importance-reweighted by |K|^delta.
#
# Output: figures/tilt_fix.pdf
#
# Data: R/results/cone_interior_diagnostics.rds (built by
# R/scripts/cone_interior_diagnostics.R).

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
IN_RDS       <- file.path(PROJECT_ROOT, "R", "results", "cone_interior_diagnostics.rds")
OUT_PDF      <- file.path(PROJECT_ROOT, "figures", "tilt_fix.pdf")

obj <- readRDS(IN_RDS)
tilt_tab <- obj$tilt_only

# At q=20 continuous the PD sample size is too small (~20 draws); IS at large
# delta degenerates onto a few draws. Drop unreliable cells from the figure.
tilt_tab <- subset(tilt_tab, !(variant == "continuous" & q == 20))

tilt_tab$variant_label <- factor(tilt_tab$variant,
                                  levels = c("discrete", "continuous"),
                                  labels = c("Discrete (bgms-style)",
                                             "Continuous (Wang-style)"))

theme_paper <- theme_minimal(base_size = 10) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = 9),
    panel.grid.minor   = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 10, face = "plain"),
    plot.title         = element_text(size = 10, face = "plain")
  )

palette_q <- c("5" = "#a6cee3", "10" = "#1f78b4", "15" = "#b2df8a", "20" = "#33a02c")

fig <- ggplot(tilt_tab,
              aes(x = delta, y = lmin_median, colour = factor(q), fill = factor(q),
                  shape = factor(q))) +
  geom_ribbon(aes(ymin = lmin_q05, ymax = lmin_q95), alpha = 0.15, colour = NA) +
  geom_line() +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 0, y = 0.05, label = expression(lambda[min] == 0.05),
           hjust = -0.05, vjust = -0.4, size = 2.8, colour = "grey40") +
  facet_wrap(~ variant_label, ncol = 2) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 2)) +
  scale_colour_manual(values = palette_q, name = "q") +
  scale_fill_manual(values = palette_q, name = "q") +
  scale_shape_manual(values = c(16, 17, 15, 18), name = "q") +
  labs(x = expression(tilt~exponent~delta),
       y = expression(lambda[min] * "(K) given PD (tilt-weighted)"),
       title = expression("Tilt-only effect on the bulk of the PD-conditional prior")) +
  theme_paper

ggplot2::ggsave(OUT_PDF, plot = fig, width = 9, height = 3.7, units = "in")
cat(sprintf("Wrote %s\n", OUT_PDF))
