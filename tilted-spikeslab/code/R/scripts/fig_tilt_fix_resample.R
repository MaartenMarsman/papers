# R/scripts/fig_tilt_fix_resample.R
#
# Figure 2: the tilt's effect on the bulk of the PD-conditional prior, using
# rejection-sampled (not importance-reweighted) draws from the delta-tilted
# prior. Data: R/results/cone_interior_tilt_resample.rds.

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
IN_RDS       <- file.path(PROJECT_ROOT, "R", "results", "cone_interior_tilt_resample.rds")
OUT_PDF      <- file.path(PROJECT_ROOT, "figures", "tilt_fix.pdf")

obj <- readRDS(IN_RDS)
tab <- obj$summary

# Drop cells with too few accepted draws for stable quantile estimation.
min_accepted <- 25L
tab$keep <- tab$n_accepted >= min_accepted
df <- tab[tab$keep, ]
df$variant_label <- factor(df$variant,
                            levels = c("discrete", "continuous"),
                            labels = c("Discrete", "Continuous"))
df$q_factor <- factor(df$q, levels = sort(unique(tab$q)))

theme_paper <- theme_minimal(base_size = 9) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = 8),
    panel.grid.minor   = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 9, face = "plain"),
    plot.title         = element_text(size = 9, face = "plain")
  )

palette_q <- c(`5` = "#a6cee3", `10` = "#1f78b4", `15` = "#b2df8a",
               `20` = "#33a02c", `30` = "#fb9a99", `50` = "#e31a1c")
shapes_q  <- c(`5` = 16, `10` = 17, `15` = 15, `20` = 18, `30` = 4, `50` = 3)

fig <- ggplot(df,
              aes(x = delta, y = median,
                  colour = q_factor, fill = q_factor, shape = q_factor)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.10, colour = NA) +
  geom_line() +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 0, y = 0.05, label = expression(lambda[min] == 0.05),
           hjust = -0.05, vjust = -0.4, size = 2.6, colour = "grey40") +
  facet_wrap(~ variant_label, ncol = 2) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 2)) +
  scale_colour_manual(values = palette_q, name = "q") +
  scale_fill_manual(values = palette_q, name = "q") +
  scale_shape_manual(values = shapes_q, name = "q") +
  labs(x = expression(tilt~exponent~delta),
       y = expression(lambda[min] * "(K) given PD, delta-tilted prior"),
       title = "Tilt-only effect on the bulk of the PD-conditional prior (rejection-sampled)") +
  theme_paper

ggplot2::ggsave(OUT_PDF, plot = fig, width = 9, height = 3.7, units = "in")
cat(sprintf("Wrote %s\n", OUT_PDF))
cat("\nCells included:\n")
print(tab[, c("variant", "q", "delta", "n_proposals", "n_pd", "n_accepted",
              "median", "keep")])
