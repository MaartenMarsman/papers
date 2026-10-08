# R/scripts/fig_tilt_fix_expanded.R
#
# Figure 2: the tilt's effect on the bulk of the PD-conditional prior.
# Two panels (discrete, continuous), x = delta, y = median lambda_min with
# 5-95 percent ribbon, line / point per q value. Reads from the expanded
# cone_interior dataset and falls back to the old narrower sweep if absent.
#
# Uses sigma = 1 (the worst-case operating regime) for the main figure.

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
IN_RDS       <- file.path(PROJECT_ROOT, "R", "results", "cone_interior_expanded.rds")
OUT_PDF      <- file.path(PROJECT_ROOT, "figures", "tilt_fix.pdf")

obj <- readRDS(IN_RDS)

# Build the tilt-summary long table at sigma = 1
rows <- list()
min_ess <- 50  # IS effective sample size threshold below which we drop the cell
for (r in obj$results) {
  if (abs(r$sigma - 1) > 1e-12) next
  for (ts in r$tilt_summary) {
    if (is.na(ts$median) || is.null(ts$ess) || ts$ess < min_ess) next
    rows[[length(rows) + 1L]] <- data.frame(
      variant = r$variant,
      q       = r$q,
      delta   = ts$delta,
      median  = ts$median,
      q05     = ts$q05,
      q95     = ts$q95,
      ess     = ts$ess,
      stringsAsFactors = FALSE
    )
  }
}
df <- do.call(rbind, rows)
df$variant_label <- factor(df$variant,
                            levels = c("discrete", "continuous"),
                            labels = c("Discrete", "Continuous"))
df$q_factor <- factor(df$q, levels = sort(unique(df$q)))

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
       y = expression(lambda[min] * "(K) given PD (tilt-weighted)"),
       title = "Tilt-only effect on the bulk of the PD-conditional prior") +
  theme_paper

ggplot2::ggsave(OUT_PDF, plot = fig, width = 9, height = 3.7, units = "in")
cat(sprintf("Wrote %s\n", OUT_PDF))
