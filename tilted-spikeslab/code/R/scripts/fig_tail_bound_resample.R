# R/scripts/fig_tail_bound_resample.R
#
# Figure 4: empirical tail bound on lambda_min(Theta) under the delta-tilted
# prior, using rejection-sampled (not importance-reweighted) draws conditional
# on positive-definiteness. Pr[lambda_min <= x] is the empirical CDF on
# accepted draws (i.i.d.); no importance weights. Theory (Proposition 3.5)
# predicts log-log slope = delta + 1 asymptotically as x -> 0.

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
IN_RDS       <- file.path(PROJECT_ROOT, "R", "results", "cone_interior_tilt_resample.rds")
OUT_PDF      <- file.path(PROJECT_ROOT, "figures", "tail_bound.pdf")

obj <- readRDS(IN_RDS)
x_grid <- 10^seq(-3, -0.3, length.out = 20)  # 0.001 .. 0.5

min_accepted <- 500L  # need enough draws to estimate the small-x tail

rows <- list()
for (r in obj$results) {
  if (!(r$q %in% c(5L, 10L, 20L))) next            # representative q values
  if (length(r$lmin) < min_accepted) next

  lmin <- r$lmin
  pr_below <- vapply(x_grid, function(x) mean(lmin <= x), numeric(1L))
  keep <- pr_below > 0
  if (!any(keep)) next

  rows[[length(rows) + 1L]] <- data.frame(
    variant  = r$variant,
    q        = r$q,
    delta    = r$delta,
    x        = x_grid[keep],
    pr_below = pr_below[keep],
    stringsAsFactors = FALSE
  )
}
df <- do.call(rbind, rows)
df$variant_label <- factor(df$variant,
                            levels = c("discrete", "continuous"),
                            labels = c("Discrete", "Continuous"))
df$delta_key <- sprintf("delta=%g", df$delta)

palette_delta <- c("delta=0"   = "#a6cee3",
                   "delta=0.5" = "#1f78b4",
                   "delta=1"   = "#33a02c",
                   "delta=2"   = "#e31a1c")
shapes_delta  <- c("delta=0"   = 16,
                   "delta=0.5" = 17,
                   "delta=1"   = 15,
                   "delta=2"   = 18)

# Reference lines y = C * x^(delta + 1) anchored at the smallest x (so the
# reference is an upper envelope on the empirical tail across the panel).
ref_df <- do.call(rbind, by(df, list(df$variant, df$q, df$delta), function(sub) {
  # Anchor at the largest empirical value (which dominates the upper-bound
  # constant on this x grid).
  C_panel <- max(sub$pr_below / sub$x^(sub$delta[1] + 1))
  x_vals <- seq(min(sub$x), max(sub$x), length.out = 50)
  data.frame(
    variant = sub$variant[1],
    q       = sub$q[1],
    delta   = sub$delta[1],
    delta_key = sub$delta_key[1],
    x = x_vals,
    pr_ref = C_panel * x_vals^(sub$delta[1] + 1),
    stringsAsFactors = FALSE
  )
}))
ref_df$variant_label <- factor(ref_df$variant,
                                levels = c("discrete", "continuous"),
                                labels = c("Discrete", "Continuous"))

theme_paper <- theme_minimal(base_size = 9) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = 8),
    panel.grid.minor   = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 9, face = "plain"),
    plot.title         = element_text(size = 9, face = "plain")
  )

fig <- ggplot(df,
              aes(x = x, y = pr_below, colour = delta_key, shape = delta_key)) +
  geom_line(data = ref_df,
            aes(x = x, y = pr_ref, colour = delta_key, group = delta_key),
            linetype = "dashed", alpha = 0.6, inherit.aes = FALSE) +
  geom_point(size = 2) +
  facet_grid(variant_label ~ q,
             labeller = labeller(.cols = function(x) sprintf("q = %s", x))) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = palette_delta, name = expression(delta),
                      labels = c("0", "0.5", "1", "2")) +
  scale_shape_manual(values = shapes_delta, name = expression(delta),
                     labels = c("0", "0.5", "1", "2")) +
  labs(x = expression(x),
       y = expression(Pr * "[" * lambda[min] * "(K) " <= x * "  given PD, delta-tilted]"),
       title = expression("Empirical tail bound vs the " * x^(delta+1) * " upper envelope (dashed)")) +
  theme_paper

ggplot2::ggsave(OUT_PDF, plot = fig, width = 9, height = 5.5, units = "in")
cat(sprintf("Wrote %s\n", OUT_PDF))
