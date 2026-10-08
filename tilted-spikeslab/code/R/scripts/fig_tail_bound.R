# R/scripts/fig_tail_bound.R
#
# Figure 3: empirical tail bound on lambda_min(Theta) under the tilted prior.
# For each (variant, q) at sigma = 1, plot log Pr[lambda_min < x] vs log x for
# delta in {0, 0.5, 1, 2}, computed by importance-reweighting the PD draws by
# |K|^delta. Theory (Proposition 3.5) predicts the log-log slope = delta + 1.
# Reference lines x^(delta+1) are overlaid for each delta.

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
IN_RDS       <- file.path(PROJECT_ROOT, "R", "results", "cone_interior_expanded.rds")
OUT_PDF      <- file.path(PROJECT_ROOT, "figures", "tail_bound.pdf")

obj <- readRDS(IN_RDS)

x_grid <- 10^seq(-3, -0.3, length.out = 20)  # 0.001 .. 0.5

rows <- list()
for (r in obj$results) {
  if (abs(r$sigma - 1) > 1e-12) next                # only sigma = 1
  if (!(r$q %in% c(5L, 10L, 20L))) next             # representative q values
  if (r$n_pd < 200L) next                           # need enough PD samples

  lmin    <- r$lmin_pd_samples
  log_det <- r$log_det_pd_samples
  if (length(lmin) < 200L) next

  for (d in c(0, 0.5, 1, 2)) {
    log_w <- d * log_det
    log_w <- log_w - max(log_w)
    w <- exp(log_w)
    w_norm <- w / sum(w)

    pr_below <- vapply(x_grid, function(x) sum(w_norm * (lmin <= x)), numeric(1L))
    keep <- pr_below > 0
    if (!any(keep)) next

    rows[[length(rows) + 1L]] <- data.frame(
      variant  = r$variant,
      q        = r$q,
      delta    = d,
      x        = x_grid[keep],
      pr_below = pr_below[keep],
      stringsAsFactors = FALSE
    )
  }
}
df <- do.call(rbind, rows)
df$variant_label <- factor(df$variant,
                            levels = c("discrete", "continuous"),
                            labels = c("Discrete", "Continuous"))
df$delta_label <- factor(sprintf("delta=%g", df$delta),
                          levels = c("delta=0", "delta=0.5", "delta=1", "delta=2"),
                          labels = c(expression(delta == 0),
                                     expression(delta == 0.5),
                                     expression(delta == 1),
                                     expression(delta == 2)))

theme_paper <- theme_minimal(base_size = 9) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = 8),
    panel.grid.minor   = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 9, face = "plain"),
    plot.title         = element_text(size = 9, face = "plain")
  )

palette_delta <- c("delta=0"   = "#a6cee3",
                   "delta=0.5" = "#1f78b4",
                   "delta=1"   = "#33a02c",
                   "delta=2"   = "#e31a1c")
shapes_delta  <- c("delta=0"   = 16,
                   "delta=0.5" = 17,
                   "delta=1"   = 15,
                   "delta=2"   = 18)
df$delta_key <- sprintf("delta=%g", df$delta)

# Reference lines y = C * x^(delta + 1) anchored at the largest x to make the
# slope visible. One reference per (variant, q, delta) panel.
ref_df <- do.call(rbind, by(df, list(df$variant, df$q, df$delta), function(sub) {
  x_anchor <- max(sub$x)
  y_anchor <- max(sub$pr_below[sub$x == x_anchor])
  x_vals <- seq(min(sub$x), max(sub$x), length.out = 50)
  data.frame(
    variant = sub$variant[1],
    q       = sub$q[1],
    delta   = sub$delta[1],
    delta_key = sub$delta_key[1],
    x = x_vals,
    pr_ref = y_anchor * (x_vals / x_anchor)^(sub$delta[1] + 1),
    stringsAsFactors = FALSE
  )
}))

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
       y = expression(Pr * "[" * lambda[min] * "(K) " <= x * "  given PD, tilt-weighted]"),
       title = expression("Empirical tail bound vs the " * x^(delta+1) * " reference (dashed)")) +
  theme_paper

ggplot2::ggsave(OUT_PDF, plot = fig, width = 9, height = 5.5, units = "in")
cat(sprintf("Wrote %s\n", OUT_PDF))
