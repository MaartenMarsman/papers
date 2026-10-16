# 24-fig-calibration.R
# PAV-adjusted calibration plots for the worked example (Robustness section),
# one panel per symptom, following the visual-predictive-check
# recommendations of Sailynoja et al. (arXiv 2503.01509) and the CORP
# reliability diagram (isotonic recalibration via pool-adjacent-violators).
#
# The model's person-level predictions are the conditional distributions
# P(x_ij | rest of person i's responses); the conditional distribution of
# one variable given the rest is the model's own regression unit, so
# calibration is assessed per variable (pooling across variables could hide
# compensating miscalibration). For every person and category threshold k
# the model issues P(x_ij <= k | rest) and the data record the outcome; PAV
# (isotonic regression) estimates observed frequency as a monotone function
# of predicted probability. A calibrated model tracks the diagonal. The 95%
# consistency band resamples each person's CATEGORY from its predicted
# distribution (the four threshold events are nested, so events must not be
# resampled independently) and shows the wander of a calibrated-by-
# construction model.
#
# Reads tables/fit-wenchuan.rds (from R/01-worked-example.R).
# Outputs: figures/fig-calibration.pdf, tables/calibration-summary.csv

library(bgms)

set.seed(123)

fit <- readRDS("tables/fit-wenchuan.rds")
data("Wenchuan", package = "bgms")
wenchuan <- Wenchuan[stats::complete.cases(Wenchuan), ]

# method = "posterior-mean" plugs the posterior mean of the parameters into
# the conditional distributions rather than averaging predictions over draws.
probs <- predict(fit, newdata = wenchuan, type = "probabilities",
                 method = "posterior-mean")

pav_curve <- function(p, y, grid) {
  ir <- stats::isoreg(p, y)
  stats::approx(ir$x[ir$ord], ir$yf, xout = grid, method = "linear",
                rule = 2, ties = "ordered")$y
}

grid <- seq(0, 1, length.out = 101)
n_rep <- 200

curves <- list()
summary_rows <- list()
for (v in names(probs)) {
  pm <- probs[[v]]                                   # n x 5 category probs
  n_cat <- ncol(pm)
  n_thresh <- n_cat - 1
  cum <- t(apply(pm, 1, cumsum))[, -n_cat]           # n x 4 cumulative
  # assumes categories coded 1 to 5 in the data and predict() columns 0 to 4
  obs_cat <- wenchuan[, v] - 1
  p <- as.vector(cum)
  y <- as.integer(outer(obs_cat, seq_len(n_thresh) - 1, `<=`))

  obs_curve <- pav_curve(p, y, grid)
  band <- replicate(n_rep, {
    cat_star <- apply(pm, 1, function(pr) sample.int(n_cat, 1, prob = pr) - 1)
    pav_curve(p, as.integer(outer(cat_star, seq_len(n_thresh) - 1, `<=`)), grid)
  })
  lo <- apply(band, 1, stats::quantile, 0.025)
  hi <- apply(band, 1, stats::quantile, 0.975)

  curves[[v]] <- data.frame(variable = v, grid = grid, pav = obs_curve,
                            lo = lo, hi = hi)
  summary_rows[[v]] <- data.frame(
    variable = v,
    mean_dev = mean(abs(obs_curve - grid)),
    max_dev = max(abs(obs_curve - grid)),
    share_outside_band = mean(obs_curve < lo | obs_curve > hi)
  )
}

curve_data <- do.call(rbind, curves)
summary_df <- do.call(rbind, summary_rows)
summary_df <- summary_df[order(-summary_df$max_dev), ]
utils::write.csv(summary_df, "tables/calibration-summary.csv", row.names = FALSE)

cat("per-variable calibration (worst first):\n")
print(transform(summary_df, mean_dev = round(mean_dev, 3),
                max_dev = round(max_dev, 3),
                share_outside_band = round(share_outside_band, 2)),
      row.names = FALSE)

# --- Figure: 17 small multiples, styled like the other check figures -----------

slab_fill <- "#C9DCEE"
dark_blue <- "#1F3A63"

fig <- ggplot2::ggplot(curve_data, ggplot2::aes(x = grid)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                       fill = slab_fill, alpha = 0.9) +
  ggplot2::geom_abline(intercept = 0, slope = 1, color = "grey60",
                       linewidth = 0.3) +
  ggplot2::geom_line(ggplot2::aes(y = pav), color = dark_blue,
                     linewidth = 0.55) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  ggplot2::scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  ggplot2::facet_wrap(~variable, ncol = 6) +
  ggplot2::labs(
    x = "Predicted probability",
    y = "Observed frequency (PAV-adjusted)"
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 8.5),
    axis.text = ggplot2::element_text(size = 7),
    panel.spacing = ggplot2::unit(0.55, "lines")
  )

ggplot2::ggsave("figures/fig-calibration.pdf", fig, width = 9, height = 5.4)
