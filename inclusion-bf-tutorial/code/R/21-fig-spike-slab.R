# 21-fig-spike-slab.R
# Concept figure for Section 2: what model-averaged estimation of one edge
# weight looks like. Every edge starts from the same two-part prior (spike:
# probability .5 of exactly zero; slab: standard-normal weight if present).
# The posterior is the same kind of two-part distribution, and the data
# decide how the mass is divided. Three edges of the worked example span the
# three verdicts.
#
# Style follows the JASP prior/posterior plots (R Graph Compendium 11.2):
# filled smooth slab curve, lollipop spike,
# no gridlines, no y-axis numbers. Heights remain commensurable within each
# panel: the spike stem is the probability of exactly zero (printed above
# it), and the slab curve is the probability per weight bin of width .01,
# so spike and slab show how the posterior mass is actually divided.
#
# Reads tables/fit-wenchuan.rds (from R/01-worked-example.R).
# Output: figures/fig-spike-slab.pdf; prints the model-averaged 95% credible
# interval of the focal edge for the manuscript.

set.seed(123)

fit <- readRDS("tables/fit-wenchuan.rds")
samples <- bgms::extract_pairwise_interactions(fit)
bf <- bgms::extract_inclusion_bf(fit)

edges <- list(
  c("intrusion", "concen"),
  c("dreams", "upset"),
  c("flash", "hyper")
)

binwidth <- 0.01
xmax <- 0.35
grid_x <- seq(-xmax, xmax, length.out = 512)

fmt_p <- function(p) {
  if (p < 0.01) "<.01" else sub("^0", "", sprintf("%.2f", p))
}

# Prior panel: spike .5 at zero, slab = half a standard-normal. The slab is
# drawn on the density scale over its own wider axis so the bell is
# visible; the posterior panels zoom to the +-.35 weight range
# the data concentrate on. The y-axis carries no numbers, so the per-panel
# unit is a display choice.
prior_title <- "prior\nidentical for every edge"
prior_x <- seq(-3.2, 3.2, length.out = 512)
curves <- list(data.frame(
  panel = prior_title, x = prior_x,
  y = 0.5 * stats::dnorm(prior_x)
))
spikes <- data.frame(panel = prior_title, p = 0.5, label = fmt_p(0.5))

for (e in edges) {
  draws <- samples[, paste(e[1], e[2], sep = "-")]
  stopifnot(all(abs(draws) < xmax))
  p_spike <- mean(draws == 0)
  verdict <- if (bf[e[1], e[2]] > 10) "presence" else
    if (bf[e[1], e[2]] < 1 / 10) "absence" else "undecided"
  bf_txt <- ifelse(bf[e[1], e[2]] >= 1, sprintf("%.1f", bf[e[1], e[2]]),
                   sprintf("%.2f", bf[e[1], e[2]]))
  title <- sprintf("%s-%s\n%s, BF = %s", e[1], e[2], verdict, bf_txt)
  dens <- stats::density(draws[draws != 0], from = -xmax, to = xmax, n = 512)
  curves[[length(curves) + 1]] <- data.frame(
    panel = title, x = dens$x, y = dens$y * (1 - p_spike) * binwidth
  )
  spikes <- rbind(spikes,
                  data.frame(panel = title, p = p_spike, label = fmt_p(p_spike)))
  cat(sprintf("%s: P(excluded) = %.3f, conditional weight range [%.2f, %.2f]\n",
              gsub("\n", " | ", title), p_spike,
              min(draws[draws != 0]), max(draws[draws != 0])))
}

curve_data <- do.call(rbind, curves)
curve_data$panel <- factor(curve_data$panel, levels = unique(curve_data$panel))
spikes$panel <- factor(spikes$panel, levels = levels(curve_data$panel))

slab_fill <- "#C9DCEE"
slab_line <- "#5B93C4"
spike_col <- "#1F3A63"

fig <- ggplot2::ggplot(curve_data, ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_area(fill = slab_fill, color = slab_line, linewidth = 0.8) +
  ggplot2::geom_segment(
    data = spikes,
    ggplot2::aes(x = 0, xend = 0, y = 0, yend = p),
    color = spike_col, linewidth = 1.3, inherit.aes = FALSE
  ) +
  ggplot2::geom_point(
    data = spikes,
    ggplot2::aes(x = 0, y = p),
    color = spike_col, size = 2.4, inherit.aes = FALSE
  ) +
  ggplot2::geom_text(
    data = spikes,
    ggplot2::aes(x = 0, y = p, label = label),
    vjust = -0.9, size = 3.4, color = spike_col, inherit.aes = FALSE
  ) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
  ggplot2::facet_wrap(~panel, nrow = 1, scales = "free") +
  ggplot2::labs(x = "Edge weight", y = "Probability") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.line.y = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 10.5, lineheight = 1.1),
    panel.spacing = ggplot2::unit(1.1, "lines")
  )

ggplot2::ggsave(
  filename = "figures/fig-spike-slab.pdf",
  plot = fig,
  width = 9,
  height = 2.7
)

# Model-averaged 95% credible interval for the focal edge (manuscript text).
focal <- samples[, "intrusion-concen"]
cat("focal edge: model-averaged mean", round(mean(focal), 3),
    "| 95% CrI", paste(round(stats::quantile(focal, c(0.025, 0.975)), 3),
                       collapse = " to "), "\n")
