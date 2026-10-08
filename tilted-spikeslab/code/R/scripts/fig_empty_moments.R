# R/scripts/fig_empty_moments.R
#
# Figure 4: empty-graph moment verification. For each q in {5, 20, 50} and
# delta in {0, 0.5, 1, 2}, compare empirical and closed-form values of three
# representative moments: E[theta_ll], E[||Theta||_F^2], E[log|Theta|].
# Data: R/results/empty_graph_moments.rds.

suppressPackageStartupMessages({
  library(ggplot2)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
IN_RDS       <- file.path(PROJECT_ROOT, "R", "results", "empty_graph_moments.rds")
OUT_PDF      <- file.path(PROJECT_ROOT, "figures", "empty_moments.pdf")

obj <- readRDS(IN_RDS)
tab <- obj$summary
metrics_keep <- c("E_theta_ll", "E_frob_sq", "E_log_det")
tab <- tab[tab$metric %in% metrics_keep, ]
tab$metric <- factor(tab$metric,
                      levels = metrics_keep,
                      labels = c(expression(E * "[" * theta[ll] * "]"),
                                 expression(E * "[||" * Theta * "||"[F]^{2} * "]"),
                                 expression(E * "[log|" * Theta * "|]")))
tab$q_label <- factor(sprintf("q=%d", tab$q),
                      levels = c("q=5", "q=20", "q=50"))

theme_paper <- theme_minimal(base_size = 9) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = 8),
    panel.grid.minor   = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 9, face = "plain"),
    plot.title         = element_text(size = 9, face = "plain")
  )

palette_q <- c("q=5" = "#a6cee3", "q=20" = "#1f78b4", "q=50" = "#33a02c")
shapes_q  <- c("q=5" = 16, "q=20" = 17, "q=50" = 15)

fig <- ggplot(tab,
              aes(x = delta, colour = q_label, shape = q_label, group = q_label)) +
  geom_line(aes(y = closed_form)) +
  geom_point(aes(y = empirical), size = 2) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3, labeller = label_parsed) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 2)) +
  scale_colour_manual(values = palette_q, name = "q") +
  scale_shape_manual(values = shapes_q, name = "q") +
  labs(x = expression(tilt~exponent~delta),
       y = "moment value",
       title = "Empty-graph moments: closed-form (lines) vs empirical (points)") +
  theme_paper

ggplot2::ggsave(OUT_PDF, plot = fig, width = 9, height = 3.2, units = "in")
cat(sprintf("Wrote %s\n", OUT_PDF))
