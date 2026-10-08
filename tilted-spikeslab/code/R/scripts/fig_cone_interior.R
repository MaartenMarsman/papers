# R/scripts/fig_cone_interior.R
#
# SUPERSEDED (pre-pivot). The manuscript's figure fig:cone-interior is produced
# by fig_cone_interior_expanded.R from cone_interior_expanded.rds. This script
# draws the older three-panel version (q <= 20, sigma = 1, both variants, with
# an eigenvalue-floor annotation) and now writes to a separate file.
#
# Three-panel figure for the cone-interior diagnostic
# (R/results/cone_interior_diagnostics.rds).
#
# Panel A: unconditional Pr[K positive-definite] vs q, log-y, both variants.
# Panel B: median + IQR of lambda_min(K) on PD-conditioned draws vs q,
#          plus the Pr[lambda_min < 0.05 | PD] line.
# Panel C: effective diagonal mean (PD-conditioned) vs q, with the nominal
#          mean alpha/beta = 1 as a reference line.
#
# Output: figures/cone_interior_legacy.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
IN_RDS       <- file.path(PROJECT_ROOT, "R", "results", "cone_interior_diagnostics.rds")
OUT_PDF      <- file.path(PROJECT_ROOT, "figures", "cone_interior_legacy.pdf")

obj <- readRDS(IN_RDS)

# Build a long-form data frame from the per-cell results
rows <- lapply(obj$results, function(r) {
  data.frame(
    variant         = r$variant,
    q               = r$q,
    pd_rate         = r$pd_rate,
    pr_lt05_given_pd = r$pr_lmin_lt05_given_pd,
    lmin_pd_median  = r$lmin_pd_median,
    lmin_pd_q05     = r$lmin_pd_q05,
    lmin_pd_q95     = r$lmin_pd_q95,
    diag_mean_pd    = r$diag_mean_pd,
    stringsAsFactors = FALSE
  )
})
df <- do.call(rbind, rows)
df$variant <- factor(df$variant,
                     levels = c("discrete", "continuous"),
                     labels = c("Discrete (bgms-style)", "Continuous (Wang-style)"))

theme_paper <- theme_minimal(base_size = 10) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background   = element_blank(),
    plot.title         = element_text(size = 10, face = "plain"),
    plot.tag           = element_text(size = 11, face = "bold")
  )

palette_variants <- c(
  "Discrete (bgms-style)"   = "#1f78b4",
  "Continuous (Wang-style)" = "#e31a1c"
)

# Panel A: Pr[PD] vs q (log scale)
pA <- ggplot(df, aes(x = q, y = pd_rate, colour = variant, shape = variant)) +
  geom_line() +
  geom_point(size = 2) +
  scale_y_log10(labels = function(x) sprintf("%g%%", x * 100),
                limits = c(1e-7, 1),
                breaks = c(1e-6, 1e-4, 1e-2, 1)) +
  scale_x_continuous(breaks = c(5, 10, 15, 20)) +
  scale_colour_manual(values = palette_variants) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "q",
       y = expression(Pr * "[K positive-definite]"),
       title = "A. Unconditional PD rate") +
  theme_paper

# Panel B: median + IQR of lambda_min on PD draws + boundary concentration
pB <- ggplot(df, aes(x = q, colour = variant, fill = variant, shape = variant)) +
  geom_ribbon(aes(ymin = lmin_pd_q05, ymax = lmin_pd_q95),
              alpha = 0.15, colour = NA) +
  geom_line(aes(y = lmin_pd_median)) +
  geom_point(aes(y = lmin_pd_median), size = 2) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 5, y = 0.05, label = expression(epsilon == 0.05),
           hjust = 0, vjust = -0.5, size = 2.8, colour = "grey40") +
  scale_x_continuous(breaks = c(5, 10, 15, 20)) +
  scale_colour_manual(values = palette_variants) +
  scale_fill_manual(values = palette_variants) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "q",
       y = expression(lambda[min] * "(K) given PD"),
       title = "B. Boundary concentration") +
  theme_paper

# Panel C: diagonal mean (PD-conditioned) vs nominal
pC <- ggplot(df, aes(x = q, y = diag_mean_pd, colour = variant, shape = variant)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 5, y = 1, label = "nominal mean",
           hjust = 0, vjust = 1.4, size = 2.8, colour = "grey40") +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(5, 10, 15, 20)) +
  scale_colour_manual(values = palette_variants) +
  scale_shape_manual(values = c(16, 17)) +
  ylim(0.9, 1.4) +
  labs(x = "q",
       y = expression(E * "[" * Theta[ll] * " given PD]"),
       title = "C. Diagonal-mean distortion") +
  theme_paper

# Extract a single shared legend from any panel
extract_legend <- function(p) {
  g <- ggplotGrob(p + theme(legend.position = "bottom"))
  leg <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  leg
}
shared_legend <- extract_legend(pA)

strip_legend <- function(p) p + theme(legend.position = "none")

fig <- gridExtra::grid.arrange(
  strip_legend(pA),
  strip_legend(pB),
  strip_legend(pC),
  shared_legend,
  ncol = 3,
  nrow = 2,
  layout_matrix = rbind(c(1, 2, 3), c(4, 4, 4)),
  heights = c(10, 1)
)

ggplot2::ggsave(OUT_PDF, plot = fig, width = 9, height = 3.5, units = "in")
cat(sprintf("Wrote %s\n", OUT_PDF))
