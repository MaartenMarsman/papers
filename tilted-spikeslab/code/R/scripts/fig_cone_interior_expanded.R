# R/scripts/fig_cone_interior_expanded.R
#
# Figure fig:cone-interior of the manuscript: cone-interior collapse across q and sigma, discrete spike-and-slab
# variant only. Two panels: (A) unconditional PD rate (log10 scale, with
# below-MC-resolution markers), (B) median smallest eigenvalue given PD.
# Colour encodes the slab scale sigma. One message per panel; the continuous
# variant and the diagonal-inflation panel were dropped in the 2026-07 redesign.
#
# Style:  R/scripts/fig_style.R (shared theme, palette, export contract)
# Data:   R/results/cone_interior_expanded.rds (cached; do not recompute)
# Output: figures/cone_interior.pdf

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
source(file.path(PROJECT_ROOT, "R", "scripts", "fig_style.R"))

IN_RDS  <- file.path(PROJECT_ROOT, "R", "results", "cone_interior_expanded.rds")
OUT_PDF <- file.path(PROJECT_ROOT, "figures", "cone_interior.pdf")

obj <- readRDS(IN_RDS)
tab <- obj$summary
tab <- tab[tab$variant == "discrete", ]

# Drop cells with n_pd < 30 from panel B (median unreliable)
tab$lmin_pd_median[tab$n_pd < 30L] <- NA_real_

# Keep as NA the cells with 0 PD draws (can't claim a rate below 1/n_draws).
tab$pd_rate_plot <- ifelse(tab$pd_rate > 0, tab$pd_rate, NA_real_)

tab$sigma_label <- factor(tab$sigma, levels = c(0.3, 1, 3))

# Below-resolution markers: for cells with 0 PD draws, plot an open triangle at
# the MC-resolution floor 1/n_draws to indicate the line is truncated.
floor_df <- tab[tab$pd_rate == 0, ]
floor_df$pd_rate_floor <- 1 / floor_df$n_draws

# Shared aesthetic layers for the two panels.
layers_common <- list(
  geom_line(linewidth = LW_MAIN),
  geom_point(size = PT_MAIN),
  scale_x_continuous(breaks = c(5, 10, 20, 30, 50)),
  scale_colour_manual(values = PAL_SIGMA, name = expression(sigma)),
  guides(colour = guide_legend(nrow = 1,
                               override.aes = list(linewidth = LW_MAIN))),
  theme_fig()
)

aes_common <- aes(x = q, colour = sigma_label, group = sigma_label)

pA <- ggplot(tab[!is.na(tab$pd_rate_plot), ],
             modifyList(aes_common, aes(y = pd_rate_plot))) +
  layers_common +
  geom_point(data = floor_df,
             aes(x = q, y = pd_rate_floor, colour = sigma_label,
                 group = sigma_label),
             shape = 25, fill = "white", size = PT_MAIN, stroke = 0.5,
             inherit.aes = FALSE, show.legend = FALSE) +
  annotate("text", x = 50, y = 10^-6.4, label = "below MC resolution",
           hjust = 1, vjust = 1, size = FS_NOTE / .pt, colour = "grey40") +
  scale_y_log10(breaks = 10^c(0, -2, -4, -6),
                labels = scales::label_log()) +
  labs(x = LAB_Q, y = expression(Pr * "[" * K~"positive definite]"))

pB <- ggplot(tab[!is.na(tab$lmin_pd_median), ],
             modifyList(aes_common, aes(y = lmin_pd_median))) +
  geom_hline(yintercept = 0, linewidth = LW_GUIDE, colour = "grey60") +
  layers_common +
  labs(x = LAB_Q,
       y = expression("median"~lambda[min](K)~"given PD"))

fig <- (pA | pB) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_fig(OUT_PDF, fig, height = fig_height(n_rows = 1L, n_legend_rows = 1L))
