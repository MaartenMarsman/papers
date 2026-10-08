# R/scripts/fig_delta_scaling.R
#
# Calibrating the tilt exponent, in the standardized frame (exponential
# diagonal, eta = 1). Panel A: median lambda_min(K) vs q under the rules
# delta = 0, 1, 2, log q (trend only; shared delta palette with
# figures/specified_vs_realized.pdf). Panel B: the tilt required at each q for
# a fixed target median (targets 0.1, 0.5, 1.0), by numerically inverting the
# median-vs-delta relation (full delta grid, including the (q-2)/2 anchors),
# with a single dashed reference curve delta = 0.5 log q labeled in-plot.
#
# Style:  R/scripts/fig_style.R (shared theme, palette, export contract)
# Data:   R/results/delta_scaling_gwishart.rds (cached; do not recompute)
# Output: figures/delta_scaling.pdf

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
source(file.path(PROJECT_ROOT, "R", "scripts", "fig_style.R"))

IN_RDS  <- file.path(PROJECT_ROOT, "R", "results", "delta_scaling_gwishart.rds")
OUT_PDF <- file.path(PROJECT_ROOT, "figures", "delta_scaling.pdf")

obj <- readRDS(IN_RDS)
tab <- obj$summary
tab <- tab[tab$n_accepted >= 25, ]

# Rules shown in panel A; the (q-2)/2 anchors stay in `tab` for the panel-B
# inversion but are not drawn. Colours come from PAL_RULE, which matches the
# delta encoding of figures/specified_vs_realized.pdf (0 grey, 1 blue,
# 2 green, log rule vermillion).
RULES_SHOWN <- c("const0", "const1", "const2", "logq")
tab_shown <- tab[tab$rule %in% RULES_SHOWN, ]
tab_shown$rule <- factor(tab_shown$rule, levels = RULES_SHOWN)

# --- Panel A: median lambda_min vs q (trend only, linear axis) -------------
pA <- ggplot(tab_shown, aes(x = q, y = median, colour = rule, group = rule)) +
  geom_line(linewidth = LW_MAIN) +
  geom_point(size = PT_MAIN) +
  scale_x_continuous(breaks = c(5, 10, 20, 30, 50)) +
  scale_colour_manual(values = PAL_RULE[RULES_SHOWN], name = NULL,
                      labels = LABELS_RULE[RULES_SHOWN]) +
  guides(colour = guide_legend(order = 1, nrow = 1)) +
  labs(x = LAB_Q,
       y = expression("median"~lambda[min](K))) +
  theme_fig()

# --- Panel B: required delta(q) for target medians -------------------------
# Interpolate median(lambda_min) vs delta within each q, using the full rule
# grid (including the (q-2)/2 anchors) as anchor points, and invert for
# targets in TARGET_LIST.
TARGET_LIST <- c(0.1, 0.5, 1.0)

req_rows <- list()
q_vals <- sort(unique(tab$q))
for (q in q_vals) {
  sub <- tab[tab$q == q, ]
  if (nrow(sub) < 3) next
  ord <- order(sub$delta)
  d_v <- sub$delta[ord]; m_v <- sub$median[ord]
  for (target in TARGET_LIST) {
    if (m_v[1] >= target) {
      d_req <- d_v[1]   # already above target with smallest delta
    } else if (max(m_v) < target) {
      d_req <- NA_real_
    } else {
      f <- splinefun(d_v, m_v, method = "monoH.FC")
      d_req <- uniroot(function(x) f(x) - target, c(d_v[1], d_v[length(d_v)]))$root
    }
    req_rows[[length(req_rows) + 1L]] <- data.frame(
      q = q, target = target, delta_req = d_req,
      stringsAsFactors = FALSE)
  }
}
req_df <- do.call(rbind, req_rows)
req_df$target_label <- factor(sprintf("%.1f", req_df$target),
                              levels = sprintf("%.1f", TARGET_LIST))

# Sequential ramp for the targets (light = small target, dark = large),
# consistent with the pal_q() mako ramp used elsewhere.
pal_target <- stats::setNames(
  viridisLite::viridis(length(TARGET_LIST), option = "mako",
                       begin = 0.75, end = 0.30),
  sprintf("%.1f", TARGET_LIST))

# Single closed-form reference: delta = 0.5 log q, labeled in-plot.
qq <- seq(5, 50, by = 1L)
ref_df <- data.frame(q = qq, delta = 0.5 * log(qq))

pB <- ggplot(req_df[!is.na(req_df$delta_req), ],
             aes(x = q, y = delta_req, colour = target_label,
                 group = target_label)) +
  geom_line(data = ref_df,
            aes(x = q, y = delta),
            inherit.aes = FALSE, colour = "black", linetype = "22",
            linewidth = 0.35) +
  annotate("text", x = 50, y = 0.5 * log(50), label = "0.5~log~italic(q)",
           parse = TRUE, hjust = 1, vjust = -0.8, size = FS_NOTE / .pt,
           colour = "grey20") +
  geom_line(linewidth = LW_MAIN) +
  geom_point(size = PT_MAIN) +
  scale_x_continuous(breaks = c(5, 10, 20, 30, 50)) +
  scale_colour_manual(values = pal_target, name = "target median") +
  guides(colour = guide_legend(order = 2, nrow = 1)) +
  coord_cartesian(ylim = c(0, 3.2)) +
  labs(x = LAB_Q, y = expression("required tilt"~delta)) +
  theme_fig()

fig <- (pA | pB) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom", legend.box = "horizontal")

save_fig(OUT_PDF, fig, height = fig_height(n_rows = 1L, n_legend_rows = 1L))

cat("\nRequired delta(q) for each target:\n")
print(req_df, row.names = FALSE)
