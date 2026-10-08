# R/scripts/fig_partial_correlation.R
#
# The induced prior on the partial correlation rho_{ij} | gamma_{ij} = 1, i.e.
# the "prior you get" on the interpretable conditional-dependence parameter.
# Reads R/results/partial_correlation.rds (from partial_correlation_sweep.R).
#
# Panel A: density of rho on (-1, 1), coloured by q, at the untilted prior
#          (delta = 0) vs the dimension-adaptive tilt (delta ~ 0.5 log q).
#          The point: at delta = 0 the induced rho prior narrows (drifts) with
#          q; the adaptive tilt holds it roughly fixed.
# Panel B: SD(rho) vs q, delta = 0 vs adaptive.
#
# Style:  R/scripts/fig_style.R (shared theme, palette, export contract)
# Data:   R/results/partial_correlation.rds (cached; do not recompute)
# Output: figures/partial_correlation_prior.pdf
#         tables/partial_correlation_prior.tex

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
source(file.path(PROJECT_ROOT, "R", "scripts", "fig_style.R"))

IN_RDS  <- file.path(PROJECT_ROOT, "R", "results", "partial_correlation.rds")
OUT_PDF <- file.path(PROJECT_ROOT, "figures", "partial_correlation_prior.pdf")
OUT_TEX <- file.path(PROJECT_ROOT, "tables", "partial_correlation_prior.tex")

obj <- readRDS(IN_RDS)
results <- obj$results
tab <- obj$summary

# For each q, tag delta = 0 as "untilted" and the delta closest to 0.5 log q
# as "adaptive".
qs <- sort(unique(sapply(results, function(r) r$q)))
pick <- function(q, which) {
  cand <- Filter(function(r) r$q == q && length(r$rho) > 30L, results)
  if (length(cand) == 0) return(NULL)
  ds <- sapply(cand, function(r) r$delta)
  if (which == "untilted") {
    idx <- which.min(abs(ds - 0))
  } else {
    idx <- which.min(abs(ds - 0.5 * log(q)))
  }
  cand[[idx]]
}

# Tilt labels shared by facet strips (panel A) and the colour legend (panel B).
TILT_LEVELS <- c("untilted", "adaptive")
TILT_PARSED <- c(
  untilted = "untilted~(delta == 0)",
  adaptive = "adaptive~(delta %~~% 0.5~log~italic(q))"
)
labels_tilt <- c(
  untilted = expression(untilted~(delta == 0)),
  adaptive = expression(adaptive~(delta %~~% 0.5~log~italic(q)))
)

# --- Panel A: rho density on (-1, 1) ---------------------------------------
kde_row <- function(r, label) {
  d <- density(r$rho, from = -1, to = 1, n = 512, bw = "SJ")
  data.frame(q = r$q, delta = r$delta, tilt = label, x = d$x, density = d$y)
}
densA <- list()
for (q in qs) {
  ru <- pick(q, "untilted"); ra <- pick(q, "adaptive")
  if (!is.null(ru)) densA[[length(densA) + 1L]] <- kde_row(ru, "untilted")
  if (!is.null(ra)) densA[[length(densA) + 1L]] <- kde_row(ra, "adaptive")
}
densA <- do.call(rbind, densA)
densA$q_factor <- factor(densA$q, levels = qs)
densA$tilt_lab <- factor(TILT_PARSED[densA$tilt], levels = TILT_PARSED)

pA <- ggplot(densA, aes(x = x, y = density, colour = q_factor)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60",
             linewidth = LW_GUIDE) +
  geom_line(linewidth = LW_MAIN) +
  facet_wrap(~ tilt_lab, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = pal_q(qs), name = expression(italic(q))) +
  guides(colour = guide_legend(order = 1, nrow = 1)) +
  labs(x = expression(rho[ij] ~ "given" ~ gamma[ij] == 1),
       y = "density") +
  coord_cartesian(xlim = c(-1, 1)) +
  theme_fig()

# --- Panel B: SD(rho) vs q -------------------------------------------------
sumB <- list()
for (q in qs) {
  for (w in TILT_LEVELS) {
    r <- pick(q, w); if (is.null(r)) next
    sumB[[length(sumB) + 1L]] <- data.frame(
      q = q, tilt = w,
      rho_sd = sd(r$rho),
      p_gt = mean(abs(r$rho) > 0.3))
  }
}
sumB <- do.call(rbind, sumB)
sumB$tilt <- factor(sumB$tilt, levels = TILT_LEVELS)

pB <- ggplot(sumB, aes(x = q, y = rho_sd, colour = tilt, group = tilt)) +
  geom_line(linewidth = LW_MAIN) +
  geom_point(size = PT_MAIN) +
  scale_x_continuous(breaks = qs) +
  scale_colour_manual(values = PAL_TILT, name = NULL, labels = labels_tilt) +
  guides(colour = guide_legend(order = 2, nrow = 1)) +
  labs(x = LAB_Q,
       y = expression("SD(" * rho[ij] ~ "given" ~ gamma[ij] == 1 * ")")) +
  theme_fig()

fig <- (pA | pB) +
  plot_layout(guides = "collect", widths = c(1.7, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom", legend.box = "vertical")

save_fig(OUT_PDF, fig, height = fig_height(n_rows = 1L, n_legend_rows = 2L))

# --- LaTeX summary table ---------------------------------------------------
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
st <- tab[order(tab$q, tab$delta), ]
lines <- c(
  "\\begin{tabular}{rrrrrr}",
  "\\hline",
  "$q$ & $\\delta$ & SD($\\rho$) & mean($\\rho$) & $\\Pr(|\\rho|>0.3)$ & $\\Pr(|\\rho|>0.5)$ \\\\",
  "\\hline",
  apply(st, 1, function(row) sprintf("%d & %s & %s & %s & %s & %s \\\\",
        as.integer(row["q"]), fmt(row["delta"], 2),
        fmt(row["rho_sd"]), fmt(row["rho_mean"]),
        fmt(row["p_abs_gt_0.3"]), fmt(row["p_abs_gt_0.5"]))),
  "\\hline",
  "\\end{tabular}"
)
writeLines(lines, OUT_TEX)
cat(sprintf("Wrote %s\n", OUT_TEX))

cat("\nInduced partial-correlation prior summary:\n")
print(st[, c("q", "delta", "n", "rho_sd", "rho_mean", "p_abs_gt_0.3")],
      row.names = FALSE)
