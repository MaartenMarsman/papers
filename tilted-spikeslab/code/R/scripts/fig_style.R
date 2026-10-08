# R/scripts/fig_style.R
#
# Shared figure style for the manuscript figure family.
# Source this file from every fig_*.R script:
#
#   source(file.path(PROJECT_ROOT, "R", "scripts", "fig_style.R"))
#
# Export contract
# ---------------
# The manuscript is a single-column arxiv article; every figure is imported
# with \includegraphics[width=\linewidth] at a text width of ~6.5 in.
# Figures are therefore exported at their final native width FIG_WIDTH and
# imported at 1:1 scale, so all point sizes below are true printed sizes.
# Heights derive from the benchmark row height FIG_ROW_HEIGHT plus legend rows.
#
# Design registry
# ---------------
# - One theme (theme_fig) for every figure; text sizes are never retuned
#   per plot.
# - Okabe-Ito colorblind-safe palette; black/grey reserved for reference
#   curves and benchmarks.
# - Named palettes, linetypes, and shapes below are the single source of
#   truth for the semantic encodings shared across figures.
# - Multi-panel figures are composed with patchwork; panel tags A, B, C come
#   from patchwork::plot_annotation(tag_levels = "A").
# - One bottom legend per figure (patchwork guides = "collect").
# - Export via save_fig() (cairo_pdf).

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
  library(scales)
})

# --- sizing contract -------------------------------------------------------
FIG_WIDTH      <- 6.5   # in; manuscript text width (import at 1:1)
FIG_ROW_HEIGHT <- 2.15  # in; benchmark height of one panel row (excl. legend)
FIG_LEGEND_ROW <- 0.40  # in; height of one bottom legend row

fig_height <- function(n_rows = 1L, n_legend_rows = 1L) {
  n_rows * FIG_ROW_HEIGHT + n_legend_rows * FIG_LEGEND_ROW
}

# --- typography (points, at final size) ------------------------------------
FS_BASE  <- 10    # axis titles
FS_SMALL <- 8.5   # axis tick labels
FS_LEG   <- 9     # legend title and text
FS_STRIP <- 9.5   # facet strip labels
FS_TAG   <- 11    # panel tags A, B, C
FS_NOTE  <- 7.5   # in-panel annotations; use size = FS_NOTE / .pt in annotate()

theme_fig <- function() {
  theme_minimal(base_size = FS_BASE) +
    theme(
      axis.title         = element_text(size = FS_BASE),
      axis.text          = element_text(size = FS_SMALL, colour = "grey30"),
      legend.position    = "bottom",
      legend.title       = element_text(size = FS_LEG),
      legend.text        = element_text(size = FS_LEG),
      legend.key.height  = unit(10, "pt"),
      legend.key.width   = unit(18, "pt"),
      legend.margin      = margin(0, 0, 0, 0),
      legend.box.spacing = unit(3, "pt"),
      legend.spacing.y   = unit(1, "pt"),
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(linewidth = 0.25, colour = "grey88"),
      strip.background   = element_blank(),
      strip.text         = element_text(size = FS_STRIP),
      plot.tag           = element_text(size = FS_TAG, face = "bold"),
      plot.margin        = margin(3, 5, 1, 3)
    )
}

# --- palette: Okabe-Ito ----------------------------------------------------
OKABE_ITO <- c(
  orange     = "#E69F00",
  skyblue    = "#56B4E9",
  green      = "#009E73",
  yellow     = "#F0E442",
  blue       = "#0072B2",
  vermillion = "#D55E00",
  purple     = "#CC79A7",
  black      = "#000000",
  grey       = "#999999"
)

# Slab scales sigma (cone-interior figure); keyed by the numeric value as text.
PAL_SIGMA <- c(
  "0.3" = OKABE_ITO[["blue"]],
  "1"   = OKABE_ITO[["orange"]],
  "3"   = OKABE_ITO[["vermillion"]]
)

# Ordered dimension palette (dark = larger q), colorblind-safe sequential ramp.
pal_q <- function(q_levels) {
  cols <- viridisLite::viridis(length(q_levels), option = "mako",
                               begin = 0.80, end = 0.15)
  stats::setNames(cols, as.character(q_levels))
}

# Tilt rules delta(q); grey reserved for the untilted reference delta = 0.
PAL_RULE <- c(
  const0 = OKABE_ITO[["grey"]],
  const1 = OKABE_ITO[["blue"]],
  const2 = OKABE_ITO[["green"]],
  logq   = OKABE_ITO[["vermillion"]],
  linear = OKABE_ITO[["purple"]]
)
LABELS_RULE <- c(
  const0 = expression(delta == 0),
  const1 = expression(delta == 1),
  const2 = expression(delta == 2),
  logq   = expression(delta == log~italic(q)),
  linear = expression(delta == (italic(q) - 2) / 2)
)

# Untilted vs dimension-adaptive tilt (partial-correlation figure).
PAL_TILT <- c(
  untilted = OKABE_ITO[["blue"]],
  adaptive = OKABE_ITO[["vermillion"]]
)

# Curve roles in specified-vs-realized figures; black is the reference.
PAL_CURVE <- c(
  specified = OKABE_ITO[["black"]],
  window    = OKABE_ITO[["skyblue"]],
  realized  = OKABE_ITO[["vermillion"]]
)

# --- shared secondary encodings --------------------------------------------
# Discrete vs continuous spike-and-slab variants: solid vs dashed, filled vs
# open points (the current manuscript shows the discrete variant only).
LT_VARIANT    <- c(Discrete = "solid", Continuous = "22")
SHAPE_VARIANT <- c(Discrete = 16, Continuous = 1)

# Default linewidth and point-size tiers.
LW_MAIN  <- 0.55   # primary curves
LW_REF   <- 0.55   # reference curves (specified slab, rule overlays)
LW_GUIDE <- 0.30   # guide lines (hline/vline markers)
PT_MAIN  <- 1.6    # data points on trend lines

# --- shared axis labels ----------------------------------------------------
LAB_Q <- expression("dimension"~italic(q))

# --- export helper ---------------------------------------------------------
# Prefer cairo_pdf; fall back to the standard pdf device when cairo is not
# compiled in (e.g. macOS without XQuartz). The default sans family maps to
# Helvetica, a base PDF font, so the fallback remains portable.
.CAIRO_OK <- local({
  ok <- FALSE
  f <- tempfile(fileext = ".pdf")
  tryCatch({
    suppressWarnings(grDevices::cairo_pdf(f))
    grDevices::dev.off()
    ok <- file.exists(f) && file.size(f) > 0
  }, error = function(e) NULL)
  unlink(f)
  ok
})

save_fig <- function(path, plot, height, width = FIG_WIDTH) {
  dev <- if (.CAIRO_OK) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(path, plot = plot, width = width, height = height,
                  units = "in", device = dev)
  cat(sprintf("Wrote %s (%.2f x %.2f in)\n", path, width, height))
}
