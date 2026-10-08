# R/scripts/fig_realized_conditional_slab.R
#
# Mechanism figure for "the slab you specify vs the slab you get", at the
# CONDITIONAL level (no Monte Carlo over the marginal). Fix every entry of K
# except k_{ij}. Then |K| is a downward parabola in k_{ij}, vanishing at the
# two PD boundaries, so the exact conditional of k_{ij} given the rest is
#
#   p(k | rest, gamma=1)  ∝  exp(-k^2 / 2 sigma^2)  *  (1 - (k - m)^2 / r^2)^delta,
#                            |k - m| < r,   r = sqrt(a c).
#
# i.e. the SPECIFIED Gaussian slab N(0, sigma^2) times a compactly supported
# "PD window" centred at m with half-width r and taper exponent delta. This
# script draws that product against the specified Gaussian for delta in {0,1,2}
# at a representative rest-block geometry (m, r) read from marginal_offdiag.rds.
#
# Style:  R/scripts/fig_style.R (shared theme, palette, export contract)
# Data:   R/results/marginal_offdiag.rds (cached; closed-form curves fresh)
# Output: figures/realized_conditional_slab.pdf

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
source(file.path(PROJECT_ROOT, "R", "scripts", "fig_style.R"))

IN_RDS  <- file.path(PROJECT_ROOT, "R", "results", "marginal_offdiag.rds")
OUT_PDF <- file.path(PROJECT_ROOT, "figures", "realized_conditional_slab.pdf")

SIGMA <- 1.0
Q_REP <- 10L      # representative dimension for the illustrative geometry

# --- representative rest-block geometry (m, r) -----------------------------
# Use median centre m and median half-width r = sqrt(a c) among on-edge draws
# at q = Q_REP (delta = 1 draws if present). Fall back to defaults if the
# results file is unavailable.
m_rep <- 0.0
r_rep <- 1.5
if (file.exists(IN_RDS)) {
  obj <- readRDS(IN_RDS)
  cand <- Filter(function(r) isTRUE(r$q == Q_REP), obj$results)
  # prefer the delta = 1 cell
  d1 <- Filter(function(r) isTRUE(r$delta == 1), cand)
  src <- if (length(d1) > 0) d1[[1]] else if (length(cand) > 0) cand[[1]] else NULL
  if (!is.null(src) && !is.null(src$m_on) && !is.null(src$a_on) && !is.null(src$c_on)) {
    rr <- sqrt(src$a_on * src$c_on)
    m_rep <- median(src$m_on, na.rm = TRUE)
    r_rep <- median(rr[is.finite(rr)], na.rm = TRUE)
  }
}
cat(sprintf("Representative geometry at q = %d: m = %.3f, r = %.3f\n",
            Q_REP, m_rep, r_rep))

# --- densities -------------------------------------------------------------
deltas <- c(0, 1, 2)
kgrid  <- seq(-3.5, 3.5, length.out = 801L)

normalise <- function(x, f) {
  # trapezoidal normalisation to unit integral
  ok <- is.finite(f) & f >= 0
  f[!ok] <- 0
  area <- sum((f[-1] + f[-length(f)]) / 2 * diff(x))
  if (area <= 0) return(f)
  f / area
}

specified <- dnorm(kgrid, 0, SIGMA)   # the slab you specify

rows <- list()
for (d in deltas) {
  u2 <- (kgrid - m_rep)^2 / r_rep^2
  window <- ifelse(u2 < 1, (1 - u2)^d, 0)          # PD window (unnormalised)
  realized <- normalise(kgrid, dnorm(kgrid, 0, SIGMA) * window)
  window_n <- normalise(kgrid, window)
  rows[[length(rows) + 1L]] <- data.frame(
    k = kgrid,
    delta = d,
    specified = specified,
    window = window_n,
    realized = realized
  )
}
df <- do.call(rbind, rows)

long <- rbind(
  data.frame(k = df$k, delta = df$delta, density = df$specified,
             curve = "specified"),
  data.frame(k = df$k, delta = df$delta, density = df$window,
             curve = "window"),
  data.frame(k = df$k, delta = df$delta, density = df$realized,
             curve = "realized")
)
long$curve <- factor(long$curve, levels = c("specified", "window", "realized"))
long$delta_lab <- factor(sprintf("delta == %g", long$delta),
                         levels = sprintf("delta == %g", deltas))

labels_curve <- c(
  specified = expression("specified slab" ~ N(0, sigma^2)),
  window    = "PD window",
  realized  = "realized conditional slab"
)
lt_curve <- c(specified = "dashed", window = "solid", realized = "solid")

p <- ggplot(long, aes(x = k, y = density, colour = curve, linetype = curve)) +
  geom_vline(xintercept = c(m_rep - r_rep, m_rep + r_rep),
             linetype = "dotted", colour = "grey60", linewidth = LW_GUIDE) +
  geom_line(linewidth = LW_MAIN) +
  facet_wrap(~ delta_lab, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = PAL_CURVE, labels = labels_curve, name = NULL) +
  scale_linetype_manual(values = lt_curve, labels = labels_curve, name = NULL) +
  labs(
    x = expression(k[ij] ~ "given the rest," ~ gamma[ij] == 1),
    y = "density") +
  coord_cartesian(xlim = c(-3, 3)) +
  theme_fig()

save_fig(OUT_PDF, p, height = fig_height(n_rows = 1L, n_legend_rows = 1L) + 0.15)
