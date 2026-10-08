# R/scripts/fig_specified_vs_realized.R
#
# "The prior you specify vs the prior you get", at one representative
# dimension (q = 10), merging the former marginal_specified_vs_realized and
# partial_correlation_prior figures into a single two-panel figure.
#
# Panel A: specified slab N(0, sigma^2) (black dashed) vs the realized
#          marginal density of an included off-diagonal entry k_ij | gamma = 1
#          at q = 10, for delta = 0, 1, 2. Realized densities via the Chen
#          (1994) conditional marginal density estimator (CMDE), reusing the
#          cached (a, c, m) rest-block geometry.
# Panel B: induced prior density on the partial correlation rho_ij | gamma = 1
#          at q = 10: untilted (delta = 0), adaptive (delta = 0.5 log q =
#          1.151 at q = 10), and delta = 2. Zero guide line.
#
# One shared colour legend for delta across both panels (delta = 1 appears in
# panel A only, delta = 0.5 log q in panel B only).
#
# Style:  R/scripts/fig_style.R (shared theme, palette, export contract)
# Data:   R/results/marginal_offdiag.rds      (cached; do not recompute)
#         R/results/partial_correlation.rds   (cached; do not recompute)
# Output: figures/specified_vs_realized.pdf

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
source(file.path(PROJECT_ROOT, "R", "scripts", "fig_style.R"))

IN_MARG <- file.path(PROJECT_ROOT, "R", "results", "marginal_offdiag.rds")
IN_PCOR <- file.path(PROJECT_ROOT, "R", "results", "partial_correlation.rds")
OUT_PDF <- file.path(PROJECT_ROOT, "figures", "specified_vs_realized.pdf")

Q_SHOW <- 10L

# --- shared delta encoding -------------------------------------------------
# One colour scale across both panels; keys ordered by the delta value at
# q = 10 (0 < 1 < 0.5 log 10 = 1.15 < 2). Grey stays reserved for the
# untilted reference delta = 0; the adaptive rule reuses the logq rule colour.
DELTA_LEVELS <- c("d0", "d1", "dlog", "d2")
PAL_DELTA <- c(
  d0   = PAL_RULE[["const0"]],
  d1   = PAL_RULE[["const1"]],
  dlog = PAL_RULE[["logq"]],
  d2   = PAL_RULE[["const2"]]
)
LABELS_DELTA <- c(
  d0   = expression(delta == 0),
  d1   = expression(delta == 1),
  dlog = expression(delta == 0.5~log~italic(q)),
  d2   = expression(delta == 2)
)

scale_delta <- function() {
  scale_colour_manual(values = PAL_DELTA, labels = LABELS_DELTA,
                      limits = DELTA_LEVELS, drop = FALSE, name = NULL)
}

# --- Panel A: specified vs realized marginal slab at q = 10 ----------------
marg <- readRDS(IN_MARG)
SIGMA <- marg$params$sigma

# CMDE machinery (as in fig_marginal_specified_vs_realized.R).
cond_density_unn <- function(y, a, c, m, delta, sigma) {
  r2 <- a * c
  u2 <- (y - m)^2 / r2
  inside <- (u2 < 1) & (r2 > 0)
  ifelse(inside, (1 - u2)^delta * exp(-y^2 / (2 * sigma^2)), 0)
}
sample_normalisers <- function(a_vec, c_vec, m_vec, delta, sigma, n_grid = 80L) {
  n <- length(a_vec); r <- sqrt(a_vec * c_vec); Zs <- numeric(n)
  for (i in seq_len(n)) {
    if (!is.finite(r[i]) || r[i] <= 0) { Zs[i] <- NA_real_; next }
    ys <- seq(m_vec[i] - r[i], m_vec[i] + r[i], length.out = n_grid)
    vals <- cond_density_unn(ys, a_vec[i], c_vec[i], m_vec[i], delta, sigma)
    Zs[i] <- sum((vals[-1] + vals[-n_grid]) / 2) * (ys[2] - ys[1])
  }
  Zs
}
cmde_density <- function(x_grid, a_vec, c_vec, m_vec, delta, sigma, n_norm_grid = 80L) {
  Zs <- sample_normalisers(a_vec, c_vec, m_vec, delta, sigma, n_norm_grid)
  ok <- is.finite(Zs) & Zs > 0
  a_vec <- a_vec[ok]; c_vec <- c_vec[ok]; m_vec <- m_vec[ok]; Zs <- Zs[ok]
  n <- length(a_vec); dens <- numeric(length(x_grid))
  for (j in seq_along(x_grid)) {
    contribs <- cond_density_unn(rep(x_grid[j], n), a_vec, c_vec, m_vec, delta, sigma) / Zs
    dens[j] <- mean(contribs)
  }
  dens
}

x_grid <- seq(-4, 4, length.out = 401L)
dx <- x_grid[2] - x_grid[1]
specified <- dnorm(x_grid, 0, SIGMA)
trapz <- function(y, dx) sum((y[-1] + y[-length(y)]) / 2) * dx

DELTA_KEY_A <- c("0" = "d0", "1" = "d1", "2" = "d2")
densA_rows <- list()
for (r in marg$results) {
  if (r$q != Q_SHOW) next
  if (is.null(r$a_on) || length(r$k12_on_edge) < 50L) next
  key <- unname(DELTA_KEY_A[as.character(r$delta)])
  if (is.na(key)) next
  dens <- cmde_density(x_grid, r$a_on, r$c_on, r$m_on, r$delta, SIGMA)
  area <- trapz(dens, dx); if (area > 0) dens <- dens / area  # guard
  densA_rows[[length(densA_rows) + 1L]] <- data.frame(
    delta_key = key, delta = r$delta, x = x_grid, density = dens)
}
densA <- do.call(rbind, densA_rows)
densA$delta_key <- factor(densA$delta_key, levels = DELTA_LEVELS)
spec_df <- data.frame(x = x_grid, density = specified)

pA <- ggplot(densA, aes(x = x, y = density, colour = delta_key,
                        group = delta_key)) +
  geom_line(data = spec_df,
            aes(x = x, y = density, linetype = "specified"),
            inherit.aes = FALSE, colour = PAL_CURVE[["specified"]],
            linewidth = LW_REF) +
  # show.legend trains colour key glyphs for all four delta levels so the two
  # panel legends are identical and patchwork collects them into one; it is
  # restricted to colour so the dashed linetype key stays intact.
  geom_line(linewidth = LW_MAIN, show.legend = c(colour = TRUE)) +
  scale_delta() +
  scale_linetype_manual(values = c(specified = "dashed"), name = NULL,
                        labels = expression("specified slab" ~ N(0, sigma^2))) +
  guides(linetype = guide_legend(
           order = 1,
           override.aes = list(linetype = "dashed",
                               colour = PAL_CURVE[["specified"]],
                               linewidth = LW_REF)),
         colour = guide_legend(order = 2, nrow = 1)) +
  labs(x = expression(k[ij] ~ "given" ~ gamma[ij] == 1),
       y = "density") +
  coord_cartesian(xlim = c(-3, 3)) +
  theme_fig()

# --- Panel B: induced partial-correlation prior at q = 10 ------------------
pcor <- readRDS(IN_PCOR)
delta_adapt <- 0.5 * log(Q_SHOW)

cand <- Filter(function(r) r$q == Q_SHOW && length(r$rho) > 30L, pcor$results)
ds <- sapply(cand, function(r) r$delta)
pick_nearest <- function(target) cand[[which.min(abs(ds - target))]]

sel <- list(d0 = pick_nearest(0), dlog = pick_nearest(delta_adapt),
            d2 = pick_nearest(2))

# Smooth KDE on (-1, 1): (i) symmetrize the draws (the prior on rho is
# exactly symmetric, so appending -rho is variance reduction, not
# distortion); (ii) boundary reflection at -1 and 1; (iii) bandwidth
# widened to 1.75x Sheather-Jones for visually smooth curves.
BW_INFLATE <- 1.75
kde_rho <- function(rho) {
  rho_s <- c(rho, -rho)
  h <- BW_INFLATE * stats::bw.SJ(rho_s)
  aug <- c(rho_s, -2 - rho_s, 2 - rho_s)   # reflect at both boundaries
  d <- density(aug, bw = h, from = -1, to = 1, n = 512)
  list(x = d$x, y = 3 * d$y, h = h)
}

densB_rows <- list()
for (key in names(sel)) {
  r <- sel[[key]]
  d <- kde_rho(r$rho)
  densB_rows[[length(densB_rows) + 1L]] <- data.frame(
    delta_key = key, delta = r$delta, x = d$x, density = d$y)
  cat(sprintf("Panel B: key %s uses delta = %.3f (n = %d, bw = %.3f)\n",
              key, r$delta, length(r$rho), d$h))
}
densB <- do.call(rbind, densB_rows)
densB$delta_key <- factor(densB$delta_key, levels = DELTA_LEVELS)

pB <- ggplot(densB, aes(x = x, y = density, colour = delta_key,
                        group = delta_key)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60",
             linewidth = LW_GUIDE) +
  geom_line(linewidth = LW_MAIN, show.legend = c(colour = TRUE)) +
  scale_delta() +
  guides(colour = guide_legend(order = 2, nrow = 1)) +
  labs(x = expression(rho[ij] ~ "given" ~ gamma[ij] == 1),
       y = "density") +
  coord_cartesian(xlim = c(-1, 1)) +
  theme_fig()

fig <- (pA | pB) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom", legend.box = "vertical")

save_fig(OUT_PDF, fig, height = fig_height(n_rows = 1L, n_legend_rows = 2L))
