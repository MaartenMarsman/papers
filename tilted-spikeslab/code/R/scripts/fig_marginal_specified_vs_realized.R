# R/scripts/fig_marginal_specified_vs_realized.R
#
# "The slab you specify vs the slab you get", at the MARGINAL level. The
# specified slab is N(0, sigma^2). The realized marginal of k_{ij} | gamma=1 is
# the average of the compactly supported conditional slabs over the rest-block
# geometry; it is estimated by the Chen (1994) conditional marginal density
# estimator (CMDE), reusing (a, c, m) from R/results/marginal_offdiag.rds.
#
# Adds: (i) the specified N(0, sigma^2) reference on the density panel;
#       (ii) quantitative gap vs (q, delta): realized/specified variance ratio,
#            total variation, and KL(realized || specified).
#
# Style:  R/scripts/fig_style.R (shared theme, palette, export contract)
# Data:   R/results/marginal_offdiag.rds (cached; do not recompute)
# Output: figures/marginal_specified_vs_realized.pdf
#         tables/marginal_specified_vs_realized.tex

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
source(file.path(PROJECT_ROOT, "R", "scripts", "fig_style.R"))

IN_RDS  <- file.path(PROJECT_ROOT, "R", "results", "marginal_offdiag.rds")
OUT_PDF <- file.path(PROJECT_ROOT, "figures", "marginal_specified_vs_realized.pdf")
OUT_TEX <- file.path(PROJECT_ROOT, "tables", "marginal_specified_vs_realized.tex")

obj <- readRDS(IN_RDS)
results <- obj$results
SIGMA <- obj$params$sigma

# --- CMDE machinery (from fig_marginal_offdiag.R) --------------------------
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

dens_rows <- list(); gap_rows <- list()
for (r in results) {
  if (is.null(r$a_on) || length(r$k12_on_edge) < 50L) next
  dens <- cmde_density(x_grid, r$a_on, r$c_on, r$m_on, r$delta, SIGMA)
  area <- trapz(dens, dx); if (area > 0) dens <- dens / area  # guard

  dens_rows[[length(dens_rows) + 1L]] <- data.frame(
    q = r$q, delta = r$delta, x = x_grid, density = dens)

  # quantitative gap between realized (CMDE) and specified (Gaussian)
  var_real <- var(r$k12_on_edge)
  tv <- 0.5 * trapz(abs(dens - specified), dx)
  pos <- dens > 0 & specified > 0
  kl <- trapz(ifelse(pos, dens * log(dens / specified), 0), dx)
  gap_rows[[length(gap_rows) + 1L]] <- data.frame(
    q = r$q, delta = r$delta,
    var_ratio = var_real / SIGMA^2,
    tv = tv, kl = kl, n = length(r$k12_on_edge))
}
dens_df <- do.call(rbind, dens_rows)
gap_df  <- do.call(rbind, gap_rows)
gap_df  <- gap_df[order(gap_df$delta, gap_df$q), ]

# --- density overlay panel (facet by delta, colour by q) -------------------
q_levels <- sort(unique(dens_df$q))
dens_df$q_factor <- factor(dens_df$q, levels = q_levels)
dens_df$delta_lab <- factor(sprintf("delta == %g", dens_df$delta),
                            levels = sprintf("delta == %g", sort(unique(dens_df$delta))))
spec_df <- data.frame(x = x_grid, density = specified)

p <- ggplot(dens_df, aes(x = x, y = density, colour = q_factor)) +
  geom_line(data = spec_df,
            aes(x = x, y = density, linetype = "specified"),
            inherit.aes = FALSE, colour = PAL_CURVE[["specified"]],
            linewidth = LW_REF) +
  geom_line(linewidth = LW_MAIN) +
  facet_wrap(~ delta_lab, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = pal_q(q_levels),
                      name = expression(italic(q))) +
  scale_linetype_manual(values = c(specified = "dashed"), name = NULL,
                        labels = expression("specified slab" ~ N(0, sigma^2))) +
  guides(linetype = guide_legend(order = 1),
         colour = guide_legend(order = 2, nrow = 1)) +
  labs(x = expression(k[ij] ~ "given" ~ gamma[ij] == 1),
       y = "density") +
  coord_cartesian(xlim = c(-3, 3)) +
  theme_fig()

save_fig(OUT_PDF, p, height = fig_height(n_rows = 1L, n_legend_rows = 1L) + 0.15)

# --- LaTeX table of the gap ------------------------------------------------
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
lines <- c(
  "\\begin{tabular}{rrrrr}",
  "\\hline",
  "$q$ & $\\delta$ & Var ratio & TV & KL \\\\",
  "\\hline",
  apply(gap_df, 1, function(row) sprintf("%d & %s & %s & %s & %s \\\\",
        as.integer(row["q"]), fmt(row["delta"], 2),
        fmt(row["var_ratio"]), fmt(row["tv"]), fmt(row["kl"]))),
  "\\hline",
  "\\end{tabular}"
)
writeLines(lines, OUT_TEX)
cat(sprintf("Wrote %s\n", OUT_TEX))

cat("\nSpecified-vs-realized gap (var ratio = realized/sigma^2):\n")
print(gap_df, row.names = FALSE)
