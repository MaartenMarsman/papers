# R/scripts/fig_marginal_offdiag.R
#
# Figure for Section 3.4 "Conditional and marginal structure of single
# entries": marginal density of theta_{ij} | gamma_{ij} = 1 and marginal
# variance vs q. Reads R/results/marginal_offdiag.rds.
#
# Panel A uses the conditional marginal density estimator (CMDE) of Chen
# (1994) -- a Rao-Blackwellised marginal estimator: for each draw we have
# the partial-residual scalars (a_i, c_i, m_i), so the conditional density
# of theta_{ij} given the rest is closed-form, and the CMDE averages these
# conditional densities over draws to give a much smoother estimate of the
# marginal than a kernel density estimate.

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
IN_RDS  <- file.path(PROJECT_ROOT, "R", "results", "marginal_offdiag.rds")
OUT_PDF <- file.path(PROJECT_ROOT, "figures", "marginal_offdiag.pdf")

obj <- readRDS(IN_RDS)
tab <- obj$summary
results <- obj$results
SIGMA <- obj$params$sigma

# Conditional density p(theta_{ij} = y | a, c, m, delta, sigma)
# ~ (1 - ((y - m)^2 / (a*c)))^delta * exp(-y^2/(2 sigma^2))
# on |y - m| < sqrt(a*c). Returns unnormalised density.
cond_density_unn <- function(y, a, c, m, delta, sigma) {
  r2 <- a * c
  u2 <- (y - m)^2 / r2
  inside <- (u2 < 1) & (r2 > 0)
  out <- numeric(length(y) * max(length(a), 1L))
  if (length(a) == 1L) {
    out <- ifelse(inside,
                  (1 - u2)^delta * exp(-y^2 / (2 * sigma^2)),
                  0)
  } else {
    # vectorised: each row of u2 is one (y, a, c, m) combination
    out <- ifelse(inside,
                  (1 - u2)^delta * exp(-y^2 / (2 * sigma^2)),
                  0)
  }
  out
}

# Per-sample normaliser Z_i = int_{m-r}^{m+r} (1 - u^2)^delta * exp(-y^2 / (2 sigma^2)) dy.
# Use a fixed-grid trapezoidal rule per sample for speed.
sample_normalisers <- function(a_vec, c_vec, m_vec, delta, sigma,
                                n_grid = 80L) {
  n <- length(a_vec)
  r <- sqrt(a_vec * c_vec)
  Zs <- numeric(n)
  for (i in seq_len(n)) {
    if (!is.finite(r[i]) || r[i] <= 0) { Zs[i] <- NA_real_; next }
    ys <- seq(m_vec[i] - r[i], m_vec[i] + r[i], length.out = n_grid)
    vals <- cond_density_unn(ys, a_vec[i], c_vec[i], m_vec[i], delta, sigma)
    Zs[i] <- sum((vals[-1] + vals[-n_grid]) / 2) * (ys[2] - ys[1])
  }
  Zs
}

# CMDE: hat_p(x) = mean over samples of p(x | a_i, c_i, m_i) / Z_i restricted
# to samples whose support contains x.
cmde_density <- function(x_grid, a_vec, c_vec, m_vec, delta, sigma,
                         n_norm_grid = 80L) {
  Zs <- sample_normalisers(a_vec, c_vec, m_vec, delta, sigma, n_norm_grid)
  ok <- is.finite(Zs) & Zs > 0
  a_vec <- a_vec[ok]; c_vec <- c_vec[ok]; m_vec <- m_vec[ok]; Zs <- Zs[ok]
  r <- sqrt(a_vec * c_vec)
  n <- length(a_vec)
  dens <- numeric(length(x_grid))
  for (j in seq_along(x_grid)) {
    x <- x_grid[j]
    contribs <- cond_density_unn(rep(x, n), a_vec, c_vec, m_vec, delta, sigma) / Zs
    dens[j] <- mean(contribs)
  }
  dens
}

# Panel A: CMDE-estimated marginal density at delta = 1, across q.
panel_delta <- 1
x_grid <- seq(-3, 3, length.out = 201L)

dens_rows <- list()
for (r in results) {
  if (r$delta != panel_delta) next
  if (length(r$k12_on_edge) < 30L) next
  if (is.null(r$a_on) || is.null(r$c_on) || is.null(r$m_on)) {
    # Fall back to kernel density if (a, c, m) absent.
    kd <- density(r$k12_on_edge, from = min(x_grid), to = max(x_grid), n = length(x_grid))
    dens_rows[[length(dens_rows) + 1L]] <- data.frame(
      q = r$q, x = x_grid, density = approx(kd$x, kd$y, xout = x_grid, rule = 2)$y,
      method = "KDE", stringsAsFactors = FALSE)
  } else {
    dens <- cmde_density(x_grid, r$a_on, r$c_on, r$m_on, panel_delta, SIGMA)
    dens_rows[[length(dens_rows) + 1L]] <- data.frame(
      q = r$q, x = x_grid, density = dens,
      method = "CMDE", stringsAsFactors = FALSE)
  }
}
dens_df <- do.call(rbind, dens_rows)
dens_df$q_factor <- factor(dens_df$q, levels = sort(unique(dens_df$q)))

palette_q <- c(`5` = "#a6cee3", `10` = "#1f78b4", `15` = "#b2df8a",
               `20` = "#33a02c", `30` = "#e31a1c")

theme_paper <- theme_minimal(base_size = 9) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = 8),
    panel.grid.minor   = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 9, face = "plain"),
    plot.title         = element_text(size = 9, face = "plain")
  )

pA_title <- if (any(dens_df$method == "CMDE")) {
  sprintf("A. Marginal density at delta = %g (Chen 1994 CMDE)", panel_delta)
} else {
  sprintf("A. Marginal density at delta = %g (KDE)", panel_delta)
}

pA <- ggplot(dens_df, aes(x = x, y = density, colour = q_factor)) +
  geom_line(linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = palette_q, name = "q") +
  labs(x = expression(theta[ij] * "  given  " * gamma[ij] == 1),
       y = "density",
       title = pA_title) +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  theme_paper

# Panel B: marginal variance vs q for delta in {0, 1, 2}, with bootstrap 95% CI.
boot_var_ci <- function(x, n_boot = 1000L, seed = 1L) {
  set.seed(seed)
  if (length(x) < 10L) return(c(NA_real_, NA_real_))
  vs <- replicate(n_boot, var(sample(x, length(x), replace = TRUE)))
  as.numeric(quantile(vs, c(0.025, 0.975)))
}

var_df <- do.call(rbind, lapply(results, function(r) {
  on <- r$k12_on_edge
  if (length(on) < 10L) return(NULL)
  ci <- boot_var_ci(on, n_boot = 1000L, seed = as.integer(100 * r$q + 10 * r$delta + 1L))
  data.frame(q = r$q, delta = r$delta, var = var(on),
             var_lo = ci[1], var_hi = ci[2],
             n = length(on), stringsAsFactors = FALSE)
}))
var_df$delta_factor <- factor(sprintf("delta=%g", var_df$delta),
                               levels = c("delta=0", "delta=1", "delta=2"))

palette_delta <- c("delta=0" = "#1f78b4",
                   "delta=1" = "#33a02c",
                   "delta=2" = "#e31a1c")
labels_delta <- c(
  "delta=0" = expression(delta == 0),
  "delta=1" = expression(delta == 1),
  "delta=2" = expression(delta == 2)
)

pB <- ggplot(var_df,
             aes(x = q, y = var, colour = delta_factor, group = delta_factor)) +
  geom_ribbon(aes(ymin = var_lo, ymax = var_hi, fill = delta_factor),
              alpha = 0.15, colour = NA) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(5, 10, 15, 20, 30)) +
  scale_y_log10() +
  scale_colour_manual(values = palette_delta, name = NULL, labels = labels_delta) +
  scale_fill_manual(values = palette_delta, name = NULL, labels = labels_delta) +
  labs(x = "q",
       y = expression("Var(" * theta[ij] * " | " * gamma[ij] == 1 * ")"),
       title = "B. Marginal variance vs q (95% bootstrap CI)") +
  theme_paper

fig <- gridExtra::grid.arrange(pA, pB, ncol = 2)
ggplot2::ggsave(OUT_PDF, plot = fig, width = 9.2, height = 3.7, units = "in")
cat(sprintf("Wrote %s\n", OUT_PDF))

cat("\nMarginal variance by (q, delta):\n")
print(tab[, c("q", "delta", "n_on_edge", "mean", "var")], row.names = FALSE)
