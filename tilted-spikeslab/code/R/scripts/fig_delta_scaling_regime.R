# Regime-explicit delta-scaling figure + inversion analysis.
# Consumes R/results/delta_scaling_regime_sweep.rds (delta_scaling_regime_sweep.R).
# Produces figures/delta_scaling_regime.pdf and
# R/results/delta_scaling_regime_inversion.rds, and prints the numbers used in
# sec:choosing-delta (criterion-1 closed form, fidelity at the geometric
# delta*, censoring report).

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
source(file.path(PROJECT_ROOT, "R", "scripts", "fig_style.R"))

res <- readRDS(file.path(PROJECT_ROOT, "R", "results", "delta_scaling_regime_sweep.rds"))

## ---- aggregate over graph replicates ----
agg_key <- interaction(res$q, res$axis, res$level, res$eta, res$delta, drop = TRUE)
agg <- do.call(rbind, lapply(split(res, agg_key), function(d) {
  data.frame(q = d$q[1], axis = d$axis[1], level = d$level[1], eta = d$eta[1],
             delta = d$delta[1], lmin_med = mean(d$lmin_med),
             var_ratio = mean(d$var_ratio), rho_sd = mean(d$rho_sd),
             dbar = mean(d$dbar_realized), m = mean(d$m),
             split_gap = mean(abs(d$lmin_second_half - d$lmin_first_half)))
}))
rownames(agg) <- NULL

cat(sprintf("max split-half gap in median lambda_min: %.4f\n", max(agg$split_gap)))

## ---- invert lambda_min median for delta*(target) ----
invert_delta <- function(delta, y, target) {
  o <- order(delta)
  delta <- delta[o]; y <- y[o]
  if (max(y) < target) return(NA_real_)      # censored: not reached at delta = 24
  if (y[1] >= target) return(delta[1])
  i <- which(y >= target)[1]
  delta[i - 1] + (target - y[i - 1]) / (y[i] - y[i - 1]) * (delta[i] - delta[i - 1])
}

targets <- c(0.25, 0.5)
inv_key <- interaction(agg$q, agg$axis, agg$level, agg$eta, drop = TRUE)
inv <- do.call(rbind, lapply(split(agg, inv_key), function(d) {
  do.call(rbind, lapply(targets, function(tg) {
    data.frame(q = d$q[1], axis = d$axis[1], level = d$level[1], eta = d$eta[1],
               target = tg, dbar = mean(d$dbar), m = mean(d$m),
               delta_star = invert_delta(d$delta, d$lmin_med, tg))
  }))
}))
rownames(inv) <- NULL
saveRDS(list(agg = agg, inversion = inv),
        file.path(PROJECT_ROOT, "R", "results", "delta_scaling_regime_inversion.rds"))

cat("\n== delta* (target median lambda_min = 0.5), by regime ==\n")
print(subset(inv, target == 0.5)[order(subset(inv, target == 0.5)$axis,
                                       subset(inv, target == 0.5)$level,
                                       subset(inv, target == 0.5)$eta,
                                       subset(inv, target == 0.5)$q), ],
      row.names = FALSE)
cat(sprintf("\ncensored cells (delta* > 24): %d of %d\n",
            sum(is.na(inv$delta_star)), nrow(inv)))

## ---- criterion 1 (per-edge constant), closed form ----
c0_fun <- function(delta, eta) {
  # c(0) = E[(1 + 2 eta^2 / T)^{-1/2}], T ~ Gamma(delta + 1, 1); log-scale
  # normalized integrand with t = (delta + 1) s so the Gamma peak sits at s ~ 1
  # for every delta
  s_int <- function(s) {
    t <- (delta + 1) * s
    exp(delta * log(t) - t - lgamma(delta + 1) + log(delta + 1)) /
      sqrt(1 + 2 * eta^2 / t)
  }
  stats::integrate(s_int, 0, Inf, rel.tol = 1e-11)$value
}
delta_star_c0 <- function(eta, kappa) {
  stats::uniroot(function(d) c0_fun(d, eta) - (1 - kappa), c(0.001, 400))$root
}
cat("\n== criterion 1: delta* for c(0) >= 1 - kappa ==\n")
for (eta in c(0.5, 1, 2)) for (kappa in c(0.1, 0.05)) {
  cat(sprintf("eta = %.1f  kappa = %.2f  delta* = %.1f  (asymptote eta^2/kappa = %.1f)\n",
              eta, kappa, delta_star_c0(eta, kappa), eta^2 / kappa))
}

## ---- fidelity at the geometric delta* (criterion 3 readout) ----
cat("\n== rho SD and variance ratio at the delta grid point nearest delta*(0.5), eta = 1 ==\n")
sub <- subset(inv, target == 0.5 & eta == 1 & !is.na(delta_star))
for (i in seq_len(nrow(sub))) {
  d <- subset(agg, q == sub$q[i] & axis == sub$axis[i] & level == sub$level[i] & eta == 1)
  j <- which.min(abs(d$delta - sub$delta_star[i]))
  cat(sprintf("q = %2d  %s = %.2f: delta* = %5.2f  ->  rho SD = %.3f, var ratio = %.3f\n",
              sub$q[i], sub$axis[i], sub$level[i], sub$delta_star[i],
              d$rho_sd[j], d$var_ratio[j]))
}

## ---- figure: delta* regime map ----
inv_plot <- subset(inv, target == 0.5)
inv_plot$level_f <- factor(inv_plot$level)
inv_plot$eta_f <- factor(inv_plot$eta,
                         labels = c("eta == 0.5", "eta == 1", "eta == 2"))
ref <- expand.grid(q = seq(10, 50, by = 1),
                   eta_f = levels(inv_plot$eta_f), stringsAsFactors = FALSE)
ref$delta_ref <- 0.5 * log(ref$q)

PAL_LEVEL <- c(OKABE_ITO[["blue"]], OKABE_ITO[["orange"]], OKABE_ITO[["vermillion"]])

panel_regime <- function(dat, axis_name, legend_title, level_labels) {
  d <- subset(dat, axis == axis_name)
  d$level_f <- droplevels(d$level_f)
  ggplot(d, aes(x = q, y = delta_star, colour = level_f, group = level_f)) +
    geom_line(data = ref, aes(x = q, y = delta_ref),
              inherit.aes = FALSE, linetype = "22",
              colour = OKABE_ITO[["grey"]], linewidth = LW_REF) +
    geom_line(linewidth = LW_MAIN, na.rm = TRUE) +
    geom_point(size = PT_MAIN, na.rm = TRUE) +
    facet_wrap(~eta_f, nrow = 1, labeller = label_parsed) +
    scale_colour_manual(values = stats::setNames(PAL_LEVEL, levels(d$level_f)),
                        labels = level_labels, name = legend_title) +
    scale_x_continuous(breaks = c(10, 20, 50)) +
    labs(x = LAB_Q, y = expression("required tilt"~delta^"*")) +
    theme_fig()
}

pA <- panel_regime(inv_plot, "degree",
                   expression("expected degree"~bar(d)),
                   c("2", "4", "8"))
pB <- panel_regime(inv_plot, "density",
                   expression("edge density"~rho),
                   c("0.05", "0.10", "0.20"))

fig <- (pA / pB) + plot_annotation(tag_levels = "A")
save_fig(file.path(PROJECT_ROOT, "figures", "delta_scaling_regime.pdf"),
         fig, height = fig_height(n_rows = 2L, n_legend_rows = 2L))
