# 30-ppc-diagnostic-ladder.R
# The full diagnostic ladder behind the (deliberately simple) predictive
# check of the robustness section. None of this appears in the manuscript:
# the interpretation requires estimator theory the tutorial cannot assume.
#
# Ladder (each statistic is compared with its posterior predictive
# distribution from 500 replicated datasets):
#   1. sufficient statistics: variable means and pairwise cross-products,
#   2. decomposition of the cross-product gap into a product-of-means part
#      and a covariance part,
#   3. per-category tilt and the minimal-responder corner,
#   4. conditionals-at-observed-data versus the implied joint,
#   5. third central cross-moments and sum-score skewness,
#   6. sum-score variance and kurtosis.
#
# Reads tables/fit-wenchuan.rds (from R/01-worked-example.R).
# Outputs: console numbers + figures/fig-ppc-third-order.pdf (the
# third-order panel).

library(bgms)

set.seed(123)
fit <- readRDS("tables/fit-wenchuan.rds")
data("Wenchuan", package = "bgms")
wen <- as.matrix(Wenchuan[stats::complete.cases(Wenchuan), ])
n <- nrow(wen)
p <- ncol(wen)
lt <- lower.tri(matrix(0, p, p))

reps <- simulate(fit, nsim = n, method = "posterior-sample",
                 ndraws = 500, seed = 123, display_progress = "none")

# --- 1. Sufficient statistics --------------------------------------------------

obs_mu <- colMeans(wen)
rep_mu <- vapply(reps, colMeans, obs_mu)
cat(sprintf("means: %d of %d below predictive median; mean gap %+.3f\n",
            sum(obs_mu < apply(rep_mu, 1, stats::median)), p,
            mean(obs_mu - rowMeans(rep_mu))))

crossm <- function(x) (crossprod(as.matrix(x)) / nrow(x))[lt]
obs_cm <- crossm(wen)
rep_cm <- vapply(reps, crossm, obs_cm)
pred_cm <- rowMeans(rep_cm)
gap <- obs_cm - pred_cm
cat(sprintf("cross-products: cor(obs, pred) %.4f; gap %+.3f (sd %.3f); %d of %d below median\n",
            stats::cor(obs_cm, pred_cm), mean(gap), stats::sd(gap),
            sum(obs_cm < apply(rep_cm, 1, stats::median)), length(obs_cm)))

# --- 2. Decomposition: product-of-means part + covariance part -----------------

cova <- function(x) {
  x <- as.matrix(x)
  stats::cov(x)[lt] * (nrow(x) - 1) / nrow(x)
}
obs_cv <- cova(wen)
rep_cv <- vapply(reps, cova, obs_cv)
pmpart <- function(m) outer(m, m)[lt]
obs_pm <- pmpart(obs_mu)
rep_pm <- vapply(seq_along(reps), function(i) pmpart(rep_mu[, i]), obs_pm)
cat(sprintf("  product-of-means part: %+.3f; covariance part: %+.3f (%d of %d above median)\n",
            mean(obs_pm - rowMeans(rep_pm)), mean(obs_cv - rowMeans(rep_cv)),
            sum(obs_cv > apply(rep_cv, 1, stats::median)), length(obs_cv)))
musum <- outer(obs_mu, obs_mu, `+`)[lt]
cat(sprintf("  gap ~ (mu_i + mu_j): fitted slope %.3f against -delta = %.3f\n",
            stats::coef(stats::lm(gap ~ musum))[2],
            -mean(rowMeans(rep_mu) - obs_mu)))

# --- 3. Category tilt and the minimal-responder corner -------------------------

for (k in sort(unique(as.vector(wen)))) {
  ob <- colMeans(wen == k)
  rp <- vapply(reps, function(x) colMeans(as.matrix(x) == k), ob)
  cat(sprintf("  category %d: observed above predictive median on %2d of %d variables (gap %+.4f)\n",
              k, sum(ob > apply(rp, 1, stats::median)), p, mean(ob - rowMeans(rp))))
}
sums <- rowSums(wen)
rsums <- vapply(reps, function(x) rowSums(as.matrix(x)), sums)
# Minimal responders: sum score at most five points above the minimum of 17
# (17 items coded 1 to 5); a fixed cutoff chosen for this check.
ob <- mean(sums <= 22)
rp <- colMeans(rsums <= 22)
cat(sprintf("minimal responders (sum <= 22): observed %.3f, replicates %.3f [%.3f, %.3f]\n",
            ob, mean(rp), stats::quantile(rp, 0.025), stats::quantile(rp, 0.975)))

# --- 4. Conditionals at observed data versus the implied joint -----------------

pr <- predict(fit, newdata = as.data.frame(wen), type = "probabilities")
cond_mu <- vapply(seq_len(p), function(j) {
  pj <- as.matrix(pr[[j]])
  cats <- suppressWarnings(as.numeric(colnames(pj)))
  if (any(is.na(cats))) cats <- sort(unique(wen[, j]))
  mean(pj %*% cats)
}, numeric(1))
cat(sprintf("conditionals at observed data: mean gap %+.4f (max |gap| %.4f)\n",
            mean(obs_mu - cond_mu), max(abs(obs_mu - cond_mu))))
cat(sprintf("implied joint:                 mean gap %+.4f\n",
            mean(obs_mu - rowMeans(rep_mu))))

# --- 5. Third order ------------------------------------------------------------

tri <- utils::combn(p, 3)
m3 <- function(x) {
  z <- scale(as.matrix(x), center = TRUE, scale = FALSE)
  apply(tri, 2, function(t) mean(z[, t[1]] * z[, t[2]] * z[, t[3]]))
}
obs3 <- m3(wen)
rep3 <- vapply(reps, m3, obs3)
cat(sprintf("third cross-moments: %d of %d above predictive median; %d outside 95%% intervals\n",
            sum(obs3 > apply(rep3, 1, stats::median)), ncol(tri),
            sum(obs3 < apply(rep3, 1, stats::quantile, 0.025) |
                obs3 > apply(rep3, 1, stats::quantile, 0.975))))
skew <- function(s) mean((s - mean(s))^3) / stats::sd(s)^3
rp <- vapply(seq_along(reps), function(i) skew(rsums[, i]), 1)
cat(sprintf("sum-score skewness: observed %.3f, predictive %.3f [%.3f, %.3f]\n",
            skew(sums), mean(rp), stats::quantile(rp, 0.025), stats::quantile(rp, 0.975)))

# --- 6. Fourth order -----------------------------------------------------------

rpv <- vapply(seq_along(reps), function(i) stats::var(rsums[, i]), 1)
cat(sprintf("sum-score variance: observed %.1f, predictive %.1f [%.1f, %.1f]; gap %.1f vs 2*sum(cov gap) = %.1f\n",
            stats::var(sums), mean(rpv), stats::quantile(rpv, 0.025),
            stats::quantile(rpv, 0.975), stats::var(sums) - mean(rpv),
            2 * sum(obs_cv - rowMeans(rep_cv))))
kurt <- function(s) mean((s - mean(s))^4) / stats::sd(s)^4
rpk <- vapply(seq_along(reps), function(i) kurt(rsums[, i]), 1)
cat(sprintf("sum-score kurtosis: observed %.3f, predictive %.3f [%.3f, %.3f]\n",
            kurt(sums), mean(rpk), stats::quantile(rpk, 0.025), stats::quantile(rpk, 0.975)))

# --- Third-order figure (same style as the other check figures) ---------------

slab_line <- "#5B93C4"
dark_blue <- "#1F3A63"
lo3 <- apply(rep3, 1, stats::quantile, 0.025)
hi3 <- apply(rep3, 1, stats::quantile, 0.975)
fig3 <- ggplot2::ggplot(
  data.frame(obs = obs3, lo = lo3, hi = hi3, mean = rowMeans(rep3),
             miss = obs3 < lo3 | obs3 > hi3),
  ggplot2::aes(x = mean, y = obs)
) +
  ggplot2::geom_abline(intercept = 0, slope = 1, color = "grey60") +
  ggplot2::geom_segment(ggplot2::aes(x = lo, xend = hi, yend = obs),
                        color = slab_line, alpha = 0.25) +
  ggplot2::geom_point(ggplot2::aes(color = miss), size = 1.0) +
  ggplot2::scale_color_manual(values = c(`FALSE` = dark_blue, `TRUE` = "#A66A1E"),
                              guide = "none") +
  ggplot2::labs(
    x = "Posterior predictive third central cross-moment (95% interval)",
    y = "Observed third central cross-moment",
    title = sprintf("Third-order cross-moments (%d triples): the misfit skips this order",
                    ncol(tri))
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5))
ggplot2::ggsave("figures/fig-ppc-third-order.pdf", fig3, width = 6.5, height = 5)
cat("wrote figures/fig-ppc-third-order.pdf\n")
