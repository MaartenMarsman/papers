# 23-run-ppc.R
# Posterior predictive check for the worked example (Robustness section):
# simulate() draws one replicated dataset per posterior draw, so parameter
# and structure uncertainty propagate into the replicates.
#
# The manuscript shows only the sum-score check: a global
# statistic no parameter targets, interpretable without estimator theory.
# The margin and cross-product statistics computed below are diagnostics
# whose reading depends on the pseudolikelihood's anchoring; they feed
# the full ladder in R/30-ppc-diagnostic-ladder.R and stay out of the figure.
#
# Reads tables/fit-wenchuan.rds (from R/01-worked-example.R).
# Outputs: figures/fig-ppc.pdf (sum-score panel only), tables/ppc-coverage.csv

library(bgms)

set.seed(123)

fit <- readRDS("tables/fit-wenchuan.rds")
data("Wenchuan", package = "bgms")
wenchuan <- Wenchuan[stats::complete.cases(Wenchuan), ]

n <- nrow(wenchuan)
p <- ncol(wenchuan)
# 500 posterior draws, one replicated dataset of n rows per draw.
n_draws <- 500

t0 <- Sys.time()
reps <- simulate(fit, nsim = n, method = "posterior-sample",
                 ndraws = n_draws, seed = 123, display_progress = "none")
cat("simulate:", round(as.numeric(Sys.time() - t0, units = "secs")), "s,",
    length(reps), "replicated datasets\n")

# --- Statistic 1: item margins (category proportions per variable) -------------

categories <- sort(unique(as.vector(as.matrix(wenchuan))))
margins <- function(x) {
  sapply(seq_len(ncol(x)), function(j)
    sapply(categories, function(k) mean(x[, j] == k)))
}
obs_margin <- margins(as.matrix(wenchuan))
rep_margin <- vapply(reps, margins, obs_margin)
margin_mean <- apply(rep_margin, c(1, 2), mean)
margin_lo <- apply(rep_margin, c(1, 2), stats::quantile, 0.025)
margin_hi <- apply(rep_margin, c(1, 2), stats::quantile, 0.975)
margin_cover <- mean(obs_margin >= margin_lo & obs_margin <= margin_hi)
cat(sprintf("margins: %.0f of %d category proportions inside 95%% intervals (%.1f%%)\n",
            margin_cover * length(obs_margin), length(obs_margin),
            100 * margin_cover))

# Pattern layer of the margins: the per-variable mean response
obs_mu <- colMeans(wenchuan)
rep_mu <- vapply(reps, colMeans, obs_mu)
mu_median <- apply(rep_mu, 1, stats::median)
cat(sprintf("variable means above predictive median: %d of %d; mean gap %.3f\n",
            sum(obs_mu > mu_median), p, mean(obs_mu - rowMeans(rep_mu))))

# --- Statistic 2: pairwise sufficient statistics (mean cross-products) ---------

lt <- lower.tri(matrix(0, p, p))
crossmoment <- function(x) {
  x <- as.matrix(x)
  (crossprod(x) / nrow(x))[lt]
}
obs_cm <- crossmoment(wenchuan)
rep_cm <- vapply(reps, crossmoment, obs_cm)
cm_mean <- rowMeans(rep_cm)
cm_median <- apply(rep_cm, 1, stats::median)
cm_lo <- apply(rep_cm, 1, stats::quantile, 0.025)
cm_hi <- apply(rep_cm, 1, stats::quantile, 0.975)
cm_cover <- obs_cm >= cm_lo & obs_cm <= cm_hi
cat(sprintf("pairwise cross-products: %d of %d inside 95%% intervals\n",
            sum(cm_cover), length(obs_cm)))
cat(sprintf("pairwise cross-products above predictive median: %d of %d; mean gap %.3f\n",
            sum(obs_cm > cm_median), length(obs_cm), mean(obs_cm - cm_mean)))

# Decomposition of the cross-product gap into a product-of-means part
# (margins) and a centered-covariance part (dependence), so the two sources
# of misfit can be read separately.
cova <- function(x) {
  x <- as.matrix(x)
  stats::cov(x)[lt] * (nrow(x) - 1) / nrow(x)
}
obs_cv <- cova(wenchuan)
rep_cv <- vapply(reps, cova, obs_cv)
cv_median <- apply(rep_cv, 1, stats::median)
pmpart <- function(m) outer(m, m)[lt]
obs_pm <- pmpart(obs_mu)
rep_pm <- vapply(seq_along(reps), function(i) pmpart(rep_mu[, i]), obs_pm)
cat(sprintf("  of which product-of-means part: mean gap %.3f\n",
            mean(obs_pm - rowMeans(rep_pm))))
cat(sprintf("  of which covariance part: %d of %d above predictive median; mean gap %.3f\n",
            sum(obs_cv > cv_median), length(obs_cv), mean(obs_cv - rowMeans(rep_cv))))

# --- Statistic 3: sum-score distribution ---------------------------------------

score_range <- (p * min(categories)):(p * max(categories))
score_dist <- function(x) {
  s <- factor(rowSums(x), levels = score_range)
  as.numeric(table(s)) / nrow(x)
}
obs_score <- score_dist(as.matrix(wenchuan))
rep_score <- vapply(reps, score_dist, obs_score)
score_lo <- apply(rep_score, 1, stats::quantile, 0.025)
score_hi <- apply(rep_score, 1, stats::quantile, 0.975)
score_cover <- obs_score >= score_lo & obs_score <= score_hi
cat(sprintf("sum-score distribution: %d of %d score values inside 95%% intervals\n",
            sum(score_cover), length(obs_score)))

utils::write.csv(
  data.frame(
    statistic = c("category proportions", "pairwise cross-products", "sum-score values"),
    total = c(length(obs_margin), length(obs_cm), length(obs_score)),
    covered = c(round(margin_cover * length(obs_margin)), sum(cm_cover),
                sum(score_cover))
  ),
  "tables/ppc-coverage.csv", row.names = FALSE
)

# --- Figure: two panels, styled like fig-spike-slab ----------------------------

slab_fill <- "#C9DCEE"
slab_line <- "#5B93C4"
dark_blue <- "#1F3A63"

score_df <- data.frame(score = score_range, obs = obs_score,
                       lo = score_lo, hi = score_hi)
panel_score <- ggplot2::ggplot(score_df, ggplot2::aes(x = score)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                       fill = slab_fill, color = slab_line, linewidth = 0.4) +
  ggplot2::geom_point(ggplot2::aes(y = obs), color = dark_blue, size = 1.1) +
  ggplot2::labs(
    x = paste("Sum score across the", p, "symptoms"),
    y = "Proportion of participants"
  )

theme_ppc <- ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5))

ggplot2::ggsave("figures/fig-ppc.pdf", panel_score + theme_ppc,
                width = 6.5, height = 3.4)
