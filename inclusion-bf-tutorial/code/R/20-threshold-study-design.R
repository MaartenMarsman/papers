# Known-truth calibration study, step 1: calibration constants, generating
# networks, data sets, and the fit manifest.
# Design: the simulation appendix of the manuscript.
#
# Three arms, each with its own calibration source:
#   ordinal (stage 1)  Wenchuan reference fit in tables/fit-wenchuan.rds
#   binary  (stage 2)  a bgm() fit of the ADHD data shipped with bgms
#   ggm     (stage 2)  documented, no shipped continuous data set exists
#
# Writes tables/threshold-study/datasets/*.rds (one per data set),
# tables/threshold-study/calibration.rds, and
# tables/threshold-study/manifest.csv (one row per fit, test and reference).

source("R/threshold-study-helpers.R")

paths <- study_paths()

num_variables <- 15L
num_categories_ordinal <- 4L    # 5-point ordinal: 4 categories above the base
densities <- c(low = 0.2, mid = 0.4, high = 0.6)
sample_sizes <- c(200L, 1000L)
data_seeds_stage1 <- c(1L, 2L)
data_seeds_stage2 <- 1L
sampler_seeds <- 1:5

# Test-fit settings: the two the tutorial discusses (Table 1 defaults and the
# long rerun). Chain count is fixed at 4 because split-Rhat, the agreement
# test, and the summed ESS are all functions of it, so thresholds calibrated
# at another chain count would not apply to readers following Table 1.
test_settings <- list(
  default = c(warmup = 2000L, iter = 2000L),
  long = c(warmup = 5000L, iter = 10000L)
)

# Reference runs: chain count is free here (nothing is graded against it and
# warmup parallelizes), so they buy precision with chains.
reference_setting <- c(warmup = 10000L, iter = 25000L)
reference_chains <- 12L
test_chains <- 4L

# Target scale for the nonzero interactions of the ORDINAL arm, in the
# data-preferred band of the prior sensitivity analysis (RMS of included
# interactions on Wenchuan: 0.18 [0.148, 0.218]). The calibration gives a
# Wenchuan pool RMS of 0.207, so this is a 3% adjustment that pins the arm to
# a fixed constant.
#
# The binary arm is NOT rescaled to this number. An interaction enters the
# linear predictor as theta * x, so a 0/1 variable needs roughly four times
# the coefficient of a 0-4 ordinal variable for the same coupling; the ADHD
# fit's conditional weights come out near 0.72 for that reason. Forcing them
# onto the ordinal scale would leave the binary arm's edges undetectable and
# would no longer be a fit-calibrated arm. Each arm keeps its own
# data-preferred scale.
target_weight_rms <- 0.20

# --- Calibration: ordinal arm from the Wenchuan reference fit -------------------

cat("=== ordinal calibration (Wenchuan) ===\n")

data("Wenchuan", package = "bgms")
wenchuan <- Wenchuan[stats::complete.cases(Wenchuan), ] - 1L

wenchuan_fit <- readRDS("tables/fit-wenchuan.rds")
wenchuan_summary <- summary(wenchuan_fit)
wenchuan_raw <- S7::prop(wenchuan_fit, "raw_samples")

# Conditional (slab) interaction strengths: the mean weight over the draws in
# which the edge is included. The unconditional column would shrink every
# weight toward zero by its inclusion probability.
conditional_weight <- local({
  n_chains <- as.integer(wenchuan_raw$nchains)
  num_edges <- length(wenchuan_raw$parameter_names$pairwise)
  totals <- numeric(num_edges)
  counts <- numeric(num_edges)
  for (ch in seq_len(n_chains)) {
    included <- wenchuan_raw$indicator[[ch]] == 1
    values <- wenchuan_raw$pairwise[[ch]]
    values[!included] <- 0
    totals <- totals + colSums(values)
    counts <- counts + colSums(included)
  }
  ifelse(counts > 0, totals / counts, NA_real_)
})

decisive <- wenchuan_summary$indicator$mean > 0.9 & !is.na(conditional_weight)
weight_pool_raw <- abs(conditional_weight[decisive])
p_negative_ordinal <- mean(conditional_weight[decisive] < 0)
weight_pool_rms_raw <- sqrt(mean(weight_pool_raw^2))
weight_pool_ordinal <- weight_pool_raw * target_weight_rms / weight_pool_rms_raw

cat(sprintf("  decisive edges: %d; conditional |weight| RMS %.3f -> %.3f\n",
            sum(decisive), weight_pool_rms_raw, target_weight_rms))
cat(sprintf("  negative sign share: %.3f\n", p_negative_ordinal))

# The study uses 15 variables and Wenchuan has 17, so a fixed subset of items
# supplies the thresholds and the marginal targets. Their own thresholds keep
# the simulated response distributions item-specific rather than uniform.
set.seed(20260730L)
ordinal_items <- sort(sample(ncol(wenchuan), num_variables))
wenchuan_thresholds <- S7::prop(wenchuan_fit,
                                "posterior_mean_main")[ordinal_items, ,
                                                       drop = FALSE]
wenchuan_targets <- t(apply(wenchuan[, ordinal_items], 2, function(v) {
  prop.table(table(factor(v, levels = 0:num_categories_ordinal)))
}))

rm(wenchuan_fit, wenchuan_raw)
invisible(gc())

# --- Calibration: binary arm from a bgm() fit of the ADHD data -----------------

cat("=== binary calibration (ADHD) ===\n")

data("ADHD", package = "bgms")
adhd <- ADHD[, setdiff(colnames(ADHD), grep("group", colnames(ADHD),
                                            ignore.case = TRUE, value = TRUE))]
adhd <- adhd[stats::complete.cases(adhd), ]

adhd_fit_path <- file.path(paths$root, "fit-adhd-calibration.rds")
if (file.exists(adhd_fit_path)) {
  adhd_calibration <- readRDS(adhd_fit_path)
} else {
  set.seed(123)
  adhd_fit <- bgms::bgm(
    x = adhd, variable_type = "ordinal",
    iter = 10000L, warmup = 5000L,
    edge_selection = TRUE, edge_prior = bgms::bernoulli_prior(0.5),
    interaction_prior = bgms::normal_prior(scale = 1),
    chains = 4L, cores = 4L, seed = 123,
    display_progress = "none", verbose = FALSE
  )
  adhd_summary <- summary(adhd_fit)
  adhd_raw <- S7::prop(adhd_fit, "raw_samples")
  adhd_conditional <- local({
    num_edges <- length(adhd_raw$parameter_names$pairwise)
    totals <- numeric(num_edges)
    counts <- numeric(num_edges)
    for (ch in seq_len(as.integer(adhd_raw$nchains))) {
      included <- adhd_raw$indicator[[ch]] == 1
      values <- adhd_raw$pairwise[[ch]]
      values[!included] <- 0
      totals <- totals + colSums(values)
      counts <- counts + colSums(included)
    }
    ifelse(counts > 0, totals / counts, NA_real_)
  })
  adhd_calibration <- list(
    conditional_weight = adhd_conditional,
    pip = adhd_summary$indicator$mean,
    thresholds = S7::prop(adhd_fit, "posterior_mean_main"),
    endorsement = colMeans(adhd)
  )
  saveRDS(adhd_calibration, adhd_fit_path)
  rm(adhd_fit, adhd_raw)
  invisible(gc())
}

adhd_decisive <- adhd_calibration$pip > 0.9 &
  !is.na(adhd_calibration$conditional_weight)
adhd_pool_raw <- abs(adhd_calibration$conditional_weight[adhd_decisive])
p_negative_binary <- mean(adhd_calibration$conditional_weight[adhd_decisive] < 0)
adhd_pool_rms_raw <- sqrt(mean(adhd_pool_raw^2))
weight_pool_binary <- adhd_pool_raw   # kept at the ADHD scale; see above

cat(sprintf("  decisive edges: %d; conditional |weight| RMS %.3f (kept)\n",
            sum(adhd_decisive), adhd_pool_rms_raw))
cat(sprintf("  negative sign share: %.3f\n", p_negative_binary))

# Binary variables are ordinal with num_categories = 1: one threshold each.
set.seed(20260731L)
binary_items <- sort(sample(ncol(adhd), num_variables))
binary_thresholds <- adhd_calibration$thresholds[binary_items, , drop = FALSE]
binary_endorsement <- adhd_calibration$endorsement[binary_items]
binary_targets <- cbind(1 - binary_endorsement, binary_endorsement)
cat(sprintf("  endorsement rates of the 15 sampled items: %.2f - %.2f\n",
            min(binary_endorsement), max(binary_endorsement)))

# --- Calibration: GGM arm (documented; no continuous data set ships) -----------

# Psychological GGMs carry partial correlations in roughly the 0.1-0.3 range.
# No continuous data set ships with bgms, so this arm is calibrated by
# assertion rather than to a fit: nonzero partial correlations are drawn
# uniformly on [0.1, 0.3] with the same negative-sign share as the ordinal
# arm, and the off-diagonal block is shrunk toward zero until the precision
# matrix is positive definite (the shrink factor is stored per network).
ggm_pcor_range <- c(0.10, 0.30)
p_negative_ggm <- p_negative_ordinal

# --- Build networks and data sets ----------------------------------------------

cat("=== networks and data sets ===\n")

network_seed <- function(arm_index, density_index, data_seed) {
  as.integer(70000L + arm_index * 1000L + density_index * 100L + data_seed)
}

dataset_rows <- NULL
network_rows <- NULL

build_ordinal_like <- function(arm, arm_index, weight_pool, p_negative,
                               thresholds, targets, num_categories,
                               data_seeds) {
  for (d in seq_along(densities)) {
    for (data_seed in data_seeds) {
      net_seed <- network_seed(arm_index, d, data_seed)
      network <- simulate_network(num_variables, densities[d], weight_pool,
                                  p_negative, seed = net_seed)
      cat(sprintf("  %s density %.1f seed %d: %d edges, calibrating...\n",
                  arm, densities[d], data_seed, sum(network$true_edge)))
      calibrated <- calibrate_thresholds(
        thresholds = thresholds, pairwise = network$pairwise,
        target_props = targets, num_categories = num_categories,
        seed = net_seed, verbose = FALSE
      )
      cat(sprintf("    final max|dprop| = %.4f\n", calibrated$deviation_final))
      network_rows <<- rbind(network_rows, data.frame(
        arm = arm, density = unname(densities[d]), data_seed = data_seed,
        network_seed = net_seed, num_edges = sum(network$true_edge),
        weight_rms = sqrt(mean(network$true_weight[network$true_edge]^2)),
        calibration_deviation = calibrated$deviation_final,
        shrink = NA_real_, stringsAsFactors = FALSE
      ))
      for (n in sample_sizes) {
        dataset_id <- sprintf("%s-d%02d-s%d-n%d", arm,
                              round(100 * densities[d]), data_seed, n)
        x <- bgms::simulate_mrf(
          num_states = n, num_variables = num_variables,
          num_categories = num_categories, pairwise = network$pairwise,
          main = calibrated$thresholds, variable_type = "ordinal",
          iter = 2000L, seed = net_seed + n
        )
        colnames(x) <- sprintf("V%02d", seq_len(num_variables))
        saveRDS(list(dataset_id = dataset_id, arm = arm,
                     density = unname(densities[d]), n = n,
                     data_seed = data_seed, x = x,
                     pairwise = network$pairwise,
                     thresholds = calibrated$thresholds,
                     true_weight = network$true_weight,
                     true_edge = network$true_edge),
                file.path(paths$datasets, paste0(dataset_id, ".rds")))
        dataset_rows <<- rbind(dataset_rows, data.frame(
          dataset_id = dataset_id, arm = arm,
          density = unname(densities[d]), n = n, data_seed = data_seed,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

build_ordinal_like("ordinal", 1L, weight_pool_ordinal, p_negative_ordinal,
                   wenchuan_thresholds, wenchuan_targets,
                   num_categories_ordinal, data_seeds_stage1)

build_ordinal_like("binary", 2L, weight_pool_binary, p_negative_binary,
                   binary_thresholds, binary_targets, 1L, data_seeds_stage2)

for (d in seq_along(densities)) {
  for (data_seed in data_seeds_stage2) {
    net_seed <- network_seed(3L, d, data_seed)
    network <- simulate_ggm_network(num_variables, densities[d],
                                    ggm_pcor_range, p_negative_ggm,
                                    seed = net_seed)
    cat(sprintf("  ggm density %.1f seed %d: %d edges, shrink %.3f\n",
                densities[d], data_seed, sum(network$true_edge),
                network$shrink))
    network_rows <- rbind(network_rows, data.frame(
      arm = "ggm", density = unname(densities[d]), data_seed = data_seed,
      network_seed = net_seed, num_edges = sum(network$true_edge),
      weight_rms = sqrt(mean(network$true_weight[network$true_edge]^2)),
      calibration_deviation = NA_real_, shrink = network$shrink,
      stringsAsFactors = FALSE
    ))
    for (n in sample_sizes) {
      dataset_id <- sprintf("ggm-d%02d-s%d-n%d", round(100 * densities[d]),
                            data_seed, n)
      x <- bgms::simulate_mrf(
        num_states = n, num_variables = num_variables, num_categories = 1L,
        pairwise = network$precision, main = rep(0, num_variables),
        variable_type = "continuous", seed = net_seed + n
      )
      colnames(x) <- sprintf("V%02d", seq_len(num_variables))
      saveRDS(list(dataset_id = dataset_id, arm = "ggm",
                   density = unname(densities[d]), n = n,
                   data_seed = data_seed, x = x,
                   pairwise = network$precision, thresholds = NULL,
                   true_weight = network$true_weight,
                   true_edge = network$true_edge, shrink = network$shrink),
              file.path(paths$datasets, paste0(dataset_id, ".rds")))
      dataset_rows <- rbind(dataset_rows, data.frame(
        dataset_id = dataset_id, arm = "ggm",
        density = unname(densities[d]), n = n, data_seed = data_seed,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# --- Manifest -------------------------------------------------------------------

stage_of_arm <- c(ordinal = 1L, binary = 2L, ggm = 2L)

manifest <- NULL
for (i in seq_len(nrow(dataset_rows))) {
  row <- dataset_rows[i, ]
  for (setting_name in names(test_settings)) {
    setting <- test_settings[[setting_name]]
    for (sampler_seed in sampler_seeds) {
      fit_id <- sprintf("%s-%s-seed%d", row$dataset_id, setting_name,
                        sampler_seed)
      manifest <- rbind(manifest, data.frame(
        fit_id = fit_id, stage = stage_of_arm[[row$arm]], arm = row$arm,
        kind = "test", dataset_id = row$dataset_id, density = row$density,
        n = row$n, data_seed = row$data_seed, setting = setting_name,
        iter = unname(setting["iter"]), warmup = unname(setting["warmup"]),
        chains = test_chains, sampler_seed = sampler_seed,
        stringsAsFactors = FALSE
      ))
    }
  }
  manifest <- rbind(manifest, data.frame(
    fit_id = sprintf("%s-reference", row$dataset_id),
    stage = stage_of_arm[[row$arm]], arm = row$arm, kind = "reference",
    dataset_id = row$dataset_id, density = row$density, n = row$n,
    data_seed = row$data_seed, setting = "reference",
    iter = unname(reference_setting["iter"]),
    warmup = unname(reference_setting["warmup"]),
    chains = reference_chains, sampler_seed = 99L,
    stringsAsFactors = FALSE
  ))
}

utils::write.csv(manifest, file.path(paths$root, "manifest.csv"),
                 row.names = FALSE)
utils::write.csv(network_rows, file.path(paths$root, "networks.csv"),
                 row.names = FALSE)

saveRDS(list(
  target_weight_rms = target_weight_rms,
  ordinal = list(weight_pool = weight_pool_ordinal,
                 weight_rms_raw = weight_pool_rms_raw,
                 p_negative = p_negative_ordinal,
                 items = ordinal_items,
                 thresholds = wenchuan_thresholds,
                 targets = wenchuan_targets),
  binary = list(weight_pool = weight_pool_binary,
                weight_rms_raw = adhd_pool_rms_raw,
                p_negative = p_negative_binary,
                items = binary_items,
                thresholds = binary_thresholds,
                targets = binary_targets),
  ggm = list(pcor_range = ggm_pcor_range, p_negative = p_negative_ggm)
), file.path(paths$root, "calibration.rds"))

cat("\n=== manifest ===\n")
print(table(manifest$arm, manifest$kind))
cat(sprintf("%d data sets, %d fits (%d test, %d reference)\n",
            nrow(dataset_rows), nrow(manifest),
            sum(manifest$kind == "test"), sum(manifest$kind == "reference")))
