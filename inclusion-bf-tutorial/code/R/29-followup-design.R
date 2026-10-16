# Follow-up design: gated binary rerun cells and a near-boundary GGM arm.
# Design: the simulation appendix of the manuscript (follow-up stage).
#
# BINARY. The gate of R/27-binary-field-gate.R chose a coupling scale per
# density against the ADHD field statistics. Only the cells whose chosen scale
# differs from the 1.0 the first arm ran at are regenerated; their superseded
# datasets and checkpoints
# are archived rather than deleted so the change stays auditable. The kept
# cells' fit_ids and dataset_ids are untouched, so nothing downstream has to
# know which cells were rerun.
#
# GGM. A new arm designed for ERROR YIELD rather than grid symmetry: verdict
# errors live within 0.25 of a threshold, so fits are spent on settings that
# put many edges THERE instead of on settled edges that can never produce one.
# The existing GGM reference runs give the map from partial correlation to
# Bayes factor, so each sample size gets the band that straddles a verdict
# boundary. The bands and densities below are design constants read off the
# stage-2 GGM reference results:
#
#   n = 200   true edges cross log10 BF = 1 near |pcor| 0.22, and true
#             NON-edges already sit near the absence boundary, so both edge
#             classes are near-boundary at once
#   n = 500   interpolated band
#   n = 1000  true edges cross near |pcor| 0.09
#
# The band is specified AFTER the positive-definiteness shrinkage, which
# multiplies every partial correlation down, so the requested band is divided
# by the realized shrink factor and the network rebuilt until the realized
# band lands where it was asked to.
#
# Output: regenerated/added datasets, manifest rows, and
# tables/threshold-study-followup-design.csv

source("R/threshold-study-helpers.R")

paths <- study_paths()
num_variables <- 15L
sampler_seeds_binary <- 1:5
sampler_seeds_ggm <- 1:10      # more seeds: edge-fits are what yield errors
test_settings <- list(default = c(warmup = 2000L, iter = 2000L),
                      long = c(warmup = 5000L, iter = 10000L))
reference_setting <- c(warmup = 10000L, iter = 25000L)
reference_chains <- 12L
test_chains <- 4L

archive <- file.path(paths$root, "superseded")
if (!dir.exists(archive)) dir.create(archive, recursive = TRUE)

calibration <- readRDS(file.path(paths$root, "calibration.rds"))
chosen <- readRDS(file.path(paths$root, "binary-field-gate-chosen.rds"))
manifest <- utils::read.csv(file.path(paths$root, "manifest.csv"),
                            stringsAsFactors = FALSE)

design_log <- NULL

# --- Binary: regenerate only the cells the gate rejected -------------------------

cat("=== binary rerun cells ===\n")
for (i in seq_len(nrow(chosen))) {
  density <- chosen$density[i]
  scale <- chosen$scale[i]
  if (isTRUE(chosen$keep[i])) {
    cat(sprintf("density %.1f: scale %.3f is within 10%% of 1 -- KEPT, checkpoints untouched\n",
                density, scale))
    design_log <- rbind(design_log, data.frame(
      arm = "binary", density = density, scale = scale, action = "kept",
      stringsAsFactors = FALSE))
    next
  }
  existing <- file.path(paths$datasets,
                        sprintf("binary-d%02d-s1-n200.rds", round(100 * density)))
  if (file.exists(existing) &&
      isTRUE(all.equal(readRDS(existing)$coupling_scale, scale))) {
    cat(sprintf("density %.1f: already regenerated at scale %.3f -- skipping\n",
                density, scale))
    design_log <- rbind(design_log, data.frame(
      arm = "binary", density = density, scale = scale,
      action = "regenerated", stringsAsFactors = FALSE))
    next
  }
  cat(sprintf("density %.1f: scale %.3f -- regenerating\n", density, scale))

  network_seed <- 72000L + round(100 * density)
  network <- simulate_network(
    num_variables, density,
    weight_pool = calibration$binary$weight_pool * scale,
    p_negative = calibration$binary$p_negative, seed = network_seed
  )
  calibrated <- calibrate_thresholds(
    thresholds = calibration$binary$thresholds, pairwise = network$pairwise,
    target_props = calibration$binary$targets, num_categories = 1L,
    seed = network_seed, verbose = FALSE
  )

  for (n in c(200L, 1000L)) {
    dataset_id <- sprintf("binary-d%02d-s1-n%d", round(100 * density), n)
    old_path <- file.path(paths$datasets, paste0(dataset_id, ".rds"))
    if (file.exists(old_path)) {
      file.rename(old_path, file.path(archive, paste0(dataset_id, ".rds")))
    }
    for (fit_id in manifest$fit_id[manifest$dataset_id == dataset_id]) {
      old_fit <- fit_checkpoint_path(fit_id)
      if (file.exists(old_fit)) {
        file.rename(old_fit, file.path(archive, paste0(fit_id, ".rds")))
      }
    }
    x <- bgms::simulate_mrf(n, num_variables, 1L, network$pairwise,
                            calibrated$thresholds, "ordinal", iter = 2000L,
                            seed = network_seed + n)
    colnames(x) <- sprintf("V%02d", seq_len(num_variables))
    saveRDS(list(dataset_id = dataset_id, arm = "binary", density = density,
                 n = n, data_seed = 1L, x = x, pairwise = network$pairwise,
                 thresholds = calibrated$thresholds,
                 true_weight = network$true_weight,
                 true_edge = network$true_edge, coupling_scale = scale),
            file.path(paths$datasets, paste0(dataset_id, ".rds")))
    row_sums <- rowSums(x)
    correlations <- stats::cor(x)
    cat(sprintf("  %s: constant rows %.3f, mean |phi| %.3f\n", dataset_id,
                mean(row_sums == 0 | row_sums == num_variables),
                mean(abs(correlations[upper.tri(correlations)]))))
    design_log <- rbind(design_log, data.frame(
      arm = "binary", density = density, scale = scale, action = "regenerated",
      stringsAsFactors = FALSE))
  }
}

# --- GGM: a near-boundary arm ----------------------------------------------------

cat("\n=== GGM near-boundary arm ===\n")
# Density is chosen PER SAMPLE SIZE, not held fixed. The n = 200 band needs
# partial correlations near 0.26, and at density 0.4 a 15-variable precision
# matrix with entries that large is not positive definite: the shrinkage
# would bite hard and differ between seeds, making the two seeds
# incomparable. Density 0.3 admits the whole band with no shrinkage at all.
ggm_bands <- list(`200` = list(band = c(0.17, 0.26), density = 0.30),
                  `500` = list(band = c(0.11, 0.17), density = 0.40),
                  `1000` = list(band = c(0.075, 0.115), density = 0.40))
ggm_rows <- NULL

for (n_label in names(ggm_bands)) {
  n <- as.integer(n_label)
  band <- ggm_bands[[n_label]]$band
  ggm_density <- ggm_bands[[n_label]]$density
  for (data_seed in 1:2) {
    seed <- 74000L + n + data_seed
    # Compensate for the shrinkage: ask for a band scaled up by 1 / shrink
    # until the realized post-shrink band lands where it was requested.
    requested <- band
    for (attempt in 1:6) {
      network <- simulate_ggm_network(num_variables, ggm_density, requested,
                                      calibration$ggm$p_negative, seed = seed)
      realized <- range(abs(network$true_weight[network$true_edge]))
      if (abs(realized[2] - band[2]) < 0.005) break
      requested <- requested * (band[2] / realized[2])
    }
    cat(sprintf("  n=%4d seed %d: |pcor| in [%.3f, %.3f], shrink %.3f, %d edges\n",
                n, data_seed, realized[1], realized[2], network$shrink,
                sum(network$true_edge)))
    dataset_id <- sprintf("ggmnb-n%d-s%d", n, data_seed)
    x <- bgms::simulate_mrf(n, num_variables, num_categories = 1L,
                            pairwise = network$precision,
                            main = rep(0, num_variables),
                            variable_type = "continuous", seed = seed + 7L)
    colnames(x) <- sprintf("V%02d", seq_len(num_variables))
    saveRDS(list(dataset_id = dataset_id, arm = "ggmnb", density = ggm_density,
                 n = n, data_seed = data_seed, x = x,
                 pairwise = network$precision, thresholds = NULL,
                 true_weight = network$true_weight,
                 true_edge = network$true_edge, shrink = network$shrink),
            file.path(paths$datasets, paste0(dataset_id, ".rds")))
    ggm_rows <- rbind(ggm_rows, data.frame(dataset_id = dataset_id, n = n,
                                           data_seed = data_seed,
                                           density = ggm_density,
                                           stringsAsFactors = FALSE))
    design_log <- rbind(design_log, data.frame(
      arm = "ggmnb", density = ggm_density, scale = NA_real_,
      action = "new", stringsAsFactors = FALSE))
  }
}

# --- Manifest rows for the new GGM arm (stage 3) --------------------------------

new_rows <- NULL
for (i in seq_len(nrow(ggm_rows))) {
  row <- ggm_rows[i, ]
  for (setting_name in names(test_settings)) {
    setting <- test_settings[[setting_name]]
    for (sampler_seed in sampler_seeds_ggm) {
      new_rows <- rbind(new_rows, data.frame(
        fit_id = sprintf("%s-%s-seed%d", row$dataset_id, setting_name,
                         sampler_seed),
        stage = 3L, arm = "ggmnb", kind = "test",
        dataset_id = row$dataset_id, density = row$density, n = row$n,
        data_seed = row$data_seed, setting = setting_name,
        iter = unname(setting["iter"]), warmup = unname(setting["warmup"]),
        chains = test_chains, sampler_seed = sampler_seed,
        stringsAsFactors = FALSE))
    }
  }
  new_rows <- rbind(new_rows, data.frame(
    fit_id = sprintf("%s-reference", row$dataset_id), stage = 3L,
    arm = "ggmnb", kind = "reference", dataset_id = row$dataset_id,
    density = row$density, n = row$n, data_seed = row$data_seed,
    setting = "reference", iter = unname(reference_setting["iter"]),
    warmup = unname(reference_setting["warmup"]), chains = reference_chains,
    sampler_seed = 99L, stringsAsFactors = FALSE))
}

manifest <- manifest[manifest$arm != "ggmnb", ]
manifest <- rbind(manifest, new_rows)
utils::write.csv(manifest, file.path(paths$root, "manifest.csv"),
                 row.names = FALSE)

rerun_binary <- manifest[manifest$arm == "binary" &
                           !vapply(manifest$fit_id, fit_is_done, logical(1)), ]
utils::write.csv(rerun_binary,
                 file.path(paths$root, "manifest-binary-rerun.csv"),
                 row.names = FALSE)
utils::write.csv(new_rows[new_rows$kind == "test", ],
                 file.path(paths$root, "manifest-ggmnb-test.csv"),
                 row.names = FALSE)
utils::write.csv(new_rows[new_rows$kind == "reference", ],
                 file.path(paths$root, "manifest-ggmnb-reference.csv"),
                 row.names = FALSE)
utils::write.csv(design_log, "tables/threshold-study-followup-design.csv",
                 row.names = FALSE)

cat(sprintf("\nbinary fits to rerun: %d (%d test, %d reference)\n",
            nrow(rerun_binary), sum(rerun_binary$kind == "test"),
            sum(rerun_binary$kind == "reference")))
cat(sprintf("GGM near-boundary fits to run: %d (%d test, %d reference)\n",
            nrow(new_rows), sum(new_rows$kind == "test"),
            sum(new_rows$kind == "reference")))
