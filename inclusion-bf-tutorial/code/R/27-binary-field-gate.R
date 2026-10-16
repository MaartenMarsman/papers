# Binary arm rerun, step 1: the PER-NODE FIELD GATE.
# Design: the simulation appendix of the manuscript (follow-up stage).
#
# The first binary arm calibrated the per-EDGE weight distribution to ADHD but
# imposed densities of 20/40/60% on top of it without checking the per-NODE
# field. Couplings are therefore rescaled by
#
#     s = (calibrating mean degree) / (imposed mean degree)
#
# so that the expected total coupling reaching a node matches the calibrating
# fit's. That is the prior; the GATE is what decides. Each candidate scale is
# simulated and its constant-row share and mean pairwise phi are compared with
# the same statistics measured on the 15 sampled ADHD items themselves.
#
# WHY A GATE RATHER THAN A RULE. Holding the per-edge weights fixed while the
# imposed density varies moves the per-node field in both directions: a
# sparse network can be too weakly coupled relative to the calibrating data
# (few constant rows, low phi) and a dense one too strongly coupled (most rows
# constant), so no single rescaling direction is right for every density. A
# high constant-row share in a symptom screener is also a real population
# feature (many respondents endorse nothing, some endorse most), not
# automatically a phase-transition artefact; comparing each cell with the
# statistics of the calibrating data is what tells the two apart.
#
# Output: tables/threshold-study-binary-field-gate.csv

source("R/threshold-study-helpers.R")

paths <- study_paths()
num_variables <- 15L
densities <- c(0.2, 0.4, 0.6)
gate_n <- 1000L

calibration <- readRDS(file.path(paths$root, "calibration.rds"))
adhd_calibration <- readRDS(file.path(paths$root, "fit-adhd-calibration.rds"))

data("ADHD", package = "bgms")
adhd <- ADHD[, setdiff(colnames(ADHD), grep("group", colnames(ADHD),
                                            ignore.case = TRUE, value = TRUE))]
adhd <- adhd[stats::complete.cases(adhd), ]
sampled <- as.matrix(adhd[, calibration$binary$items])

field_statistics <- function(x) {
  row_sums <- rowSums(x)
  correlations <- stats::cor(x)
  c(constant_row_share = mean(row_sums == 0 | row_sums == ncol(x)),
    mean_abs_phi = mean(abs(correlations[upper.tri(correlations)])),
    row_sum_sd = stats::sd(row_sums))
}

target <- field_statistics(sampled)
cat("=== gate targets, 15 sampled ADHD items ===\n")
print(round(target, 4))

decisive <- adhd_calibration$pip > 0.9
calibrating_degree <- 2 * sum(decisive) / ncol(adhd)
cat(sprintf("\ncalibrating mean degree: %.2f (%d decisive edges over %d items)\n",
            calibrating_degree, sum(decisive), ncol(adhd)))

# Candidate scales: the degree rule, plus a bracket around it so the gate can
# choose rather than merely confirm.
scale_multipliers <- c(0.7, 1.0, 1.4)

gate <- NULL
for (density in densities) {
  imposed_degree <- 2 * round(density * choose(num_variables, 2)) / num_variables
  degree_scale <- calibrating_degree / imposed_degree
  for (multiplier in scale_multipliers) {
    scale <- degree_scale * multiplier
    network <- simulate_network(
      num_variables, density,
      weight_pool = calibration$binary$weight_pool * scale,
      p_negative = calibration$binary$p_negative,
      seed = 72000L + round(100 * density) * 10L + round(100 * multiplier)
    )
    calibrated <- calibrate_thresholds(
      thresholds = calibration$binary$thresholds,
      pairwise = network$pairwise,
      target_props = calibration$binary$targets,
      num_categories = 1L, seed = 7L, verbose = FALSE
    )
    x <- bgms::simulate_mrf(gate_n, num_variables, 1L, network$pairwise,
                            calibrated$thresholds, "ordinal", iter = 2000L,
                            seed = 99L)
    statistics <- field_statistics(x)
    gate <- rbind(gate, data.frame(
      density = density, imposed_degree = imposed_degree,
      degree_scale = degree_scale, multiplier = multiplier, scale = scale,
      constant_row_share = unname(statistics["constant_row_share"]),
      mean_abs_phi = unname(statistics["mean_abs_phi"]),
      row_sum_sd = unname(statistics["row_sum_sd"]),
      phi_gap = unname(abs(statistics["mean_abs_phi"] - target["mean_abs_phi"])),
      constant_gap = unname(abs(statistics["constant_row_share"] -
                                  target["constant_row_share"])),
      stringsAsFactors = FALSE
    ))
    cat(sprintf("density %.1f scale %.3f (%.1fx degree rule): const %.3f phi %.3f\n",
                density, scale, multiplier,
                statistics["constant_row_share"], statistics["mean_abs_phi"]))
  }
}

# Choose per density on the summed relative gap of the two gate statistics.
gate$gate_score <- gate$phi_gap / target["mean_abs_phi"] +
  gate$constant_gap / target["constant_row_share"]
chosen <- do.call(rbind, lapply(split(gate, gate$density), function(block) {
  block[which.min(block$gate_score), ]
}))

utils::write.csv(gate, "tables/threshold-study-binary-field-gate.csv",
                 row.names = FALSE)

cat("\n=== gate table ===\n")
print(gate[, c("density", "scale", "multiplier", "constant_row_share",
               "mean_abs_phi", "gate_score")], row.names = FALSE, digits = 3)
cat(sprintf("\ntargets: constant rows %.3f, mean |phi| %.3f\n",
            target["constant_row_share"], target["mean_abs_phi"]))
cat("\n=== chosen scale per density ===\n")
print(chosen[, c("density", "scale", "constant_row_share", "mean_abs_phi",
                 "gate_score")], row.names = FALSE, digits = 3)

cat("\n=== which cells must be rerun? ===\n")
cat("The first arm ran at scale 1.0 throughout. A cell is kept only if its\n")
cat("chosen scale is within 10% of 1.0.\n")
chosen$keep <- abs(chosen$scale - 1) <= 0.10
print(chosen[, c("density", "scale", "keep")], row.names = FALSE, digits = 3)
saveRDS(chosen, file.path(paths$root, "binary-field-gate-chosen.rds"))
