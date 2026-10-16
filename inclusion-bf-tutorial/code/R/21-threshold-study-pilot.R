# Known-truth calibration study, step 2: the TIMING PILOT.
# Design: the simulation appendix of the manuscript.
# (the projection sizes the full grid before it is launched).
#
# What it measures, and why that is enough to project the grid:
#
#   test fits      one fit per arm x setting x n (12 fits), run under the
#                  production configuration of 3 concurrent workers at 4
#                  chains each, so the measured seconds already contain the
#                  contention that 12 busy cores create. Two settings per
#                  (arm, n) cell identify a two-parameter cost model
#                  t = c_warmup * warmup + c_sample * iter, which then
#                  predicts any other setting in that cell.
#
#   reference runs one TRUNCATED run per arm (12 chains, 2,000 warmup +
#                  2,500 sampling, at the largest n), alone on the machine.
#                  Comparing it against what the same cost model predicts at
#                  4 chains gives the 12-chain penalty: 5 of this machine's
#                  15 cores are performance cores, so the last chains of a
#                  12-chain run land on efficiency cores and the run does not
#                  scale for free. The full reference setting is then
#                  projected with that penalty applied.
#
# Writes tables/threshold-study/pilot-timings.csv and
# tables/threshold-study/pilot-projection.csv, and prints the projection.
#
# Usage: Rscript R/21-threshold-study-pilot.R [--stage build|collect]
#   build   (default) writes the pilot manifests, to be run by R/22 workers
#   collect reads the checkpoints back and produces the projection

source("R/threshold-study-helpers.R")

parse_args <- function(args) {
  out <- list(stage = "build")
  i <- 1L
  while (i <= length(args)) {
    out[[sub("^--", "", args[i])]] <- args[i + 1L]
    i <- i + 2L
  }
  out
}
options <- parse_args(commandArgs(trailingOnly = TRUE))
paths <- study_paths()

manifest <- utils::read.csv(file.path(paths$root, "manifest.csv"),
                            stringsAsFactors = FALSE)

pilot_density <- 0.4          # the middle density, one network per arm
pilot_data_seed <- 1L
pilot_sampler_seed <- 1L
truncated_reference <- c(warmup = 2000L, iter = 2500L)

test_pilot_path <- file.path(paths$root, "manifest-pilot-test.csv")
reference_pilot_path <- file.path(paths$root, "manifest-pilot-reference.csv")

# --- build ---------------------------------------------------------------------

if (options$stage == "build") {
  test_pilot <- manifest[manifest$kind == "test" &
                           manifest$density == pilot_density &
                           manifest$data_seed == pilot_data_seed &
                           manifest$sampler_seed == pilot_sampler_seed, ]
  # Order so that the three concurrent workers each get one arm's worth of
  # comparable work rather than all the cheap fits landing on one slot.
  test_pilot <- test_pilot[order(test_pilot$n, test_pilot$setting,
                                 test_pilot$arm), ]

  reference_pilot <- manifest[manifest$kind == "reference" &
                                manifest$density == pilot_density &
                                manifest$data_seed == pilot_data_seed &
                                manifest$n == max(manifest$n), ]
  reference_pilot$fit_id <- paste0(reference_pilot$fit_id, "-truncated")
  reference_pilot$warmup <- unname(truncated_reference["warmup"])
  reference_pilot$iter <- unname(truncated_reference["iter"])
  reference_pilot$setting <- "reference-truncated"

  utils::write.csv(test_pilot, test_pilot_path, row.names = FALSE)
  utils::write.csv(reference_pilot, reference_pilot_path, row.names = FALSE)
  cat(sprintf("pilot manifests written: %d test fits, %d truncated references\n",
              nrow(test_pilot), nrow(reference_pilot)))
  quit(save = "no")
}

# --- collect and project --------------------------------------------------------

read_timings <- function(path) {
  pilot <- utils::read.csv(path, stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(pilot)), function(i) {
    checkpoint_path <- fit_checkpoint_path(pilot$fit_id[i])
    if (!file.exists(checkpoint_path)) return(NULL)
    checkpoint <- readRDS(checkpoint_path)
    cbind(pilot[i, c("fit_id", "arm", "kind", "n", "density", "setting",
                     "iter", "warmup", "chains")],
          data.frame(seconds_fit = checkpoint$seconds_fit,
                     seconds_total = checkpoint$seconds_total))
  })
  do.call(rbind, rows)
}

density_probe_path <- file.path(paths$root, "manifest-density-probe.csv")

timings <- rbind(read_timings(test_pilot_path),
                 read_timings(reference_pilot_path))
if (is.null(timings)) stop("no pilot checkpoints found; run the workers first",
                           call. = FALSE)
timings$seconds_observables <- timings$seconds_total - timings$seconds_fit
utils::write.csv(timings, file.path(paths$root, "pilot-timings.csv"),
                 row.names = FALSE)

cat("=== pilot timings (seconds) ===\n")
print(timings[, c("arm", "n", "setting", "chains", "seconds_fit",
                  "seconds_observables")], row.names = FALSE, digits = 4)

# Cost model per (arm, n) cell from the two test settings:
#   seconds = c_warmup * warmup + c_sample * iter
test_timings <- timings[timings$kind == "test", ]
cells <- unique(test_timings[, c("arm", "n")])
cost_model <- NULL
for (k in seq_len(nrow(cells))) {
  cell <- test_timings[test_timings$arm == cells$arm[k] &
                         test_timings$n == cells$n[k], ]
  if (nrow(cell) < 2L) next
  coefficients <- qr.solve(as.matrix(cell[, c("warmup", "iter")]),
                           cell$seconds_fit)
  cost_model <- rbind(cost_model, data.frame(
    arm = cells$arm[k], n = cells$n[k],
    c_warmup = coefficients[1], c_sample = coefficients[2],
    seconds_observables = mean(cell$seconds_observables),
    stringsAsFactors = FALSE
  ))
}
cat("\n=== cost model (seconds per iteration, 4 chains, 3 concurrent) ===\n")
print(cost_model, row.names = FALSE, digits = 3)

# Density multipliers. The pilot times one density only, and run time need
# not be flat in density (the long setting in particular slows down at the
# top density), so a projection built on the middle density alone can
# understate the dense cells. The multipliers come from an optional
# density-probe manifest (ordinal arm, n = 1000, all three densities) and
# default to 1 when it is absent. They are estimated relative to the pilot
# density and applied to every arm: only the ordinal arm is probed, so the
# other two borrow its factors.
density_multiplier <- data.frame(setting = character(), density = numeric(),
                                 multiplier = numeric(),
                                 stringsAsFactors = FALSE)
if (file.exists(density_probe_path)) {
  probe_columns <- c("arm", "n", "density", "setting", "seconds_fit")
  probe <- rbind(read_timings(density_probe_path)[, probe_columns],
                 timings[timings$kind == "test" & timings$n == 1000 &
                           timings$arm == "ordinal", probe_columns])
  probe <- probe[probe$arm == "ordinal" & probe$n == 1000, ]
  for (setting_name in unique(probe$setting)) {
    block <- probe[probe$setting == setting_name, ]
    base <- block$seconds_fit[block$density == pilot_density]
    if (!length(base)) next
    density_multiplier <- rbind(density_multiplier, data.frame(
      setting = setting_name, density = block$density,
      multiplier = block$seconds_fit / base[1], stringsAsFactors = FALSE
    ))
  }
  cat("\n=== density multipliers (ordinal, n = 1000, relative to density",
      pilot_density, ") ===\n")
  print(stats::xtabs(multiplier ~ setting + density, density_multiplier),
        digits = 3)
}

lookup_multiplier <- function(setting, density) {
  # Reference runs are long runs in character (10,000 warmup, 25,000
  # sampling), so they take the long setting's density factor. That is the
  # conservative choice of the two available.
  key <- if (setting == "default") "default" else "long"
  hit <- density_multiplier$multiplier[density_multiplier$setting == key &
                                         density_multiplier$density == density]
  if (length(hit)) hit[1] else 1
}

predict_seconds <- function(arm, n, warmup, iter, density = pilot_density,
                            setting = "long") {
  row <- cost_model[cost_model$arm == arm & cost_model$n == n, ]
  (row$c_warmup * warmup + row$c_sample * iter) *
    lookup_multiplier(setting, density)
}

# 12-chain penalty: observed truncated reference against the 4-chain
# prediction for the same iteration counts.
reference_timings <- timings[timings$kind == "reference", ]
reference_timings$predicted_4chain <- vapply(seq_len(nrow(reference_timings)),
  function(i) predict_seconds(reference_timings$arm[i], reference_timings$n[i],
                              reference_timings$warmup[i],
                              reference_timings$iter[i],
                              density = pilot_density, setting = "long"),
  numeric(1))
reference_timings$chain_penalty <- reference_timings$seconds_fit /
  reference_timings$predicted_4chain
cat("\n=== 12-chain penalty (truncated reference vs 4-chain prediction) ===\n")
print(reference_timings[, c("arm", "n", "warmup", "iter", "seconds_fit",
                            "predicted_4chain", "chain_penalty")],
      row.names = FALSE, digits = 3)

penalty <- stats::setNames(reference_timings$chain_penalty,
                           reference_timings$arm)
# Extraction cost per stored draw at 12 chains, so the reference runs' own
# observables pass is projected too rather than assumed free.
reference_extraction_per_draw <- stats::setNames(
  reference_timings$seconds_observables /
    (reference_timings$iter * reference_timings$chains),
  reference_timings$arm
)

# Reference runs are cheaper at n = 200 than at n = 1000 in the same
# proportion the test fits are, so the model is applied per (arm, n) cell and
# the penalty per arm.
full_manifest <- manifest
full_manifest$predicted_seconds <- vapply(seq_len(nrow(full_manifest)),
  function(i) {
    row <- full_manifest[i, ]
    base <- predict_seconds(row$arm, row$n, row$warmup, row$iter,
                            density = row$density, setting = row$setting)
    if (row$kind == "reference") {
      base * penalty[[row$arm]] +
        reference_extraction_per_draw[[row$arm]] * row$iter * row$chains
    } else {
      base + cost_model$seconds_observables[cost_model$arm == row$arm &
                                              cost_model$n == row$n]
    }
  }, numeric(1))
if (any(full_manifest$predicted_seconds <= 0)) {
  warning("the cost model predicts non-positive times; inspect pilot-timings.csv",
          call. = FALSE)
}

# Wall clock: test fits run 3 at a time (the concurrency the pilot was timed
# under), reference runs one at a time.
projection <- NULL
for (stage in sort(unique(full_manifest$stage))) {
  for (arm in unique(full_manifest$arm[full_manifest$stage == stage])) {
    for (kind in c("test", "reference")) {
      block <- full_manifest[full_manifest$stage == stage &
                               full_manifest$arm == arm &
                               full_manifest$kind == kind, ]
      if (!nrow(block)) next
      concurrency <- if (kind == "test") 3 else 1
      projection <- rbind(projection, data.frame(
        stage = stage, arm = arm, kind = kind, fits = nrow(block),
        cpu_hours = sum(block$predicted_seconds) / 3600,
        wall_hours = sum(block$predicted_seconds) / concurrency / 3600,
        stringsAsFactors = FALSE
      ))
    }
  }
}
utils::write.csv(projection, file.path(paths$root, "pilot-projection.csv"),
                 row.names = FALSE)
# The projection sizes the grid before it runs, so it is kept with the
# results in tables/ as well as in the checkpoint tree.
utils::write.csv(projection, "tables/threshold-study-pilot-projection.csv",
                 row.names = FALSE)
utils::write.csv(timings, "tables/threshold-study-pilot-timings.csv",
                 row.names = FALSE)

cat("\n=== PROJECTION ===\n")
print(projection, row.names = FALSE, digits = 3)
cat(sprintf("\nstage 1 wall clock: %.1f h   stage 2 wall clock: %.1f h",
            sum(projection$wall_hours[projection$stage == 1]),
            sum(projection$wall_hours[projection$stage == 2])))
cat(sprintf("\nTOTAL wall clock: %.1f h over %d fits\n",
            sum(projection$wall_hours), sum(projection$fits)))
cat(sprintf("edge-fits: stage 1 %d, stage 2 %d\n",
            105L * sum(projection$fits[projection$stage == 1 &
                                         projection$kind == "test"]),
            105L * sum(projection$fits[projection$stage == 2 &
                                         projection$kind == "test"])))

# --- does the generated data span the diagnostic range? -------------------------

# The sizing argument of the design rests on every fit contributing edge-fits
# across the whole range (settled, few-flip, near-boundary, healthy). If the
# pilot's edges are all decisive and settled, the grid would buy replication
# of one regime rather than coverage, and that is worth knowing before the
# compute is spent rather than after.
pilot_test <- utils::read.csv(test_pilot_path, stringsAsFactors = FALSE)
observables <- do.call(rbind, lapply(pilot_test$fit_id, function(id) {
  path <- fit_checkpoint_path(id)
  if (!file.exists(path)) return(NULL)
  readRDS(path)$observables
}))

flip_regime <- cut(observables$total_flips,
                   breaks = c(-1, 0, 50, 500, 2000, Inf),
                   labels = c("0 (settled)", "1-50", "51-500", "501-2000",
                              ">2000"))
coverage <- as.data.frame(table(setting = observables$setting,
                                flips = flip_regime))
cat("\n=== flip-count regime coverage in the pilot fits ===\n")
print(stats::xtabs(Freq ~ setting + flips, coverage))

cat("\n=== verdict x truth in the pilot fits (test settings, not truth-graded) ===\n")
print(table(verdict = observables$verdict, true_edge = observables$true_edge))

cat("\n=== log10 BF spread by arm ===\n")
print(do.call(rbind, lapply(split(observables, observables$arm), function(d) {
  data.frame(edge_fits = nrow(d),
             q10 = stats::quantile(d$log10_bf, 0.10, na.rm = TRUE),
             median = stats::median(d$log10_bf, na.rm = TRUE),
             q90 = stats::quantile(d$log10_bf, 0.90, na.rm = TRUE),
             near_boundary = mean(abs(d$log10_bf) < 1.5, na.rm = TRUE),
             settled = mean(d$total_flips == 0))
})), digits = 3)

utils::write.csv(coverage, file.path(paths$root, "pilot-coverage.csv"),
                 row.names = FALSE)
