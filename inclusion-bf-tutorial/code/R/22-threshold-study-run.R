# Known-truth calibration study, step 3: the fit worker.
# Design: the simulation appendix of the manuscript.
#
# One process runs a slice of a manifest and checkpoints the per-edge
# observables of every fit to tables/threshold-study/fits/<fit_id>.rds.
# Fits whose checkpoint already exists are skipped, so the grid is resumable
# and partial results are analyzable at any point.
#
# Concurrency is by process, not by thread: bgm() spreads its own chains over
# `--cores`, so the test arm runs as several of these workers at 4 cores each
# and the reference arm as a single worker at 12.
#
# Usage:
#   Rscript R/22-threshold-study-run.R --manifest <csv> [--slot k --nslots K]
#                                      [--cores 4] [--nboot 99] [--label tag]

source("R/threshold-study-helpers.R")

parse_args <- function(args) {
  out <- list(manifest = file.path(study_dir, "manifest.csv"), slot = 0L,
              nslots = 1L, cores = NULL, nboot = 99L, label = "run")
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    value <- args[i + 1L]
    out[[key]] <- switch(key,
                         manifest = value,
                         label = value,
                         as.integer(value))
    i <- i + 2L
  }
  out
}

options <- parse_args(commandArgs(trailingOnly = TRUE))
paths <- study_paths()

manifest <- utils::read.csv(options$manifest, stringsAsFactors = FALSE)
mine <- manifest[(seq_len(nrow(manifest)) - 1L) %% options$nslots ==
                   options$slot, , drop = FALSE]

log_path <- file.path(paths$logs, sprintf("%s-slot%d.log", options$label,
                                          options$slot))
log_line <- function(...) {
  message <- sprintf(...)
  cat(message, file = log_path, append = TRUE)
  cat(message)
  utils::flush.console()
}

log_line("worker %s slot %d/%d: %d fits, %d cores, nboot %d\n",
         options$label, options$slot, options$nslots, nrow(mine),
         if (is.null(options$cores)) 0L else options$cores, options$nboot)

dataset_cache <- new.env(parent = emptyenv())
load_dataset <- function(dataset_id) {
  if (is.null(dataset_cache[[dataset_id]])) {
    dataset_cache[[dataset_id]] <- readRDS(
      file.path(paths$datasets, paste0(dataset_id, ".rds"))
    )
  }
  dataset_cache[[dataset_id]]
}

for (i in seq_len(nrow(mine))) {
  row <- mine[i, ]
  if (fit_is_done(row$fit_id)) {
    log_line("[%s] skip %s (checkpoint exists)\n",
             format(Sys.time(), "%H:%M:%S"), row$fit_id)
    next
  }
  dataset <- load_dataset(row$dataset_id)
  n_boot <- if (row$kind == "reference") 0L else options$nboot
  elapsed <- tryCatch(
    run_and_checkpoint(row, dataset, n_boot = n_boot, cores = options$cores),
    error = function(e) {
      log_line("[%s] FAILED %s: %s\n", format(Sys.time(), "%H:%M:%S"),
               row$fit_id, conditionMessage(e))
      NA_real_
    }
  )
  if (!is.na(elapsed)) {
    log_line("[%s] done %s in %.1f s (%d/%d)\n", format(Sys.time(), "%H:%M:%S"),
             row$fit_id, elapsed, i, nrow(mine))
  }
}

log_line("worker %s slot %d finished\n", options$label, options$slot)
