# R/scripts/Z_prior_brute_force.R
#
# Brute-force MC estimate of Z_prior(Gamma; delta, epsilon, alpha, beta, sigma)
# for both discrete and continuous spike-and-slab variants, with a hot inner
# loop in C++ (R/src/Z_prior_brute.cpp) and a 12-core parallel driver
# (mclapply) over independent (graph, delta, epsilon, variant) configurations.
#
# Phase 1.4 of the spikeslab project; companion to notes/construction.md and
# notes/propriety.md. Used as ground truth for Laplace work and as a check on
# the closed-form result for the discrete empty graph.

suppressPackageStartupMessages({
  library(Rcpp)
  library(RcppArmadillo)
  library(parallel)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
CPP_SRC      <- file.path(PROJECT_ROOT, "R", "src", "Z_prior_brute.cpp")
OUT_BRUTE    <- file.path(PROJECT_ROOT, "R", "results", "Z_prior_brute_q234.rds")
OUT_CLOSED   <- file.path(PROJECT_ROOT, "R", "results", "Z_prior_closed_empty.rds")

# ----- compile the C++ inner loop (cached after first run) -----
Rcpp::sourceCpp(CPP_SRC)

# ----- build graphs (q in {2, 3, 4}) -----
build_graphs <- function() {
  g2_empty <- matrix(0L, 2, 2)
  g2_full  <- matrix(0L, 2, 2); g2_full[1,2] <- g2_full[2,1] <- 1L

  g3_empty <- matrix(0L, 3, 3)
  g3_path  <- matrix(0L, 3, 3)
  g3_path[1,2] <- g3_path[2,1] <- 1L
  g3_path[2,3] <- g3_path[3,2] <- 1L
  g3_tri   <- matrix(1L, 3, 3); diag(g3_tri) <- 0L

  g4_empty <- matrix(0L, 4, 4)
  g4_path  <- matrix(0L, 4, 4)
  g4_path[1,2] <- g4_path[2,1] <- 1L
  g4_path[2,3] <- g4_path[3,2] <- 1L
  g4_path[3,4] <- g4_path[4,3] <- 1L
  g4_cyc   <- g4_path; g4_cyc[1,4] <- g4_cyc[4,1] <- 1L
  g4_full  <- matrix(1L, 4, 4); diag(g4_full) <- 0L

  list(
    q2_empty = g2_empty, q2_full = g2_full,
    q3_empty = g3_empty, q3_path = g3_path, q3_tri = g3_tri,
    q4_empty = g4_empty, q4_path = g4_path, q4_cyc = g4_cyc, q4_full = g4_full
  )
}

# ----- build the per-edge slab-SD matrix -----
build_slab_sd <- function(G, sigma_on, sigma_off) {
  q <- nrow(G)
  M <- matrix(0, q, q)
  if (q > 1L) {
    for (i in 1:(q - 1L)) {
      for (j in (i + 1L):q) {
        M[i, j] <- if (G[i, j] == 1L) sigma_on else sigma_off
        M[j, i] <- M[i, j]
      }
    }
  }
  M
}

# ----- one job -----
run_one_config <- function(args) {
  G <- args$G
  storage.mode(G) <- "integer"
  res <- z_prior_brute_cpp(
    G, args$alpha, args$beta, args$slab_sd,
    args$delta, args$epsilon, args$n_draws, args$seed
  )
  data.frame(
    graph    = args$gname,
    q        = nrow(G),
    variant  = args$variant,
    alpha    = args$alpha,
    beta     = args$beta,
    sigma    = args$sigma_label,
    delta    = args$delta,
    epsilon  = args$epsilon,
    z_hat    = res$z_hat,
    z_se     = res$z_se,
    log_z    = res$log_z_hat,
    n_accept = res$n_accept,
    n_draws  = res$n_draws,
    seed     = args$seed,
    stringsAsFactors = FALSE
  )
}

# ----- driver: sweep over (graph, delta, epsilon, variant) -----
run_z_prior_sweep <- function(n_draws    = 1e6,
                              alpha      = 3.0,
                              beta       = 3.0,
                              sigma      = 1.0,
                              v0         = 0.1,
                              v1         = 1.0,
                              delta_grid = c(0, 0.5, 1.0),
                              eps_grid   = c(0, 0.05, 0.1),
                              seed_base  = 42L,
                              mc_cores   = 12L) {
  graphs <- build_graphs()
  jobs   <- list()
  k      <- 0L

  for (gname in names(graphs)) {
    G <- graphs[[gname]]
    slab_sd_d <- build_slab_sd(G, sigma_on = sigma, sigma_off = 0)
    slab_sd_c <- build_slab_sd(G, sigma_on = v1,    sigma_off = v0)

    for (delta in delta_grid) {
      for (epsilon in eps_grid) {
        k <- k + 1L
        jobs[[length(jobs) + 1L]] <- list(
          gname = gname, G = G, alpha = alpha, beta = beta,
          slab_sd = slab_sd_d, sigma_label = sigma,
          delta = delta, epsilon = epsilon,
          n_draws = n_draws, seed = seed_base + 10000L * k,
          variant = "discrete"
        )
        jobs[[length(jobs) + 1L]] <- list(
          gname = gname, G = G, alpha = alpha, beta = beta,
          slab_sd = slab_sd_c, sigma_label = v1,
          delta = delta, epsilon = epsilon,
          n_draws = n_draws, seed = seed_base + 10000L * k + 1L,
          variant = "continuous"
        )
      }
    }
  }

  cat(sprintf("Running %d jobs on %d cores (n_draws=%g each)...\n",
              length(jobs), mc_cores, n_draws))
  results <- parallel::mclapply(jobs, run_one_config, mc.cores = mc_cores)

  failed <- vapply(results, function(r) inherits(r, "try-error") || is.null(r),
                   logical(1L))
  if (any(failed)) {
    stop(sprintf("%d jobs failed", sum(failed)), call. = FALSE)
  }

  out <- do.call(rbind, results)
  out$rel_se <- out$z_se / pmax(out$z_hat, .Machine$double.eps)
  out
}

# ----- closed-form discrete empty-graph normaliser (for cross-check) -----
z_prior_closed_discrete_empty <- function(q, alpha, beta, delta, epsilon) {
  log_unit_integral <- lgamma(alpha + delta) +
    pgamma(beta * epsilon, shape = alpha + delta, lower.tail = FALSE, log.p = TRUE) -
    lgamma(alpha) - delta * log(beta)
  log_z <- q * log_unit_integral
  list(log_z = log_z, z = exp(log_z))
}

closed_form_table <- function(q_grid     = c(2, 3, 4),
                              alpha      = 3.0,
                              beta       = 3.0,
                              delta_grid = c(0, 0.5, 1.0),
                              eps_grid   = c(0, 0.05, 0.1)) {
  rows <- list()
  for (q in q_grid) {
    for (delta in delta_grid) {
      for (epsilon in eps_grid) {
        cf <- z_prior_closed_discrete_empty(q, alpha, beta, delta, epsilon)
        rows[[length(rows) + 1L]] <- data.frame(
          q = q, alpha = alpha, beta = beta,
          delta = delta, epsilon = epsilon,
          log_z_closed = cf$log_z, z_closed = cf$z,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

# ----- main driver (runs when sourced or Rscript-invoked) -----
if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  set.seed(42L)
  t0 <- Sys.time()

  brute_tab <- run_z_prior_sweep(
    n_draws    = 1e6,
    alpha      = 3.0,
    beta       = 3.0,
    sigma      = 1.0,
    v0         = 0.1,
    v1         = 1.0,
    delta_grid = c(0, 0.5, 1.0),
    eps_grid   = c(0, 0.05, 0.1),
    seed_base  = 42L,
    mc_cores   = 12L
  )
  t1 <- Sys.time()
  cat(sprintf("Brute-force sweep done in %.1f s. %d rows.\n",
              as.numeric(t1 - t0, units = "secs"), nrow(brute_tab)))
  saveRDS(brute_tab, OUT_BRUTE)
  cat(sprintf("Saved brute-force results to %s\n", OUT_BRUTE))

  closed_tab <- closed_form_table(
    q_grid     = c(2, 3, 4),
    alpha      = 3.0,
    beta       = 3.0,
    delta_grid = c(0, 0.5, 1.0),
    eps_grid   = c(0, 0.05, 0.1)
  )
  saveRDS(closed_tab, OUT_CLOSED)
  cat(sprintf("Saved closed-form table to %s\n", OUT_CLOSED))

  cat("\nSummary (discrete empty graphs, brute vs closed-form):\n")
  ed <- brute_tab[brute_tab$variant == "discrete" & grepl("_empty$", brute_tab$graph), ]
  ed$q <- as.integer(sub("q([0-9])_empty", "\\1", ed$graph))
  cmp <- merge(ed, closed_tab,
               by = c("q", "alpha", "beta", "delta", "epsilon"))
  cmp$rel_err <- (cmp$z_hat - cmp$z_closed) / cmp$z_closed
  print(cmp[order(cmp$q, cmp$delta, cmp$epsilon),
            c("q", "delta", "epsilon",
              "z_hat", "z_closed", "log_z", "log_z_closed", "rel_err", "z_se")],
        row.names = FALSE, digits = 5)

  cat("\nDone.\n")
}
