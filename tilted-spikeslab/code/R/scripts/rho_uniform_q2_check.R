# R/scripts/rho_uniform_q2_check.R
# Existence check at q = 2 (single edge, no rest block): the analytic claim is
#   p(rho) ∝ (1 - rho^2)^delta  as sigma -> infinity (flat slab),
# hence beta_hat -> delta, and delta = 0 gives EXACTLY uniform rho.
# Output: printed table only.

suppressPackageStartupMessages({ library(Rcpp) })
PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
sourceCpp(file.path(PROJECT_ROOT, "R", "src", "Z_prior_gwishart_slab_rho.cpp"))
set.seed(20260724L)

fit_beta <- function(rho) {
  rho <- rho[is.finite(rho) & abs(rho) < 1]
  S <- sum(log1p(-rho^2)); n <- length(rho)
  negll <- function(b) -(b * S - n * (0.5*log(pi) + lgamma(b+1) - lgamma(b+1.5)))
  optimize(negll, c(-0.95, 40))$minimum
}

grid <- expand.grid(delta = c(0, 0.5, 1), sigma = c(1, 4, 20, 100))
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  d <- grid$delta[i]; s <- grid$sigma[i]
  o <- gwishart_slab_record_rho_cpp(q = 2L, delta = d, beta = 1, sigma = s,
        gamma_prob = 1, n_target = 20000L, n_proposal_cap = 2000000L,
        max_inner_tries = 1000L, seed = 1000L + i, condition_on_edge_01 = TRUE)
  rho <- -o$k12 / sqrt(o$k11 * o$k22); rho <- rho[is.finite(rho)]
  data.frame(delta = d, sigma = s, n = length(rho),
             rho_sd = sd(rho), beta_hat = fit_beta(rho))
}))
cat("Uniform <=> beta = 0 (SD = 0.577). Prediction: beta_hat -> delta as sigma -> Inf.\n\n")
print(res[order(res$delta, res$sigma), ], row.names = FALSE, digits = 3)
