# R/scripts/rho_uniform_calibration.R
#
# What specification yields an approximately UNIFORM prior on the partial
# correlation rho_{ij} | gamma_{ij} = 1? Two levers:
#   - delta : determinant tilt (G-Wishart df = 2 delta + 2)
#   - sigma : slab standard deviation (the narrow slab shrinks rho toward 0)
#
# For each (q, delta, sigma) we draw on-edge rho and summarise its shape by
#   (i)  SD(rho)                 [Uniform(-1,1) has SD = 1/sqrt(3) = 0.5774]
#   (ii) beta = MLE exponent of the symmetric density  f(rho) ∝ (1 - rho^2)^beta
#        beta = 0  <=> uniform;  beta > 0 <=> peaked at 0;  beta < 0 <=> piled at +-1
#   (iii) Pr(|rho| > 0.5)        [Uniform value = 0.5]
#
# Output: R/results/rho_uniform_calibration.rds

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
SRC_CPP <- file.path(PROJECT_ROOT, "R", "src", "Z_prior_gwishart_slab_rho.cpp")
OUT_RDS <- file.path(PROJECT_ROOT, "R", "results", "rho_uniform_calibration.rds")
sourceCpp(SRC_CPP)

set.seed(20260724L)
BETA <- 1

qs     <- c(10L, 20L)
deltas <- c(0, 0.5, 1, 2, 4.5)
sigmas <- c(1, 2, 4, 8, 20)

p_inc_for_q <- c(`10` = 0.20, `20` = 0.10)
n_target    <- 4000L
cap_for_q   <- c(`10` = 4000000L, `20` = 60000000L)

# MLE of beta in f(rho) ∝ (1 - rho^2)^beta on (-1,1).
# logZ(beta) = log(sqrt(pi)) + lgamma(beta+1) - lgamma(beta+1.5)
fit_beta <- function(rho) {
  rho <- rho[is.finite(rho) & abs(rho) < 1]
  S <- sum(log1p(-rho^2)); n <- length(rho)
  negll <- function(b) -(b * S - n * (0.5 * log(pi) + lgamma(b + 1) - lgamma(b + 1.5)))
  opt <- optimize(negll, interval = c(-0.95, 40))
  opt$minimum
}

cells <- list(); sc <- 1L
for (q in qs) for (d in deltas) for (s in sigmas) {
  cells[[length(cells)+1L]] <- list(q=q, delta=d, sigma=s,
    p_inc=p_inc_for_q[as.character(q)], cap=cap_for_q[as.character(q)], seed=9000L+sc)
  sc <- sc + 1L
}
cat(sprintf("Cells: %d\n", length(cells)))

run_one <- function(cell) {
  res <- gwishart_slab_record_rho_cpp(
    q=cell$q, delta=cell$delta, beta=BETA, sigma=cell$sigma,
    gamma_prob=cell$p_inc, n_target=n_target, n_proposal_cap=cell$cap,
    max_inner_tries=10000L, seed=cell$seed, condition_on_edge_01=TRUE)
  on <- res$gamma12 == 1
  rho <- -res$k12[on] / sqrt(res$k11[on] * res$k22[on])
  rho <- rho[is.finite(rho)]
  data.frame(q=cell$q, delta=cell$delta, sigma=cell$sigma, n=length(rho),
             rho_sd=sd(rho), beta_hat=fit_beta(rho),
             p_gt_0.5=mean(abs(rho) > 0.5),
             n_acc=res$n_accepted, n_prop=res$n_proposals)
}

t0 <- Sys.time()
out <- mclapply(cells, run_one, mc.cores=12L, mc.preschedule=FALSE)
cat(sprintf("Done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))
tab <- do.call(rbind, out)
saveRDS(list(summary=tab, params=list(beta=BETA)), OUT_RDS)
cat(sprintf("Wrote %s\n", OUT_RDS))
cat("\nUniform target: SD = 0.577, beta = 0, Pr(|rho|>0.5) = 0.5\n\n")
print(tab[order(tab$q, tab$delta, tab$sigma),
          c("q","delta","sigma","n","rho_sd","beta_hat","p_gt_0.5")],
      row.names=FALSE, digits=3)
