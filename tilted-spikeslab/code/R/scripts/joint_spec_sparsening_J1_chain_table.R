# Chain-based regeneration of Table tab:joint-sparsening (q = 5, p_inc = 0.5)
# and a direct measurement of the compensated specification.
#
# Motivation: the J1 brute-force normalizing constants are importance-sampling
# estimates with heavy-tailed weights at delta = 2; a 2e6-draw recheck
# (R/results/j1_delta2_recheck.log; the recheck was run interactively and its
# script was not archived) shows between-seed half-differences of 0.06 to 0.25
# expected edges there (0.06 to 0.09 at p_inc = 0.25, 0.14 to 0.25 at p_inc =
# 0.5), while the verified chain agrees between seeds
# to +-0.01. The table is therefore regenerated from long runs of the
# verified C++ chain (verify_sampler_port.R), one method for all cells.
#
# Compensation: instead of enumerating with noisy Z, run the chain directly at
# the compensated inclusion probability p_c with p_c/(1 - p_c) = 1/c(0)
# (intended odds 1 at p_inc = 0.5) and measure the recovered E[m] against 5.
#
# Output: R/results/joint_spec_sparsening_J1_chain_table.rds (+ log),
#         tables/joint_sparsening.tex

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})
PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
sourceCpp(file.path(PROJECT_ROOT, "R", "src", "spikeslab_sampler.cpp"))

c0_fun <- function(delta, eta) {
  s_int <- function(s) {
    t <- (delta + 1) * s
    exp(delta * log(t) - t - lgamma(delta + 1) + log(delta + 1)) /
      sqrt(1 + 2 * eta^2 / t)
  }
  stats::integrate(s_int, 0, Inf, rel.tol = 1e-11)$value
}

q <- 5L
p_inc <- 0.5
etas <- c(0.5, 1, 2)
deltas <- c(0, 0.5 * log(q), 1, 2)
seeds <- c(7601L, 7602L)
n_iter <- 200000L
n_burn <- 2000L

cells <- list()
for (eta in etas) for (delta in deltas) for (arm in c("plain", "compensated"))
  for (seed in seeds) {
    cells[[length(cells) + 1L]] <- list(eta = eta, delta = delta, arm = arm,
                                        seed = seed)
  }
cat(sprintf("running %d chains on 12 cores\n", length(cells)))

run_cell <- function(cell) {
  set.seed(cell$seed)
  c0 <- c0_fun(cell$delta, cell$eta)
  p <- if (cell$arm == "plain") p_inc else 1 / (1 + c0)
  ch <- sampler_chain_cpp(S = matrix(0, q, q), n_obs = 0, q = q,
                          delta = cell$delta, beta = cell$eta, sigma = 1,
                          p_inc = p, spec = 0L, log_z = numeric(0),
                          n_iter = n_iter, n_burn = n_burn,
                          refresh_every = 1L, track_masks = FALSE)
  data.frame(eta = cell$eta, delta = cell$delta, arm = cell$arm,
             seed = cell$seed, p_used = p, E_m = ch$E_m,
             drift = ch$max_sigma_drift)
}

t0 <- Sys.time()
out <- mclapply(cells, run_cell, mc.cores = 12L, mc.preschedule = FALSE)
res <- do.call(rbind, out)
cat(sprintf("done in %.1f min; max drift %.2e\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            max(res$drift)))

key <- interaction(res$eta, res$delta, res$arm, drop = TRUE)
agg <- do.call(rbind, lapply(split(res, key), function(d) {
  M <- q * (q - 1) / 2
  E_J <- mean(d$E_m)
  pbar <- E_J / M
  data.frame(eta = d$eta[1], delta = d$delta[1], arm = d$arm[1],
             E_J = E_J, halfdiff = abs(diff(d$E_m))[1] / 2,
             ratio = E_J / (M * p_inc),
             odds_factor = (pbar / (1 - pbar)) / (p_inc / (1 - p_inc)),
             c0 = c0_fun(d$delta[1], d$eta[1]))
}))
rownames(agg) <- NULL
saveRDS(agg, file.path(PROJECT_ROOT, "R", "results",
                       "joint_spec_sparsening_J1_chain_table.rds"))

plain <- agg[agg$arm == "plain", ]
plain <- plain[order(plain$eta, plain$delta), ]
comp <- agg[agg$arm == "compensated", ]
comp <- comp[order(comp$eta, comp$delta), ]

cat("\n== plain joint specification (table cells) ==\n")
print(plain[, c("eta", "delta", "E_J", "halfdiff", "ratio", "odds_factor", "c0")],
      row.names = FALSE, digits = 3)
cat(sprintf("max halfdiff (plain): %.4f\n", max(plain$halfdiff)))

cat("\n== compensated specification: recovered E[m] against 5 ==\n")
comp$rel_err <- comp$E_J / 5 - 1
print(comp[, c("eta", "delta", "E_J", "halfdiff", "rel_err")],
      row.names = FALSE, digits = 3)

## ---- write the table ----
dlab <- function(delta) {
  if (abs(delta - 0.5 * log(q)) < 1e-9) "$\\tfrac12\\log q$"
  else sprintf("%g", delta)
}
lines <- c("\\begin{tabular}{llrrrr}", "\\toprule",
  "$\\eta$ & $\\delta$ & $\\mathbb{E}_{\\mathrm{J}}[m]$ & $\\mathbb{E}_{\\mathrm{J}}[m]/\\mathbb{E}_{w}[m]$ & odds factor & $c(0)$ \\\\",
  "\\midrule")
for (eta in etas) for (delta in deltas) {
  d <- plain[plain$eta == eta & abs(plain$delta - delta) < 1e-9, ]
  lines <- c(lines, sprintf("%.1f & %s & %.2f & %.3f & %.3f & %.3f \\\\",
                            eta, dlab(delta), d$E_J, d$ratio, d$odds_factor,
                            d$c0))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(PROJECT_ROOT, "tables", "joint_sparsening.tex"))
cat("\nSaved tables/joint_sparsening.tex and R/results/joint_spec_sparsening_J1_chain_table.rds\n")
