# Verification of the coupled ratio estimator of Section sec:hier-ratio
# against (a) the analytic disconnected-pair value 1/c(0) at q = 2, (b) the
# brute-force Z ratios of the J1 table at q = 4 (delta <= 1, where that
# reference is reliable), and (c) a coupled-versus-independent variance and
# correlation comparison.
#
# Output: R/results/verify_akm_ratio.rds (+ log)

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

## ---- (a) q = 2 disconnected pair: J = 1/c(0) analytically ----
cat("== (a) q = 2 empty graph + edge (1,2): J against 1/c(0) ==\n")
res_a <- list()
for (eta in c(0.5, 1, 2)) for (delta in c(0, 1, 2)) {
  set.seed(3301L)
  J_reps <- vapply(1:32, function(r)
    akm_ratio_cpp(matrix(0L, 2, 2), 1L, 2L, delta, eta, 1, 100000L)$J,
    numeric(1))
  J_hat <- mean(J_reps)
  J_se <- sd(J_reps) / sqrt(32)
  J_true <- 1 / c0_fun(delta, eta)
  z <- (J_hat - J_true) / J_se
  res_a[[length(res_a) + 1L]] <- data.frame(
    eta = eta, delta = delta, J_hat = J_hat, J_se = J_se, J_true = J_true,
    z_score = z)
  cat(sprintf("eta=%.1f delta=%d: J_hat = %.5f (se %.5f), 1/c0 = %.5f, z = %+.2f\n",
              eta, delta, J_hat, J_se, J_true, z))
}
res_a <- do.call(rbind, res_a)

## ---- (b) q = 4 against brute-force Z ratios (delta <= 1) ----
cat("\n== (b) q = 4 against brute-force Z ratios ==\n")
zt <- readRDS(file.path(PROJECT_ROOT, "R", "results", "joint_spec_sparsening_J1_z.rds"))
edges4 <- t(utils::combn(4, 2))                # bit k-1 <-> edge row k
mask_of <- function(edge_rows) sum(bitwShiftL(1L, edge_rows - 1L))
graph_of <- function(edge_rows) {
  G <- matrix(0L, 4, 4)
  for (k in edge_rows) {
    G[edges4[k, 1], edges4[k, 2]] <- 1L
    G[edges4[k, 2], edges4[k, 1]] <- 1L
  }
  G
}
z_lookup <- function(eta, delta, mask) {
  d <- zt[zt$q == 4 & zt$eta == eta & abs(zt$delta - delta) < 1e-9 &
            zt$mask == mask, ]
  c(z = mean(d$z_hat), se = sqrt(sum(d$z_se^2)) / nrow(d))
}
# (graph, edge) test pairs: empty + (1,2); path 1-2-3 + (3,4) [disconnected
# endpoints]; path 1-2-3 + (1,3) [cycle-closing]; star at 2 + (3,4)
pairs <- list(
  list(label = "empty + (1,2)",        rows = integer(0),  i = 1L, j = 2L),
  list(label = "path12,23 + (3,4)",    rows = c(1L, 4L),   i = 3L, j = 4L),
  list(label = "path12,23 + (1,3)",    rows = c(1L, 4L),   i = 1L, j = 3L),
  list(label = "star2 + (3,4)",        rows = c(1L, 4L, 5L), i = 3L, j = 4L)
)
cells <- list()
for (eta in c(1, 2)) for (delta in c(0, 1)) for (p in pairs)
  cells[[length(cells) + 1L]] <- c(p, eta = eta, delta = delta)
run_b <- function(cell) {
  set.seed(3401L)
  Gm <- graph_of(cell$rows)
  edge_row <- which(edges4[, 1] == cell$i & edges4[, 2] == cell$j)
  m_minus <- mask_of(cell$rows)
  m_plus <- m_minus + bitwShiftL(1L, edge_row - 1L)
  zm <- z_lookup(cell$eta, cell$delta, m_minus)
  zp <- z_lookup(cell$eta, cell$delta, m_plus)
  J_true <- zm["z"] / zp["z"]
  J_true_se <- J_true * sqrt((zm["se"] / zm["z"])^2 + (zp["se"] / zp["z"])^2)
  J_reps <- vapply(1:8, function(r)
    akm_ratio_cpp(Gm, cell$i, cell$j, cell$delta, cell$eta, 1, 200000L)$J,
    numeric(1))
  data.frame(label = cell$label, eta = cell$eta, delta = cell$delta,
             J_hat = mean(J_reps), J_se = sd(J_reps) / sqrt(8),
             J_true = unname(J_true), J_true_se = unname(J_true_se))
}
res_b <- do.call(rbind, mclapply(cells, run_b, mc.cores = 12L,
                                 mc.preschedule = FALSE))
res_b$z_score <- (res_b$J_hat - res_b$J_true) /
  sqrt(res_b$J_se^2 + res_b$J_true_se^2)
print(res_b, row.names = FALSE, digits = 4)
cat(sprintf("max |z| = %.2f over %d cells\n", max(abs(res_b$z_score)),
            nrow(res_b)))

## ---- (c) coupled vs independent: variance and correlation ----
cat("\n== (c) coupled vs independent (path12,23,34 + (1,3), eta = 1, delta = 1) ==\n")
Gm <- graph_of(c(1L, 4L, 6L))                  # edges (1,2), (2,3), (3,4)
run_c <- function(coupled, R = 400L, N = 1000L) {
  set.seed(if (coupled) 3501L else 3502L)
  reps <- t(vapply(seq_len(R), function(r) {
    est <- akm_ratio_cpp(Gm, 1L, 3L, 1, 1, 1, N, coupled = coupled)
    c(est$log_J, est$mean_minus, est$mean_plus)
  }, numeric(3)))
  list(var_logJ = var(reps[, 1]), corr = cor(reps[, 2], reps[, 3]))
}
cc <- run_c(TRUE); ci <- run_c(FALSE)
cat(sprintf("coupled:     var(log J) = %.3e, corr(means) = %.3f\n",
            cc$var_logJ, cc$corr))
cat(sprintf("independent: var(log J) = %.3e, corr(means) = %.3f\n",
            ci$var_logJ, ci$corr))
cat(sprintf("variance reduction factor = %.1f\n", ci$var_logJ / cc$var_logJ))

saveRDS(list(a = res_a, b = res_b,
             c = data.frame(arm = c("coupled", "independent"),
                            var_logJ = c(cc$var_logJ, ci$var_logJ),
                            corr = c(cc$corr, ci$corr))),
        file.path(PROJECT_ROOT, "R", "results", "verify_akm_ratio.rds"))
cat("\nSaved R/results/verify_akm_ratio.rds\n")
