# J1: exact consequences of the joint specification at q = 4, 5 (Section 5.3).
#
# For every graph at q = 4 (64 graphs) and q = 5 (1024 graphs), estimate
# Z(Gamma; delta) by the validated CRN importance sampler
# (R/src/monotonicity_check.cpp), then compute the joint-specification graph
# marginal p_J(Gamma) prop. to w(Gamma) Z(Gamma) exactly (up to MC error in Z)
# for independent-Bernoulli graph priors.
#
# Factors: eta in {0.5, 1, 2}; delta in {0, 1, 2, 0.5 log q};
#          p_inc in {0.1, 0.25, 0.5} (post-processing only).
# Two full replicate runs with independent seeds quantify MC error.
#
# Estimands: E[edges] under p_J vs under w; implied average per-edge odds
# factor; per-(graph, absent edge) conditional odds factors Z(G+e)/Z(G)
# checked against the bracket [c(0), 1).
#
# Output: R/results/joint_spec_sparsening_J1.rds, R/results/joint_spec_sparsening_J1_z.rds,
#         tables/joint_sparsening_enumeration.tex (enumeration-based table; the
#         manuscript's tables/joint_sparsening.tex is written by
#         joint_spec_sparsening_J1_chain_table.R from long chain runs, because
#         the importance-sampling constants are unreliable at delta = 2).

suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
})

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
sourceCpp(file.path(PROJECT_ROOT, "R", "src", "monotonicity_check.cpp"))
OUT_Z <- file.path(PROJECT_ROOT, "R", "results", "joint_spec_sparsening_J1_z.rds")
OUT_RDS <- file.path(PROJECT_ROOT, "R", "results", "joint_spec_sparsening_J1.rds")
OUT_TEX <- file.path(PROJECT_ROOT, "tables", "joint_sparsening_enumeration.tex")

etas <- c(0.5, 1, 2)
p_incs <- c(0.1, 0.25, 0.5)
seeds <- c(71001L, 71002L)
n_draws_for_q <- c(`4` = 5e5, `5` = 2e5)

edge_list <- function(q) t(utils::combn(q, 2))
graph_from_mask <- function(mask, edges, q) {
  G <- matrix(0L, q, q)
  for (k in seq_len(nrow(edges))) {
    if (bitwAnd(mask, bitwShiftL(1L, k - 1L)) != 0L) {
      G[edges[k, 1], edges[k, 2]] <- 1L
      G[edges[k, 2], edges[k, 1]] <- 1L
    }
  }
  G
}
edge_count <- function(mask) sum(bitwAnd(mask, bitwShiftL(1L, 0:30)) != 0L)

c0_fun <- function(delta, eta) {
  s_int <- function(s) {
    t <- (delta + 1) * s
    exp(delta * log(t) - t - lgamma(delta + 1) + log(delta + 1)) /
      sqrt(1 + 2 * eta^2 / t)
  }
  stats::integrate(s_int, 0, Inf, rel.tol = 1e-11)$value
}

## ---- Z for every graph, every setting, two seeds ----

z_rows <- list()
for (q in c(4L, 5L)) {
  edges <- edge_list(q)
  M <- nrow(edges)
  masks <- 0:(2^M - 1L)
  deltas <- c(0, 1, 2, 0.5 * log(q))
  n_draws <- n_draws_for_q[as.character(q)]
  graphs <- lapply(masks, graph_from_mask, edges = edges, q = q)
  for (eta in etas) for (delta in deltas) for (seed in seeds) {
    t0 <- Sys.time()
    z <- mclapply(seq_along(masks), function(i) {
      r <- z_graph_aligned_cpp(graphs[[i]], beta = eta, sigma = 1,
                               delta = delta, n_draws = n_draws, seed = seed)
      c(r$z_hat, r$z_se)
    }, mc.cores = 12L)
    z <- do.call(rbind, z)
    z_rows[[length(z_rows) + 1L]] <- data.frame(
      q = q, eta = eta, delta = delta, seed = seed, mask = masks,
      m = vapply(masks, edge_count, integer(1)),
      z_hat = z[, 1], z_se = z[, 2])
    cat(sprintf("q=%d eta=%.1f delta=%.2f seed=%d done in %.1f s\n",
                q, eta, delta, seed,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}
zt <- do.call(rbind, z_rows)
saveRDS(zt, OUT_Z)

## ---- post-processing ----

summarize_setting <- function(d, p_inc) {
  M <- max(d$m)
  logw <- d$m * log(p_inc) + (M - d$m) * log(1 - p_inc)
  pj <- exp(logw) * d$z_hat
  pj <- pj / sum(pj)
  E_J <- sum(d$m * pj)
  E_w <- M * p_inc
  pbar <- E_J / M
  odds_factor <- (pbar / (1 - pbar)) / (p_inc / (1 - p_inc))
  data.frame(q = d$q[1], eta = d$eta[1], delta = d$delta[1], seed = d$seed[1],
             p_inc = p_inc, E_w = E_w, E_J = E_J, ratio = E_J / E_w,
             odds_factor = odds_factor)
}

key <- interaction(zt$q, zt$eta, zt$delta, zt$seed, drop = TRUE)
summ <- do.call(rbind, lapply(split(zt, key), function(d) {
  do.call(rbind, lapply(p_incs, function(p) summarize_setting(d, p)))
}))
rownames(summ) <- NULL

# replicate-seed agreement -> MC uncertainty of each functional
agg_key <- interaction(summ$q, summ$eta, summ$delta, summ$p_inc, drop = TRUE)
summ_agg <- do.call(rbind, lapply(split(summ, agg_key), function(d) {
  data.frame(q = d$q[1], eta = d$eta[1], delta = d$delta[1], p_inc = d$p_inc[1],
             E_w = d$E_w[1], E_J = mean(d$E_J),
             E_J_err = abs(diff(d$E_J)) / 2,
             ratio = mean(d$ratio),
             odds_factor = mean(d$odds_factor),
             odds_err = abs(diff(d$odds_factor)) / 2,
             c0 = c0_fun(d$delta[1], d$eta[1]))
}))
rownames(summ_agg) <- NULL

## ---- bracket check on all per-edge conditional odds factors ----

bracket_check <- function(d_seedavg) {
  q <- d_seedavg$q[1]
  edges <- edge_list(q)
  M <- nrow(edges)
  z <- d_seedavg$z_hat[order(d_seedavg$mask)]
  se <- d_seedavg$z_se[order(d_seedavg$mask)]
  c0 <- c0_fun(d_seedavg$delta[1], d_seedavg$eta[1])
  n_pairs <- 0L; n_inside <- 0L; rmin <- Inf; rmax <- -Inf
  for (mask in 0:(2^M - 1L)) {
    for (k in seq_len(M)) {
      bit <- bitwShiftL(1L, k - 1L)
      if (bitwAnd(mask, bit) == 0L) {
        i0 <- mask + 1L; i1 <- mask + bit + 1L
        r <- z[i1] / z[i0]
        tol <- 4 * r * sqrt((se[i0] / z[i0])^2 + (se[i1] / z[i1])^2)
        n_pairs <- n_pairs + 1L
        if (r >= c0 - tol && r < 1 + tol) n_inside <- n_inside + 1L
        rmin <- min(rmin, r); rmax <- max(rmax, r)
      }
    }
  }
  data.frame(q = q, eta = d_seedavg$eta[1], delta = d_seedavg$delta[1],
             c0 = c0, n_pairs = n_pairs, n_inside = n_inside,
             r_min = rmin, r_max = rmax)
}

zt_avg_key <- interaction(zt$q, zt$eta, zt$delta, drop = TRUE)
zt_avg <- do.call(rbind, lapply(split(zt, zt_avg_key), function(d) {
  s1 <- d[d$seed == seeds[1], ]; s2 <- d[d$seed == seeds[2], ]
  s1 <- s1[order(s1$mask), ]; s2 <- s2[order(s2$mask), ]
  data.frame(q = s1$q, eta = s1$eta, delta = s1$delta, mask = s1$mask,
             z_hat = (s1$z_hat + s2$z_hat) / 2,
             z_se = sqrt(s1$z_se^2 + s2$z_se^2) / 2)
}))
brackets <- do.call(rbind, lapply(split(zt_avg, interaction(zt_avg$q, zt_avg$eta,
                                                            zt_avg$delta, drop = TRUE)),
                                  bracket_check))
rownames(brackets) <- NULL

saveRDS(list(summary = summ_agg, brackets = brackets),
        OUT_RDS)

cat("\n== sparsening summary (p_inc = 0.5) ==\n")
print(subset(summ_agg, p_inc == 0.5), row.names = FALSE, digits = 3)
cat("\n== sparsening summary (p_inc = 0.25) ==\n")
print(subset(summ_agg, p_inc == 0.25), row.names = FALSE, digits = 3)
cat("\n== bracket check over all (graph, absent edge) pairs ==\n")
print(brackets, row.names = FALSE, digits = 3)

## ---- manuscript table: q = 5, p_inc = 0.5 ----

tb <- subset(summ_agg, q == 5 & p_inc == 0.5)
tb <- tb[order(tb$eta, tb$delta), ]
lines <- c(
  "\\begin{tabular}{llrrrr}",
  "\\toprule",
  "$\\eta$ & $\\delta$ & $\\mathbb{E}_{\\mathrm{J}}[m]$ & $\\mathbb{E}_{\\mathrm{J}}[m]/\\mathbb{E}_{w}[m]$ & odds factor & $c(0)$ \\\\",
  "\\midrule")
for (i in seq_len(nrow(tb))) {
  dl <- if (abs(tb$delta[i] - 0.5 * log(5)) < 1e-9) "$\\tfrac12\\log q$" else sprintf("%.0f", tb$delta[i])
  lines <- c(lines, sprintf("%.1f & %s & %.2f & %.3f & %.3f & %.3f \\\\",
                            tb$eta[i], dl, tb$E_J[i], tb$ratio[i],
                            tb$odds_factor[i], tb$c0[i]))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, OUT_TEX)
cat("\nWrote tables/joint_sparsening_enumeration.tex\n")
