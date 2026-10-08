# Check V2 of the verification appendix: the three representations of the
# normalizing constant in the monotonicity appendix (brute force, the pair
# reduction N_delta/G_delta of Lemma lem:pair-mixture, written Ntilde/Gtilde
# below, and the conditional-Gaussian form A(tau^2) of Lemma lem:cond-gauss)
# agree within Monte Carlo error.
#
# Claims checked numerically:
#   (a) Gamma-mixture identity for the pair integral:
#       Ntilde(u) = Gamma(delta+1) beta^{-2(delta+1)} int t^delta e^{-t} e^{-beta^2 u^2 / t} dt
#       equals the defining 2-D integral int int_{s1 s2 > u^2} (s1 s2 - u^2)^delta e^{-beta(s1+s2)}.
#   (b) Three estimators of Z(G-) and Z(G+) agree (common random numbers):
#       (1) brute force over all entries (Schur form, spot-checked against dense eigen),
#       (2) pair reduction: E[ w2 * Ntilde(U12) ] and E[ w2 * Gtilde(U12) ],
#       (3) conditional-Gaussian: E[ w3 * A(tau^2) ] and E[ w3 * A(tau^2 + sigma^2) ],
#       with A(v) = pref * int t^delta e^{-t} sqrt(t / (t + 2 beta^2 v)) dt.
#   (c) c(0) = A(sigma^2)/A(0) reproduces the known table (0.4918 / 0.629 / 0.709
#       at eta = 1 for delta = 0 / 0.8 / 1.6).
#   (d) Disconnected pair: estimator (3) gives the ratio c(0) exactly (tau^2 == 0).
#   (e) log A(v) is convex in v (finite second differences >= 0).
#
# Graphs at q = 4, absent edge e = (1,2), R = {3,4}:
#   A: edges {1-3, 2-3, 3-4}  (endpoints connected through node 3)
#   B: edges {1-3, 2-4}       (endpoints disconnected)
#   C: edges {1-3, 2-3, 2-4, 3-4} (two free coordinates in the b column)

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
set.seed(20260728)

beta_par <- 1
sigma_par <- 1
deltas <- c(0, 0.8, 1.6)
n_mc <- 3e5
known_c0 <- c(`0` = 0.4918, `0.8` = 0.629, `1.6` = 0.709)

## ---- Gauss-Legendre nodes (Golub-Welsch), used after t = u^2 substitution ----

gauss_legendre <- function(n, a, b) {
  i <- seq_len(n - 1)
  off <- i / sqrt(4 * i^2 - 1)
  jac <- matrix(0, n, n)
  jac[cbind(i, i + 1)] <- off
  jac[cbind(i + 1, i)] <- off
  eig <- eigen(jac, symmetric = TRUE)
  list(nodes = (b - a) / 2 * eig$values + (a + b) / 2,
       weights = (b - a) / 2 * 2 * eig$vectors[1, ]^2)
}

# int_0^inf t^delta e^{-t} g(t) dt  via t = x^2 on x in (0, 9)
gl <- gauss_legendre(400, 0, 9)
t_nodes <- gl$nodes^2

quad_weights <- function(delta) {
  2 * gl$weights * gl$nodes^(2 * delta + 1) * exp(-gl$nodes^2)
}

pref <- function(delta) gamma(delta + 1) * beta_par^(-2 * (delta + 1))

n_tilde <- function(u, delta, w_t) {
  # returns Ntilde(u) for a vector u
  kern <- exp(-beta_par^2 * outer(1 / t_nodes, u^2))
  pref(delta) * as.vector(crossprod(w_t, kern))
}

g_tilde <- function(u, delta, w_t) {
  shift <- 2 * beta_par^2 * sigma_par^2
  deficit <- sqrt(t_nodes / (t_nodes + shift))
  kern <- exp(-beta_par^2 * outer(1 / (t_nodes + shift), u^2))
  pref(delta) * as.vector(crossprod(w_t * deficit, kern))
}

a_fun <- function(v, delta, w_t) {
  # A(v) = pref * sum w_t sqrt(t / (t + 2 beta^2 v)), vectorized over v
  kern <- sqrt(t_nodes / outer(t_nodes, 2 * beta_par^2 * v, `+`))
  pref(delta) * as.vector(crossprod(w_t, kern))
}

## ---- (a) mixture identity vs defining 2-D integral ----

n_tilde_2d <- function(u, delta) {
  inner <- function(s1) {
    vapply(s1, function(s) {
      stats::integrate(function(s2) (s * s2 - u^2)^delta * exp(-beta_par * (s + s2)),
                       lower = u^2 / s, upper = Inf, rel.tol = 1e-10)$value
    }, numeric(1))
  }
  stats::integrate(inner, lower = 0, upper = Inf, rel.tol = 1e-8)$value
}

cat("== (a) Gamma-mixture identity for Ntilde ==\n")
for (delta in deltas) {
  w_t <- quad_weights(delta)
  for (u in c(0, 0.3, 1, 2)) {
    v1 <- n_tilde(u, delta, w_t)
    v2 <- n_tilde_2d(u, delta)
    cat(sprintf("delta = %.1f  u = %.1f  mixture = %.8f  2-D = %.8f  rel.err = %.2e\n",
                delta, u, v1, v2, abs(v1 / v2 - 1)))
  }
}

## ---- (c) c(0) against the known table ----

cat("\n== (c) c(0) = A(sigma^2)/A(0) vs known values (eta = 1) ==\n")
for (delta in deltas) {
  w_t <- quad_weights(delta)
  c0 <- a_fun(sigma_par^2, delta, w_t) / a_fun(0, delta, w_t)
  cat(sprintf("delta = %.1f  c(0) = %.4f  (known %.4f)\n",
              delta, c0, known_c0[[as.character(delta)]]))
}

## ---- (e) log-convexity of A ----

cat("\n== (e) log-convexity of A(v): min second difference of log A ==\n")
v_grid <- seq(0, 5, by = 0.05)
for (delta in deltas) {
  w_t <- quad_weights(delta)
  la <- log(a_fun(v_grid, delta, w_t))
  d2 <- diff(la, differences = 2)
  cat(sprintf("delta = %.1f  min d2 log A = %.3e  (must be >= 0)\n", delta, min(d2)))
}

## ---- graph definitions ----

# rest R = {3,4}; a = (k13, k14), b = (k23, k24); node 3 is index 1 in R
graphs <- list(
  A = list(edge_34 = TRUE,  a_free = c(TRUE, FALSE), b_free = c(TRUE, FALSE)),
  B = list(edge_34 = FALSE, a_free = c(TRUE, FALSE), b_free = c(FALSE, TRUE)),
  C = list(edge_34 = TRUE,  a_free = c(TRUE, FALSE), b_free = c(TRUE, TRUE))
)

## ---- (b) three estimators, common random numbers ----

run_graph <- function(g, delta, w_t) {
  k33 <- stats::rexp(n_mc, beta_par)
  k44 <- stats::rexp(n_mc, beta_par)
  k34 <- if (g$edge_34) stats::rnorm(n_mc, 0, sigma_par) else numeric(n_mc)
  a1 <- if (g$a_free[1]) stats::rnorm(n_mc, 0, sigma_par) else numeric(n_mc)
  a2 <- if (g$a_free[2]) stats::rnorm(n_mc, 0, sigma_par) else numeric(n_mc)
  b1 <- if (g$b_free[1]) stats::rnorm(n_mc, 0, sigma_par) else numeric(n_mc)
  b2 <- if (g$b_free[2]) stats::rnorm(n_mc, 0, sigma_par) else numeric(n_mc)
  k11 <- stats::rexp(n_mc, beta_par)
  k22 <- stats::rexp(n_mc, beta_par)
  x_edge <- stats::rnorm(n_mc, 0, sigma_par)

  det_kr <- k33 * k44 - k34^2
  pd_kr <- det_kr > 0                       # k33, k44 > 0 always
  # M = K_R^{-1}, explicit 2x2
  m11 <- k44 / det_kr; m12 <- -k34 / det_kr; m22 <- k33 / det_kr
  u11 <- a1^2 * m11 + 2 * a1 * a2 * m12 + a2^2 * m22
  u22 <- b1^2 * m11 + 2 * b1 * b2 * m12 + b2^2 * m22
  u12 <- a1 * (m11 * b1 + m12 * b2) + a2 * (m12 * b1 + m22 * b2)

  tilt_kr <- ifelse(pd_kr, det_kr^delta, 0)

  # (1) brute force via Schur form
  s11m <- k11 - u11; s22m <- k22 - u22
  det_sm <- s11m * s22m - u12^2
  det_sp <- s11m * s22m - (x_edge - u12)^2
  z1_minus <- ifelse(pd_kr & s11m > 0 & det_sm > 0, tilt_kr * det_sm^delta, 0)
  z1_plus <- ifelse(pd_kr & s11m > 0 & det_sp > 0, tilt_kr * det_sp^delta, 0)

  # (2) pair reduction
  w2 <- tilt_kr * beta_par^2 * exp(-beta_par * (u11 + u22))
  w2[!pd_kr] <- 0
  z2_minus <- w2 * n_tilde(u12, delta, w_t)
  z2_plus <- w2 * g_tilde(u12, delta, w_t)

  # (3) conditional-Gaussian: integrate the b column analytically
  bf <- which(g$b_free)
  ma1 <- m11 * a1 + m12 * a2               # (M a) at node 3
  ma2 <- m12 * a1 + m22 * a2               # (M a) at node 4
  if (length(bf) == 1) {
    mf <- if (bf == 1) m11 else m22
    cc <- if (bf == 1) ma1 else ma2
    dfac <- 1 / sqrt(1 + 2 * beta_par * sigma_par^2 * mf)
    tau2 <- sigma_par^2 * cc^2 / (1 + 2 * beta_par * sigma_par^2 * mf)
  } else {
    g11 <- 1 + 2 * beta_par * sigma_par^2 * m11
    g12 <- 2 * beta_par * sigma_par^2 * m12
    g22 <- 1 + 2 * beta_par * sigma_par^2 * m22
    det_g <- g11 * g22 - g12^2
    dfac <- 1 / sqrt(det_g)
    tau2 <- sigma_par^2 * (g22 * ma1^2 - 2 * g12 * ma1 * ma2 + g11 * ma2^2) / det_g
  }
  tau2[!pd_kr] <- 0        # tau2/dfac can be NaN on non-PD rest draws; those get weight 0
  dfac[!pd_kr] <- 0
  w3 <- tilt_kr * beta_par^2 * exp(-beta_par * u11) * dfac
  w3[!pd_kr] <- 0
  z3_minus <- w3 * a_fun_vec(tau2, delta, w_t)
  z3_plus <- w3 * a_fun_vec(tau2 + sigma_par^2, delta, w_t)

  list(z_minus = cbind(z1_minus, z2_minus, z3_minus),
       z_plus = cbind(z1_plus, z2_plus, z3_plus),
       spot = list(k33 = k33, k44 = k44, k34 = k34, a1 = a1, a2 = a2,
                   b1 = b1, b2 = b2, k11 = k11, k22 = k22, x_edge = x_edge))
}

# chunked A(v) over a long vector of draws (memory)
a_fun_vec <- function(v, delta, w_t) {
  out <- numeric(length(v))
  idx <- split(seq_along(v), ceiling(seq_along(v) / 2e4))
  for (ii in idx) out[ii] <- a_fun(v[ii], delta, w_t)
  out
}

# dense-matrix spot check of the Schur brute force (independent of all algebra above)
spot_check <- function(sp, g, delta, n_spot = 2000) {
  max_err <- 0
  mismatch <- 0L
  for (r in seq_len(n_spot)) {
    kk <- matrix(0, 4, 4)
    diag(kk) <- c(sp$k11[r], sp$k22[r], sp$k33[r], sp$k44[r])
    kk[1, 3] <- kk[3, 1] <- sp$a1[r]; kk[1, 4] <- kk[4, 1] <- sp$a2[r]
    kk[2, 3] <- kk[3, 2] <- sp$b1[r]; kk[2, 4] <- kk[4, 2] <- sp$b2[r]
    kk[3, 4] <- kk[4, 3] <- sp$k34[r]
    for (plus in c(FALSE, TRUE)) {
      kk[1, 2] <- kk[2, 1] <- if (plus) sp$x_edge[r] else 0
      ev <- eigen(kk, symmetric = TRUE, only.values = TRUE)$values
      val_dense <- if (min(ev) > 0) prod(ev)^delta else 0
      # recompute the Schur-form value for this draw
      m <- solve(kk[3:4, 3:4])
      uu <- kk[1:2, 3:4] %*% m %*% kk[3:4, 1:2]
      ss <- kk[1:2, 1:2] - uu
      pd <- min(eigen(kk[3:4, 3:4], symmetric = TRUE, only.values = TRUE)$values) > 0 &&
        ss[1, 1] > 0 && det(ss) > 0
      val_schur <- if (pd) (det(kk[3:4, 3:4]) * det(ss))^delta else 0
      if ((val_dense == 0) != (val_schur == 0)) mismatch <- mismatch + 1L
      if (val_dense > 0) max_err <- max(max_err, abs(val_schur / val_dense - 1))
    }
  }
  cat(sprintf("    spot check (%d draws): PD flag mismatches = %d, max rel err = %.2e\n",
              n_spot, mismatch, max_err))
}

cat("\n== (b) three estimators of Z(G-) and Z(G+), n =", format(n_mc, big.mark = ","), "==\n")
results <- list()
for (gname in names(graphs)) {
  for (delta in deltas) {
    w_t <- quad_weights(delta)
    res <- run_graph(graphs[[gname]], delta, w_t)
    est_m <- colMeans(res$z_minus)
    se_m <- apply(res$z_minus, 2, stats::sd) / sqrt(n_mc)
    est_p <- colMeans(res$z_plus)
    se_p <- apply(res$z_plus, 2, stats::sd) / sqrt(n_mc)
    ratio <- est_p / est_m
    c0 <- a_fun(sigma_par^2, delta, w_t) / a_fun(0, delta, w_t)
    cat(sprintf("\ngraph %s, delta = %.1f  [c(0) = %.4f]\n", gname, delta, c0))
    for (k in 1:3) {
      cat(sprintf("  est %d: Z- = %.6f (%.6f)   Z+ = %.6f (%.6f)   ratio = %.4f\n",
                  k, est_m[k], se_m[k], est_p[k], se_p[k], ratio[k]))
    }
    if (gname == "B") {
      cat(sprintf("  disconnected pair: est-3 ratio - c(0) = %.2e (must be ~0 exactly)\n",
                  ratio[3] - c0))
    }
    if (delta == deltas[1]) spot_check(res$spot, graphs[[gname]], delta)
    results[[paste(gname, delta, sep = "_")]] <- list(
      graph = gname, delta = delta, est_minus = est_m, se_minus = se_m,
      est_plus = est_p, se_plus = se_p, ratio = ratio, c0 = c0)
  }
}

saveRDS(results, file.path(PROJECT_ROOT, "R", "results", "verify_conditional_gaussian_route.rds"))
cat("\nSaved: R/results/verify_conditional_gaussian_route.rds\n")
