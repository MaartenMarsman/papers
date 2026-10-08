# R/scripts/verify-scale-injection.R
# SUPERSEDED (pre-pivot): tests the scale factorization of Z with an eigenvalue
# floor eps and a Gamma(3, .) diagonal, both dropped at the 2026-07-18 pivot.
# The floor-free factorization is Eq. eq:Z-factorization of the manuscript.
# Purpose: Verify the congruence factorisation of the prior normaliser Z for the
#          determinant-tilted spike-and-slab GGM prior:
#            (1) pointwise change-of-variables identity
#                  scaled_integrand(D K D) * |Jac| == |A|^delta * std_integrand(K),
#            (2) floor must be anchored to K (floor-on-Theta breaks it for non-isotropic D),
#            (3) between-model ratio: Z_scaled(G')/Z_scaled(G) == Z*_std(G')/Z*_std(G),
#                and Z_scaled(G)/Z*_std(G) == |A|^delta (graph-independent).
# Run:     Rscript R/scripts/verify-scale-injection.R
# Notes:   base R only; discrete (point-mass) variant, q = 3.
# See:     notes/scale-injection-congruence.md

set.seed(20260616)

## ---- parameters ----
q      <- 3
alpha  <- 3.0                  # diagonal Gamma shape (shared)
b      <- c(1.0, 1.5, 2.0)     # STANDARDISED diagonal rates b_l  (rate absorbs c_l^2)
delta  <- 0.7                  # determinant tilt
eps    <- 0.05                 # eigenvalue floor (in the STANDARDISED frame)

## two graphs on the same node set: complete, and one non-edge (1,3)
edges_full <- list(c(1,2), c(1,3), c(2,3))
edges_sub  <- list(c(1,2),          c(2,3))   # theta_13 = 0

## two congruence scalings D = diag(c):
c_const <- rep(sqrt(2), q)         # constant sigma = 2  -> isotropic D = sqrt(2) I
c_gprior <- c(1.2, 0.8, 1.5)       # non-scalar diagonal (g-prior style, distinct c_l)

## ---- helpers ----
build_K <- function(diagv, offmap, edges) {
  K <- diag(diagv)
  for (e in edges) { K[e[1], e[2]] <- offmap[[paste(e, collapse="-")]]; K[e[2], e[1]] <- K[e[1], e[2]] }
  K
}
lam_min <- function(M) min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)

# unnormalised integrand (density factors); floor indicator evaluated on `floor_mat`
integrand <- function(M, diagv_rate, off_list, edges, floor_mat, eps) {
  if (lam_min(floor_mat) <= eps) return(0)    # cone + floor indicator, anchored to floor_mat
  val <- det(M)^delta
  for (l in 1:q) val <- val * dgamma(M[l, l], shape = alpha, rate = diagv_rate[l])
  for (e in edges) {
    key <- paste(e, collapse = "-")
    val <- val * dnorm(M[e[1], e[2]], mean = 0, sd = off_list[[key]])
  }
  val
}

## =====================================================================
## (1) + (2): pointwise identity and floor-anchoring, coupled draws
## =====================================================================
pointwise_check <- function(edges, cvec, Ntry = 30000) {
  D <- diag(cvec); A_det <- prod(cvec^2)       # |A| = prod c_l^2
  beta <- b / cvec^2                            # real diagonal rates
  sig  <- setNames(vector("list", length(edges)), sapply(edges, paste, collapse="-"))
  for (e in edges) sig[[paste(e, collapse="-")]] <- cvec[e[1]] * cvec[e[2]]   # sd = c_i c_j
  b_std <- setNames(as.list(rep(1, length(edges))), sapply(edges, paste, collapse="-")) # unit slab
  Jac <- prod(cvec^2) * prod(sapply(edges, function(e) cvec[e[1]] * cvec[e[2]]))

  max_rel <- 0; n_kept <- 0; n_mismatch_floorTheta <- 0
  for (it in 1:Ntry) {
    diagv <- rgamma(q, shape = alpha, rate = b)                  # k_ll ~ Gamma(alpha, b_l)
    offmap <- setNames(as.list(rnorm(length(edges))), sapply(edges, paste, collapse="-"))
    K <- build_K(diagv, offmap, edges)
    if (lam_min(K) <= eps) next                                  # keep PD & floor-passing K
    n_kept <- n_kept + 1
    Theta <- D %*% K %*% D

    std  <- integrand(K, b, b_std, edges, floor_mat = K, eps = eps)            # floor on K
    scal <- integrand(Theta, beta, sig, edges, floor_mat = K, eps = eps)        # floor on K (anchored)
    lhs <- scal * Jac
    rhs <- A_det^delta * std
    max_rel <- max(max_rel, abs(lhs / rhs - 1))

    # floor anchored to Theta instead (should break for non-isotropic D)
    if (lam_min(Theta) <= eps) n_mismatch_floorTheta <- n_mismatch_floorTheta + 1
  }
  list(A_det = A_det, n_kept = n_kept, max_rel = max_rel,
       n_mismatch_floorTheta = n_mismatch_floorTheta)
}

cat("=== (1)/(2) Pointwise identity  scaled*Jac == |A|^delta * std ===\n")
for (g in list(list(nm="full", E=edges_full), list(nm="sub(no 1-3)", E=edges_sub))) {
  for (d in list(list(nm="constant (isotropic)", c=c_const),
                 list(nm="g-prior  (non-scalar)", c=c_gprior))) {
    r <- pointwise_check(g$E, d$c)
    cat(sprintf("  graph=%-12s D=%-22s |A|=%.5f  kept=%5d  max|ratio-1|=%.3e  floor-on-Theta mismatches=%d\n",
                g$nm, d$nm, r$A_det, r$n_kept, r$max_rel, r$n_mismatch_floorTheta))
  }
}

## =====================================================================
## (3) integral level: independent MC for ratio cancellation
## =====================================================================
# Z*_std(G) = E_{k ~ std prior}[ det(K)^delta * 1{lam_min(K) > eps} ]
mc_Zstd <- function(edges, N = 4e5) {
  acc <- 0
  diagM <- matrix(rgamma(N * q, shape = alpha, rate = rep(b, each = N)), nrow = N)
  offM  <- matrix(rnorm(N * length(edges)), nrow = N)
  for (i in 1:N) {
    offmap <- setNames(as.list(offM[i, ]), sapply(edges, paste, collapse="-"))
    K <- build_K(diagM[i, ], offmap, edges)
    if (lam_min(K) > eps) acc <- acc + det(K)^delta
  }
  acc / N
}
# Z_scaled(G) via independent theta-draws, floor anchored to K = D^-1 Theta D^-1
mc_Zscaled <- function(edges, cvec, N = 4e5) {
  D <- diag(cvec); Dinv <- diag(1/cvec); beta <- b / cvec^2
  sig <- sapply(edges, function(e) cvec[e[1]] * cvec[e[2]])
  acc <- 0
  diagM <- matrix(rgamma(N * q, shape = alpha, rate = rep(beta, each = N)), nrow = N)
  offM  <- matrix(rnorm(N * length(edges), sd = rep(sig, each = N)), nrow = N)
  for (i in 1:N) {
    offmap <- setNames(as.list(offM[i, ]), sapply(edges, paste, collapse="-"))
    Th <- build_K(diagM[i, ], offmap, edges)
    K  <- Dinv %*% Th %*% Dinv
    if (lam_min(K) > eps) acc <- acc + det(Th)^delta
  }
  acc / N
}

cat("\n=== (3) Integral level (independent MC), D = g-prior ===\n")
cg <- c_gprior; A_det <- prod(cg^2)
Zs_full <- mc_Zstd(edges_full);  Zsc_full <- mc_Zscaled(edges_full, cg)
Zs_sub  <- mc_Zstd(edges_sub);   Zsc_sub  <- mc_Zscaled(edges_sub,  cg)
cat(sprintf("  |A|^delta = %.5f\n", A_det^delta))
cat(sprintf("  Z_scaled/Z*_std : full = %.4f , sub = %.4f   (target |A|^delta)\n",
            Zsc_full / Zs_full, Zsc_sub / Zs_sub))
cat(sprintf("  between-model ratio  scaled: %.5f   std: %.5f   (should match)\n",
            Zsc_full / Zsc_sub, Zs_full / Zs_sub))

cat("\nDone.\n")
