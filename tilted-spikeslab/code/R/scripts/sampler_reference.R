# Reference implementation of Algorithm 1 (pure R, correctness before speed).
#
# Implements exactly the manuscript's specification:
# - within-model row-block Gibbs (Proposition prop:row-block, eq:within-conditionals):
#   psi ~ Gamma(delta + n/2 + 1, beta + S_ll/2),
#   k ~ N(-Q^{-1} S_{l,N_l}, Q^{-1}), Q = sigma^{-2} I + (2 beta + S_ll) C_{N_l,N_l},
#   k_ll = psi + k' C_{N_l,N_l} k, with C = (K_{-l,-l})^{-1}.
# - between-model edge toggle (Lemmas lem:edge-curve / lem:cofactor-read,
#   Proposition prop:toggle, eq:acceptance): coefficients read from Sigma = K^{-1},
#   proposal q(x) = F(x)/N Gaussian in x (exponential diagonal), acceptance
#   A_add = [w+/w-] [pslab(0)/q(0)] J with J = 1 (joint) or Z(G-)/Z(G+)
#   (hierarchical, supplied via log_J).
#
# Convention: for the toggled pair (i, j) the second index j plays the role of
# the last-ordered node, whose diagonal k_jj varies along the curve.
# This file defines functions only; sourced by verification scripts.
# The C++ port is R/src/spikeslab_sampler.cpp (used for all studies); see
# notes/sampler-reference-implementation.md for the verification history.

sample_row_block <- function(K, G, l, S_mat, n, delta, beta, sigma) {
  q <- nrow(K)
  Nl <- which(G[l, ] == 1L)
  Sll <- S_mat[l, l]
  psi <- stats::rgamma(1L, shape = delta + n / 2 + 1, rate = beta + Sll / 2)
  if (length(Nl) == 0L) {
    K[l, ] <- 0; K[, l] <- 0; K[l, l] <- psi
    return(K)
  }
  Cfull <- solve(K[-l, -l, drop = FALSE])
  idx <- match(Nl, seq_len(q)[-l])
  Cnn <- Cfull[idx, idx, drop = FALSE]
  Cnn <- (Cnn + t(Cnn)) / 2
  Q <- diag(1 / sigma^2, length(idx)) + (2 * beta + Sll) * Cnn
  ch <- chol((Q + t(Q)) / 2)
  mu <- -backsolve(ch, forwardsolve(t(ch), S_mat[l, Nl]))
  kvec <- mu + backsolve(ch, stats::rnorm(length(idx)))
  row_l <- numeric(q)
  row_l[Nl] <- kvec
  K[l, ] <- row_l; K[, l] <- row_l
  K[l, l] <- psi + drop(crossprod(kvec, Cnn %*% kvec))
  K
}

# One Metropolis-Hastings toggle of edge (i, j); log_J = log Z(G-) - log Z(G+)
# for the hierarchical specification, 0 for the joint specification.
toggle_edge <- function(K, G, i, j, S_mat, delta, beta, sigma, p_inc, log_J = 0) {
  Sig <- solve(K)
  D <- Sig[i, i] * Sig[j, j] - Sig[i, j]^2
  c2 <- sqrt(Sig[j, j] / D)
  phi <- -Sig[i, j] / sqrt(Sig[j, j] * D)
  c1 <- K[i, j] - c2 * phi
  c3 <- K[j, j] - phi^2
  Sjj <- S_mat[j, j]
  Sij <- S_mat[i, j]
  lam <- 1 / sigma^2 + (2 * beta + Sjj) / c2^2
  lin <- -Sij + (2 * beta + Sjj) * c1 / c2^2
  mu_x <- lin / lam
  log_q0 <- 0.5 * log(lam / (2 * pi)) - 0.5 * lam * mu_x^2
  log_pslab0 <- stats::dnorm(0, 0, sigma, log = TRUE)
  log_A <- log(p_inc / (1 - p_inc)) + log_pslab0 - log_q0 + log_J
  accepted <- FALSE
  if (G[i, j] == 0L) {
    if (log(stats::runif(1L)) < log_A) {
      xs <- mu_x + stats::rnorm(1L) / sqrt(lam)
      G[i, j] <- 1L; G[j, i] <- 1L
      K[i, j] <- xs; K[j, i] <- xs
      K[j, j] <- c3 + ((xs - c1) / c2)^2
      accepted <- TRUE
    }
  } else {
    if (log(stats::runif(1L)) < -log_A) {
      G[i, j] <- 0L; G[j, i] <- 0L
      K[i, j] <- 0; K[j, i] <- 0
      K[j, j] <- c3 + (c1 / c2)^2
      accepted <- TRUE
    }
  }
  list(K = K, G = G, accepted = accepted)
}

# Full chain: each iteration is one row-block sweep plus one toggle attempt per
# candidate pair in random order. spec = "joint" or "hierarchical"; for the
# hierarchical specification log_z must map a graph mask to log Z(Gamma)
# (a numeric vector indexed by mask + 1).
run_chain <- function(q, S_mat, n, delta, beta, sigma, p_inc, spec,
                      n_iter, n_burn, seed, log_z = NULL, G_init = NULL) {
  set.seed(seed)
  edges <- t(utils::combn(q, 2))
  M <- nrow(edges)
  G <- if (is.null(G_init)) matrix(0L, q, q) else G_init
  K <- diag((delta + 1) / beta, q)
  if (any(G == 1L)) K[G == 1L] <- 0.01
  mask_of <- function(G) {
    bits <- vapply(seq_len(M), function(k) G[edges[k, 1], edges[k, 2]], integer(1))
    sum(bits * 2^(seq_len(M) - 1L))
  }
  visits <- integer(2^M)
  n_att <- 0L; n_acc <- 0L
  for (it in seq_len(n_iter)) {
    for (l in seq_len(q)) {
      K <- sample_row_block(K, G, l, S_mat, n, delta, beta, sigma)
    }
    for (k in sample.int(M)) {
      i <- edges[k, 1]; j <- edges[k, 2]
      log_J <- 0
      if (spec == "hierarchical") {
        m0 <- mask_of(G)
        bit <- 2^(k - 1L)
        m_without <- if (G[i, j] == 1L) m0 - bit else m0
        m_with <- m_without + bit
        log_J <- log_z[m_without + 1L] - log_z[m_with + 1L]
      }
      res <- toggle_edge(K, G, i, j, S_mat, delta, beta, sigma, p_inc, log_J)
      K <- res$K; G <- res$G
      n_att <- n_att + 1L
      if (res$accepted) n_acc <- n_acc + 1L
    }
    if (it > n_burn) {
      m <- mask_of(G)
      visits[m + 1L] <- visits[m + 1L] + 1L
    }
  }
  list(visits = visits, acc_rate = n_acc / n_att, edges = edges)
}
