// R/src/Z_prior_gwishart_slab.cpp
//
// Perfect (rejection) sampler for the determinant-tilted spike-and-slab
// prior in the alpha = 1 (exponential diagonal) regime. The prior at alpha
// = 1, conditional on graph Gamma, factorises as
//
//   p(K | Gamma) propto  |K|^delta * exp(-beta * tr(K))
//                       * prod_{(i,j) in E(Gamma)} exp(-K_{ij}^2 / (2 sigma^2))
//                       * 1{K in M^+_{G(Gamma)}}
//
// which is a G-Wishart W_G(delta_W = 2*delta + 2, D = 2*beta*I) envelope
// times the bounded edge-slab factor w(K) in (0, 1]. The G-Wishart is sampled
// exactly via Atay-Kayis-Massam-Roverato Cholesky completion; we then accept
// with probability w(K). Accepted draws are i.i.d. from p(K | Gamma).
//
// For the joint sampler we draw Gamma_{ij} iid Bernoulli(gamma_prob) and use
// nested rejection: keep proposing K under the same Gamma until acceptance
// (capped). The outputs cover both the joint (Gamma, K) draws and per-Gamma
// summaries.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>

// Single AKM proposal of K given a fixed graph adjacency matrix and degree
// vector nu. Returns the upper-triangular Cholesky factor Phi (zero below
// the diagonal) such that K = Phi^T Phi / (2 beta). Phi is written into Phi_out.
static inline void akm_propose(int q,
                                double delta_W,
                                const arma::imat& E,
                                const arma::ivec& nu,
                                arma::mat& Phi_out) {
  Phi_out.zeros();
  // Free entries
  for (int i = 0; i < q; ++i) {
    Phi_out(i, i) = std::sqrt(R::rchisq(delta_W + nu[i]));
    for (int j = i + 1; j < q; ++j) {
      if (E(i, j) == 1) {
        Phi_out(i, j) = R::rnorm(0.0, 1.0);
      }
    }
  }
  // Non-free entries: completion theta_{ij} = 0 gives
  //   Phi_{ij} = - sum_{v < i} Phi_{vi} Phi_{vj} / Phi_{ii}
  for (int i = 0; i < q; ++i) {
    for (int j = i + 1; j < q; ++j) {
      if (E(i, j) == 0) {
        double s = 0.0;
        for (int v = 0; v < i; ++v) {
          s += Phi_out(v, i) * Phi_out(v, j);
        }
        Phi_out(i, j) = -s / Phi_out(i, i);
      }
    }
  }
}

// Perfect sampler for the tilted spike-and-slab prior at alpha = 1, joint
// over (Gamma, K). Gamma_{ij} iid Bernoulli(gamma_prob). Per Gamma we run a
// nested rejection on the slab factor; we cap the inner loop at
// max_inner_tries to bound wall time. The outer loop runs until either
// n_target accepted joint draws are obtained, or the total number of K
// proposals reaches n_proposal_cap.

// [[Rcpp::export]]
Rcpp::List gwishart_slab_resample_cpp(int q,
                                      double delta,
                                      double beta,
                                      double sigma,
                                      double gamma_prob,
                                      int n_target,
                                      int n_proposal_cap,
                                      int max_inner_tries,
                                      int seed) {
  Rcpp::Environment base = Rcpp::Environment::base_env();
  Rcpp::Function set_seed = base["set.seed"];
  set_seed(seed);

  const double delta_W = 2.0 * delta + 2.0;
  const double inv_2beta = 1.0 / (2.0 * beta);
  const double inv_2sigma2 = 1.0 / (2.0 * sigma * sigma);

  arma::imat E(q, q, arma::fill::zeros);
  arma::ivec nu(q);
  arma::mat Phi(q, q, arma::fill::zeros);
  arma::mat K(q, q);
  arma::vec eigvals;

  std::vector<double> lmin_acc;
  std::vector<double> log_det_acc;
  std::vector<int>    n_edges_acc;
  std::vector<int>    inner_tries_acc;
  lmin_acc.reserve(2048);
  log_det_acc.reserve(2048);
  n_edges_acc.reserve(2048);
  inner_tries_acc.reserve(2048);

  int n_accepted     = 0;
  int n_proposals    = 0;
  int n_graphs_tried = 0;
  int n_graphs_kept  = 0;

  while (n_accepted < n_target && n_proposals < n_proposal_cap) {
    // Sample a new graph
    E.zeros(); nu.zeros();
    int n_edges = 0;
    for (int i = 0; i < q - 1; ++i) {
      for (int j = i + 1; j < q; ++j) {
        if (R::unif_rand() < gamma_prob) {
          E(i, j) = 1; ++nu[i]; ++n_edges;
        }
      }
    }
    ++n_graphs_tried;

    bool kept = false;
    for (int inner = 0; inner < max_inner_tries; ++inner) {
      if (n_proposals >= n_proposal_cap) break;
      ++n_proposals;

      akm_propose(q, delta_W, E, nu, Phi);

      K = Phi.t() * Phi;
      K *= inv_2beta;
      K = arma::symmatu(K);  // enforce exact symmetry against floating-point drift

      // Slab log-weight on edges of E
      double log_w = 0.0;
      for (int i = 0; i < q - 1; ++i) {
        for (int j = i + 1; j < q; ++j) {
          if (E(i, j) == 1) {
            log_w -= K(i, j) * K(i, j) * inv_2sigma2;
          }
        }
      }
      if (std::log(R::unif_rand()) >= log_w) continue;

      // Accept this (Gamma, K) draw
      bool ok = arma::eig_sym(eigvals, K);
      if (!ok) continue;
      const double lmin = eigvals.min();
      const double log_det = arma::sum(arma::log(eigvals));

      ++n_accepted;
      lmin_acc.push_back(lmin);
      log_det_acc.push_back(log_det);
      n_edges_acc.push_back(n_edges);
      inner_tries_acc.push_back(inner + 1);
      kept = true;
      break;
    }
    if (kept) ++n_graphs_kept;
  }

  return Rcpp::List::create(
    Rcpp::Named("lmin")           = arma::vec(lmin_acc),
    Rcpp::Named("log_det")        = arma::vec(log_det_acc),
    Rcpp::Named("n_edges")        = arma::ivec(n_edges_acc),
    Rcpp::Named("inner_tries")    = arma::ivec(inner_tries_acc),
    Rcpp::Named("n_target")       = n_target,
    Rcpp::Named("n_proposals")    = n_proposals,
    Rcpp::Named("n_accepted")     = n_accepted,
    Rcpp::Named("n_graphs_tried") = n_graphs_tried,
    Rcpp::Named("n_graphs_kept")  = n_graphs_kept,
    Rcpp::Named("q")              = q,
    Rcpp::Named("delta")          = delta
  );
}

// Diagnostic: sample only from the G-Wishart envelope (no slab rejection)
// for a fixed graph. Used for empty-graph sanity checks.
// [[Rcpp::export]]
Rcpp::List gwishart_envelope_only_cpp(int q,
                                       double delta,
                                       double beta,
                                       Rcpp::IntegerMatrix E_in,
                                       int n_draws, int seed) {
  Rcpp::Environment base = Rcpp::Environment::base_env();
  Rcpp::Function set_seed = base["set.seed"];
  set_seed(seed);

  const double delta_W = 2.0 * delta + 2.0;
  const double inv_2beta = 1.0 / (2.0 * beta);

  arma::imat E(q, q, arma::fill::zeros);
  arma::ivec nu(q, arma::fill::zeros);
  for (int i = 0; i < q - 1; ++i) {
    for (int j = i + 1; j < q; ++j) {
      E(i, j) = E_in(i, j);
      if (E(i, j) == 1) ++nu[i];
    }
  }
  arma::mat Phi(q, q, arma::fill::zeros);
  arma::mat K(q, q);
  arma::vec eigvals;

  arma::vec lmin(n_draws);
  arma::vec log_det(n_draws);
  arma::mat diag_vals(n_draws, q);

  for (int t = 0; t < n_draws; ++t) {
    akm_propose(q, delta_W, E, nu, Phi);
    K = Phi.t() * Phi;
    K *= inv_2beta;
    K = arma::symmatu(K);
    arma::eig_sym(eigvals, K);
    lmin(t) = eigvals.min();
    log_det(t) = arma::sum(arma::log(eigvals));
    for (int l = 0; l < q; ++l) diag_vals(t, l) = K(l, l);
  }

  return Rcpp::List::create(
    Rcpp::Named("lmin")      = lmin,
    Rcpp::Named("log_det")   = log_det,
    Rcpp::Named("diag_vals") = diag_vals,
    Rcpp::Named("q")         = q,
    Rcpp::Named("delta")     = delta
  );
}
