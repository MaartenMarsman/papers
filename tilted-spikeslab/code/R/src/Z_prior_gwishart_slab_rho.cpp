// R/src/Z_prior_gwishart_slab_rho.cpp
//
// Variant of Z_prior_gwishart_slab_k12.cpp that ALSO records the two
// diagonals K[0,0] and K[1,1] per accepted joint (Gamma, K) draw, so that
// the induced partial correlation
//
//     rho_{1,2} = -K[0,1] / sqrt(K[0,0] * K[1,1])
//
// can be estimated on gamma_{1,2} = 1 draws. Everything else (the alpha = 1
// G-Wishart-slab perfect sampler, the condition_on_edge_01 trick, the
// partial residuals a, c and centre m) is identical to the k12 recorder.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>

static inline void akm_propose_rho(int q,
                                    double delta_W,
                                    const arma::imat& E,
                                    const arma::ivec& nu,
                                    arma::mat& Phi_out) {
  Phi_out.zeros();
  for (int i = 0; i < q; ++i) {
    Phi_out(i, i) = std::sqrt(R::rchisq(delta_W + nu[i]));
    for (int j = i + 1; j < q; ++j) {
      if (E(i, j) == 1) Phi_out(i, j) = R::rnorm(0.0, 1.0);
    }
  }
  for (int i = 0; i < q; ++i) {
    for (int j = i + 1; j < q; ++j) {
      if (E(i, j) == 0) {
        double s = 0.0;
        for (int v = 0; v < i; ++v) s += Phi_out(v, i) * Phi_out(v, j);
        Phi_out(i, j) = -s / Phi_out(i, i);
      }
    }
  }
}

// [[Rcpp::export]]
Rcpp::List gwishart_slab_record_rho_cpp(int q,
                                        double delta, double beta, double sigma,
                                        double gamma_prob,
                                        int n_target, int n_proposal_cap,
                                        int max_inner_tries, int seed,
                                        bool condition_on_edge_01 = false) {
  Rcpp::Environment base = Rcpp::Environment::base_env();
  Rcpp::Function set_seed = base["set.seed"];
  set_seed(seed);

  const double delta_W = 2.0 * delta + 2.0;
  const double inv_2beta = 1.0 / (2.0 * beta);
  const double inv_2sigma2 = 1.0 / (2.0 * sigma * sigma);

  arma::imat E(q, q, arma::fill::zeros);
  arma::ivec nu(q, arma::fill::zeros);
  arma::mat Phi(q, q, arma::fill::zeros);
  arma::mat K(q, q);

  std::vector<double> k12, k11, k22;
  std::vector<int>    gamma12;
  std::vector<double> a_vec, c_vec, m_vec;
  k12.reserve(2048); k11.reserve(2048); k22.reserve(2048);
  gamma12.reserve(2048);
  a_vec.reserve(2048); c_vec.reserve(2048); m_vec.reserve(2048);

  // Indices J = {2, ..., q-1}
  arma::uvec rest;
  if (q >= 3) {
    rest.set_size(q - 2);
    for (int k = 2; k < q; ++k) rest(k - 2) = k;
  }

  int n_accepted = 0;
  int n_proposals = 0;

  // Bookkeeping added after the 2026-10-08 code review: graphs proposed and
  // kept, and failures of the rest-block inversion. NOTE: a graph is dropped
  // when no K is accepted within max_inner_tries proposals, which favors
  // graphs with a high slab acceptance rate when that rate is below about
  // 1/max_inner_tries (observed at q >= 20 with delta >= 1); compare
  // n_graphs_kept with n_graphs_tried before using such cells.
  long long n_graphs_tried = 0, n_graphs_kept = 0, n_inv_fail = 0;
  while (n_accepted < n_target && n_proposals < n_proposal_cap) {
    ++n_graphs_tried;
    E.zeros(); nu.zeros();
    if (condition_on_edge_01) {
      E(0, 1) = 1; ++nu[0];
    }
    for (int i = 0; i < q - 1; ++i) {
      for (int j = i + 1; j < q; ++j) {
        if (condition_on_edge_01 && i == 0 && j == 1) continue;
        if (R::unif_rand() < gamma_prob) {
          E(i, j) = 1; ++nu[i];
        }
      }
    }
    for (int inner = 0; inner < max_inner_tries; ++inner) {
      if (n_proposals >= n_proposal_cap) break;
      ++n_proposals;
      akm_propose_rho(q, delta_W, E, nu, Phi);
      K = Phi.t() * Phi;
      K *= inv_2beta;
      K = arma::symmatu(K);
      double log_w = 0.0;
      for (int i = 0; i < q - 1; ++i) {
        for (int j = i + 1; j < q; ++j) {
          if (E(i, j) == 1) log_w -= K(i, j) * K(i, j) * inv_2sigma2;
        }
      }
      if (std::log(R::unif_rand()) >= log_w) continue;
      ++n_accepted;
      ++n_graphs_kept;

      // The two diagonals and the off-diagonal (the pieces of rho_{1,2}).
      k11.push_back(K(0, 0));
      k22.push_back(K(1, 1));
      k12.push_back(K(0, 1));
      gamma12.push_back(E(0, 1));

      // Partial residuals a, c and partial-covariance centre m (q >= 3).
      double a_val = K(0, 0), c_val = K(1, 1), m_val = 0.0;
      if (q >= 3) {
        arma::vec b0 = K.submat(arma::uvec({0}), rest).t();
        arma::vec b1 = K.submat(arma::uvec({1}), rest).t();
        arma::mat Krest = K.submat(rest, rest);
        arma::mat Krest_inv;
        bool ok = arma::inv_sympd(Krest_inv, Krest);
        if (!ok) ok = arma::inv(Krest_inv, Krest);   // pivoted LU fallback
        if (!ok) ++n_inv_fail;
        if (ok) {
          a_val = K(0, 0) - arma::as_scalar(b0.t() * Krest_inv * b0);
          c_val = K(1, 1) - arma::as_scalar(b1.t() * Krest_inv * b1);
          m_val = arma::as_scalar(b0.t() * Krest_inv * b1);
        } else {
          a_val = arma::datum::nan;
          c_val = arma::datum::nan;
          m_val = arma::datum::nan;
        }
      }
      a_vec.push_back(a_val);
      c_vec.push_back(c_val);
      m_vec.push_back(m_val);

      break;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("k12") = arma::vec(k12),
    Rcpp::Named("k11") = arma::vec(k11),
    Rcpp::Named("k22") = arma::vec(k22),
    Rcpp::Named("gamma12") = arma::ivec(gamma12),
    Rcpp::Named("a") = arma::vec(a_vec),
    Rcpp::Named("c") = arma::vec(c_vec),
    Rcpp::Named("m") = arma::vec(m_vec),
    Rcpp::Named("n_proposals") = n_proposals,
    Rcpp::Named("n_accepted") = n_accepted,
    Rcpp::Named("n_graphs_tried") = n_graphs_tried,
    Rcpp::Named("n_graphs_kept") = n_graphs_kept,
    Rcpp::Named("n_inv_fail") = n_inv_fail,
    Rcpp::Named("q") = q,
    Rcpp::Named("delta") = delta
  );
}
