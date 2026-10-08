// R/src/Z_prior_tilt_resample.cpp
//
// Rejection sampler from the delta-tilted spike-and-slab prior, supported on
// {K positive-definite}. Proposes K with diagonals Gamma(alpha + delta, beta)
// and off-diagonals from the spike/slab. Accepts with probability
//     r(K) = (det K / prod K_ll)^delta * 1{K positive-definite},
// which is in [0, 1] by Hadamard's inequality on the determinant of a PD
// matrix. Accepted draws are exact i.i.d. samples from the delta-tilted prior
// conditional on K positive-definite. No importance weights; each accepted
// draw counts as one effective sample.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>

// [[Rcpp::export]]
Rcpp::List z_prior_tilt_resample_cpp(int q,
                                     double alpha, double beta, double delta,
                                     int variant,           // 0 = discrete, 1 = continuous
                                     double sigma,
                                     double v0, double v1,
                                     double gamma_prob,
                                     int n_proposals, int seed) {
  arma::mat K(q, q, arma::fill::zeros);
  arma::vec eigenvalues(q);

  Rcpp::Environment base = Rcpp::Environment::base_env();
  Rcpp::Function set_seed = base["set.seed"];
  set_seed(seed);

  const double shape = alpha + delta;
  const double inv_beta = 1.0 / beta;

  std::vector<double> lmin_acc;
  std::vector<double> log_det_acc;
  std::vector<double> diag_mean_acc;
  lmin_acc.reserve(1024);
  log_det_acc.reserve(1024);
  diag_mean_acc.reserve(1024);

  int n_pd = 0;
  int n_accepted = 0;

  for (int t = 0; t < n_proposals; ++t) {
    K.zeros();
    double trace_K = 0.0;
    double sum_log_diag = 0.0;
    for (int l = 0; l < q; ++l) {
      const double k_ll = R::rgamma(shape, inv_beta);
      K(l, l) = k_ll;
      trace_K += k_ll;
      sum_log_diag += std::log(k_ll);
    }
    for (int i = 0; i < q - 1; ++i) {
      for (int j = i + 1; j < q; ++j) {
        const bool gamma_ij = (R::unif_rand() < gamma_prob);
        double v = 0.0;
        if (variant == 0) {
          v = gamma_ij ? R::rnorm(0.0, sigma) : 0.0;
        } else {
          v = R::rnorm(0.0, gamma_ij ? v1 : v0);
        }
        K(i, j) = v;
        K(j, i) = v;
      }
    }
    if (!arma::eig_sym(eigenvalues, K)) continue;
    const double lmin = eigenvalues.min();
    if (lmin <= 0.0) continue;
    ++n_pd;
    const double log_det = arma::sum(arma::log(eigenvalues));
    // Acceptance ratio r = (det K / prod K_ll)^delta in [0, 1] for PD K.
    if (delta > 0.0) {
      const double log_r = delta * (log_det - sum_log_diag);
      if (std::log(R::unif_rand()) >= log_r) continue;
    }
    ++n_accepted;
    lmin_acc.push_back(lmin);
    log_det_acc.push_back(log_det);
    diag_mean_acc.push_back(trace_K / static_cast<double>(q));
  }

  return Rcpp::List::create(
    Rcpp::Named("lmin")        = arma::vec(lmin_acc),
    Rcpp::Named("log_det")     = arma::vec(log_det_acc),
    Rcpp::Named("diag_mean")   = arma::vec(diag_mean_acc),
    Rcpp::Named("n_proposals") = n_proposals,
    Rcpp::Named("n_pd")        = n_pd,
    Rcpp::Named("n_accepted")  = n_accepted,
    Rcpp::Named("q")           = q,
    Rcpp::Named("delta")       = delta
  );
}
