// R/src/Z_prior_brute.cpp
//
// Brute-force MC estimator of the modified spike-and-slab prior normaliser
// Z_prior(Gamma; delta, epsilon, alpha, beta, sigma) under either variant.
// Variant is determined by the off-diagonal SD matrix slab_sd:
//   - Discrete : slab_sd(i,j) = 0 for (i,j) not in E(Gamma), sigma for (i,j) in E.
//   - Continuous: slab_sd(i,j) = v0 for (i,j) not in E(Gamma), v1 for (i,j) in E.
//
// Base measure (the importance-sampling proposal):
//   K_ll ~ Gamma(alpha, rate = beta),
//   K_ij ~ N(0, slab_sd(i,j)^2) for i < j (symmetrised).
//
// Importance weight:
//   w(K) = |K|^delta * 1{lambda_min(K) > epsilon}.
//
// MC estimate of Z_prior is mean of w(K) over the base draws.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List z_prior_brute_cpp(const arma::imat& G,
                             double alpha, double beta,
                             const arma::mat& slab_sd,
                             double delta, double epsilon,
                             int n_draws, int seed) {
  const int q = G.n_rows;
  arma::mat K(q, q, arma::fill::zeros);
  arma::vec eigenvalues(q);

  Rcpp::Environment base = Rcpp::Environment::base_env();
  Rcpp::Function set_seed = base["set.seed"];
  set_seed(seed);

  double sum_w  = 0.0;
  double sum_w2 = 0.0;
  long long n_accept = 0;
  const double inv_beta = 1.0 / beta;

  for (int t = 0; t < n_draws; ++t) {
    K.zeros();
    for (int l = 0; l < q; ++l) {
      K(l, l) = R::rgamma(alpha, inv_beta);
    }
    for (int i = 0; i < q - 1; ++i) {
      for (int j = i + 1; j < q; ++j) {
        const double sd = slab_sd(i, j);
        const double v  = (sd > 0.0) ? R::rnorm(0.0, sd) : 0.0;
        K(i, j) = v;
        K(j, i) = v;
      }
    }

    if (!arma::eig_sym(eigenvalues, K)) continue;

    const double lmin = eigenvalues.min();
    if (lmin > epsilon) {
      const double log_det = arma::sum(arma::log(eigenvalues));
      const double w = std::exp(delta * log_det);
      sum_w  += w;
      sum_w2 += w * w;
      ++n_accept;
    }
  }

  const double z_hat = sum_w / n_draws;
  const double var_w = (sum_w2 - sum_w * z_hat) / static_cast<double>(n_draws - 1);
  const double z_se  = std::sqrt(var_w / n_draws);

  return Rcpp::List::create(
    Rcpp::Named("z_hat")     = z_hat,
    Rcpp::Named("z_se")      = z_se,
    Rcpp::Named("log_z_hat") = (z_hat > 0.0) ? std::log(z_hat)
                                              : -std::numeric_limits<double>::infinity(),
    Rcpp::Named("n_accept")  = static_cast<double>(n_accept),
    Rcpp::Named("n_draws")   = n_draws
  );
}
