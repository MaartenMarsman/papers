// R/src/monotonicity_check.cpp
//
// Monte Carlo estimator of the tilted spike-and-slab prior normalizing
// constant Z(Gamma; delta) at alpha = 1 (exponential diagonal), with
// DRAW ALIGNMENT across graphs: every entry of the symmetric matrix is
// drawn on every iteration (diagonal Exp(beta), off-diagonal N(0, sigma^2)),
// in a fixed order, and non-edges are zeroed afterwards. Calling this
// function for two graphs with the same seed therefore uses common random
// numbers, so ratios of the returned estimates have strongly reduced
// Monte Carlo error.
//
// Convention: normalized entry densities inside f (draws come from those
// densities), so Z-hat = mean of |K|^delta * 1{K positive definite}.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List z_graph_aligned_cpp(const arma::imat& G,
                               double beta, double sigma, double delta,
                               int n_draws, int seed) {
  const int q = G.n_rows;
  arma::mat K(q, q, arma::fill::zeros);
  arma::vec eigenvalues(q);

  Rcpp::Environment base = Rcpp::Environment::base_env();
  Rcpp::Function set_seed = base["set.seed"];
  set_seed(seed);

  const double scale = 1.0 / beta;  // Exp(beta) has scale 1/beta
  double sum_w = 0.0, sum_w2 = 0.0;
  long long n_pd = 0;

  for (int t = 0; t < n_draws; ++t) {
    // Fixed draw order regardless of the graph: diagonal, then upper triangle.
    for (int l = 0; l < q; ++l) {
      K(l, l) = R::rexp(scale);
    }
    for (int i = 0; i < q - 1; ++i) {
      for (int j = i + 1; j < q; ++j) {
        const double v = R::rnorm(0.0, sigma);
        const double kept = (G(i, j) == 1) ? v : 0.0;
        K(i, j) = kept;
        K(j, i) = kept;
      }
    }

    if (!arma::eig_sym(eigenvalues, K)) continue;
    if (eigenvalues.min() > 0.0) {
      const double w = std::exp(delta * arma::sum(arma::log(eigenvalues)));
      sum_w  += w;
      sum_w2 += w * w;
      ++n_pd;
    }
  }

  const double z_hat = sum_w / n_draws;
  const double var_w = (sum_w2 - sum_w * z_hat) / static_cast<double>(n_draws - 1);

  return Rcpp::List::create(
    Rcpp::Named("z_hat") = z_hat,
    Rcpp::Named("z_se")  = std::sqrt(var_w / n_draws),
    Rcpp::Named("n_pd")  = static_cast<double>(n_pd),
    Rcpp::Named("n_draws") = n_draws
  );
}
