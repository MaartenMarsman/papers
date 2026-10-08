// R/src/Z_prior_diag.cpp
//
// Monte Carlo sampler from the UNRESTRICTED entry-wise construction (Gamma
// diagonal, spike-and-slab off-diagonals, no positive-definiteness
// restriction). Records lambda_min(K), tr(K)/q and whether the draw is
// positive definite, which gives the loss of prior mass and the boundary
// concentration of Section sec:collapse.
//
// Used by R/scripts/cone_interior_diagnostics_expanded.R (figure
// fig:cone-interior) and by the superseded cone_interior_diagnostics.R.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Diagnostic that samples (K, gamma) jointly from the base measure with
// gamma_ij ~ Bernoulli(gamma_prob) per draw (matching the Z project's
// joint-sampling diagnostic of unconditional Pr[K PD]).
//
// variant = 0 (discrete): K_ij = 0 if gamma_ij = 0; ~N(0, sigma^2) otherwise.
// variant = 1 (continuous): K_ij ~ N(0, v0^2) if gamma_ij = 0; ~N(0, v1^2) otherwise.
//
// [[Rcpp::export]]
Rcpp::List z_prior_diag_cpp(int q,
                            double alpha, double beta,
                            int variant,           // 0 = discrete, 1 = continuous
                            double sigma,          // discrete slab sd
                            double v0, double v1,  // continuous spike/slab sd
                            double gamma_prob,
                            int n_draws, int seed) {
  arma::mat K(q, q, arma::fill::zeros);
  arma::vec eigenvalues(q);

  Rcpp::Environment base = Rcpp::Environment::base_env();
  Rcpp::Function set_seed = base["set.seed"];
  set_seed(seed);

  arma::vec lmin_all(n_draws);
  arma::vec log_det_all(n_draws);
  arma::vec diag_mean_all(n_draws);
  arma::ivec pd_flag(n_draws);

  const double inv_beta = 1.0 / beta;

  for (int t = 0; t < n_draws; ++t) {
    K.zeros();
    double trace_K = 0.0;
    for (int l = 0; l < q; ++l) {
      const double k_ll = R::rgamma(alpha, inv_beta);
      K(l, l) = k_ll;
      trace_K += k_ll;
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

    diag_mean_all(t) = trace_K / static_cast<double>(q);

    if (!arma::eig_sym(eigenvalues, K)) {
      lmin_all(t)    = arma::datum::nan;
      log_det_all(t) = arma::datum::nan;
      pd_flag(t) = 0;
      continue;
    }
    const double lmin = eigenvalues.min();
    lmin_all(t) = lmin;
    pd_flag(t)  = (lmin > 0.0) ? 1 : 0;
    log_det_all(t) = (lmin > 0.0) ? arma::sum(arma::log(eigenvalues))
                                  : arma::datum::nan;
  }

  return Rcpp::List::create(
    Rcpp::Named("lmin")      = lmin_all,
    Rcpp::Named("log_det")   = log_det_all,
    Rcpp::Named("diag_mean") = diag_mean_all,
    Rcpp::Named("pd")        = pd_flag,
    Rcpp::Named("n_draws")   = n_draws,
    Rcpp::Named("q")         = q
  );
}
