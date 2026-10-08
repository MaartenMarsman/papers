// C++ implementation of Algorithm alg:sampler of the manuscript (the complete
// sampler), matching the pure-R reference in R/scripts/sampler_reference.R,
// which was verified against exact enumeration (see the verification
// appendix, checks V4 to V8). Compile with Rcpp::sourceCpp().
//
// Model: entry-wise discrete spike-and-slab prior on the precision matrix K
// with slab N(0, sigma^2), exponential diagonal with rate beta (alpha = 1 only;
// other Gamma shapes are not implemented), determinant tilt |K|^delta, and an
// independent Bernoulli(p_inc) prior on edges.
//
// Cost structure (Section sec:algorithm):
// - Sigma = K^{-1} is maintained along the chain;
//   C_{N,N} = Sigma_{N,N} - Sigma_{N,l} Sigma_{l,N} / Sigma_{ll} in O(d^2);
//   the Gaussian draw factorizes the d x d precision Q at O(d^3);
//   a symmetric rank-two Sherman-Morrison-Woodbury update refreshes Sigma in
//   O(q^2) after each row update and each accepted toggle.
// - Candidate toggles read three entries of Sigma in O(1); the proposal
//   q(x) = N(mu_x, 1/lam) is the closed Gaussian form of F(x)/N at the
//   exponential diagonal; acceptance per eq:acceptance.
// - Sigma is recomputed exactly (inv_sympd, with inv as fallback) in three
//   situations: every refresh_every sweeps (refresh_every = 0 disables this,
//   and then no drift is measured), when the Cholesky factorization in a row
//   draw fails, and when a diagonal slack psi falls below 1e-6. The reported
//   max_sigma_drift is the largest max|Sigma - K^{-1}| seen at the periodic
//   refreshes.
//
// sampler_chain_cpp(S, n_obs, q, delta, beta, sigma, p_inc, spec, log_z,
//                   n_iter, n_burn, refresh_every, track_masks,
//                   n_aux = 0, stat_thin = 0, time_blocks = false)
//   S            q x q scatter matrix Y'Y (not a covariance); the zero matrix
//                for a prior-level run.
//   n_obs        number of observations n; 0 runs the chain at the prior level.
//   q            dimension (must equal nrow(S)).
//   delta        tilt exponent.
//   beta         rate of the exponential diagonal prior. In the manuscript's
//                standardized parametrization eta = beta * sigma; the scripts
//                pass sigma = 1 and beta = eta.
//   sigma        standard deviation of the slab.
//   p_inc        edge-inclusion probability of the Bernoulli graph prior.
//   spec         0 = joint specification (J = 1);
//                1 = hierarchical specification with exact per-graph constants
//                    supplied in log_z (a vector of length 2^M, M = q(q-1)/2,
//                    indexed by graph mask + 1; feasible for q <= 5 only);
//                2 = hierarchical specification with the coupled Monte Carlo
//                    estimator of Z(G-)/Z(G+), n_aux draws per arm and per
//                    candidate (the plug-in chain of Section sec:hier-ratio).
//   log_z        the table for spec = 1; pass numeric(0) otherwise.
//   n_iter       iterations; each iteration is one row-block sweep over all
//                nodes followed by one toggle attempt for every pair in a
//                fresh random order.
//   n_burn       burn-in iterations excluded from edge_freq, m_hist, E_m and
//                the thinned statistics (acc_rate includes burn-in).
//   refresh_every exact recomputation of Sigma every this many sweeps.
//   track_masks  record visit counts per graph (requires M <= 62; memory 2^M).
//   n_aux        spec = 2: auxiliary draws per arm per candidate (required);
//                spec = 1: if > 0, additionally runs the estimator at every
//                candidate and records how often the accept/reject decision
//                would differ (disc_rate), which also advances the RNG.
//   stat_thin    if > 0, record (m, k_11, k_12, log|K|) every stat_thin-th
//                post-burn-in iteration.
//   time_blocks  if TRUE, return block timers in seconds.
//   Initial state: the empty graph with K = ((delta + 1)/beta) I.
//   Graph masks: edge k (combn(q, 2) order) is bit k of the mask.
//
// Returns a list:
//   edge_freq    q x q, post-burn-in inclusion frequency of each edge in the
//                UPPER triangle (lower triangle and diagonal are zero); these
//                are the posterior edge-inclusion probabilities.
//   m_hist       counts of the edge count m = 0..M over post-burn-in iterations.
//   E_m          mean post-burn-in edge count.
//   acc_rate     accepted / attempted toggles, including burn-in.
//   max_sigma_drift  see above.
//   mask_visits  (track_masks) visit counts per graph mask, 0-based index.
//   rec_m, rec_k11, rec_k12, rec_logdet  (stat_thin > 0) thinned statistics.
//   timing       (time_blocks) sweep_s, toggle_s, accept_s; accept_s is the
//                time of the rank-two updates after accepted toggles and is
//                nested inside toggle_s.
//   disc_rate    (spec = 1 with n_aux > 0) decision-discrepancy rate.
//   The chain does not return draws of K.
//
// akm_ratio_cpp(G_minus, i, j, delta, beta, sigma, N, coupled = TRUE,
//               return_draws = FALSE)
//   Coupled Monte Carlo estimate of Z(G-)/Z(G+) for adding the edge (i, j)
//   (1-based, i < j) to the graph G_minus (0/1 integer matrix; only the upper
//   triangle is read), with N draws per arm; coupled = FALSE draws the two arms
//   independently. Returns log_J = log Z(G-)/Z(G+), the two arm means (on the
//   linear scale, which may underflow to 0 while log_J stays finite) and, with
//   return_draws, the per-draw LOG weights (-Inf for underflowing draws).
//
// RNG: R's RNG (seed with set.seed() before calling).

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <chrono>

static const double LOG2PI = std::log(2.0 * M_PI);

// Symmetric rank-two downdate of Sigma after K' = K + U W U^t, in the
// Woodbury form that needs no inverse of W (valid for singular W, smooth as
// entries of W approach zero):
//   Sigma' = Sigma - (Sigma U) W (I + U^t Sigma U W)^{-1} (U^t Sigma).
// K' positive definite guarantees the 2x2 system is invertible.
static void smw_rank2(arma::mat& Sigma, const arma::mat& U, const arma::mat& W) {
  arma::mat SU = Sigma * U;                                  // q x 2
  arma::mat M2 = arma::eye(2, 2) + U.t() * SU * W;           // 2 x 2
  arma::mat Sol = arma::solve(M2, SU.t());                   // 2 x q
  Sigma -= SU * W * Sol;
  Sigma = 0.5 * (Sigma + Sigma.t());
}

// ---------------------------------------------------------------------------
// Coupled Monte Carlo estimator of J = Z(G-)/Z(G+) (Section 4.4).
//
// AKM-proposal representation at the exponential diagonal: with
// b = 2 delta + 2 and nu_l = #{k > l : (l, k) in E} under the natural
// ordering, the free Cholesky elements of Psi = sqrt(2 beta) Phi are
// independent under q_G: Psi_ll^2 ~ chi-squared(b + nu_l) and Psi_ls ~ N(0,1)
// for free (l, s); non-free entries complete by the zero constraints
// Psi_rs = -(1/Psi_rr) sum_{t<r} Psi_tr Psi_ts. Then
//   Z(G) = (2 pi sigma^2)^{-m/2} C_G E_q[ w{K(Psi)} f_G(Psi) ]
//   (up to the graph-independent factor beta^q of eq:Z-proposal, which cancels
//   in every ratio and is omitted here),
//   C_G  = (pi/beta)^{m/2} prod_l beta^{-(b+nu_l)/2} Gamma((b+nu_l)/2),
// with f_G = exp(-0.5 sum_nonfree Psi_rs^2), K = Psi^t Psi / (2 beta),
// w(K) = prod_edges exp(-k_ij^2 / (2 sigma^2)). For a single-edge difference
// at edge (i, j), i < j:
//   J = (2 pi sigma^2)^{1/2} (beta/sqrt(pi)) Gamma(a)/Gamma(a + 1/2) R,
//   a = delta + 1 + nu_i(G-)/2,  R = E_{q_-}[w f] / E_{q_+}[w f].
// Coupling (the method's definition, eq:ratio-estimator): shared normals for
// common free off-diagonals; the toggled edge's free normal under G+ is
// independent; the one diagonal whose degrees of freedom differ by one shares
// the lower-degree chi-squared variate, and G+ adds the square of one
// additional independent standard normal.

// log of w{K(Psi)} f_G(Psi): completes the non-free entries of Psi in place.
static double akm_log_wf(arma::mat& Psi, const arma::imat& G,
                         double beta, double sig2) {
  const int q = Psi.n_rows;
  double log_f = 0.0;
  for (int r = 0; r < q; ++r) {
    for (int s = r + 1; s < q; ++s) {
      if (G(r, s) == 1) continue;
      double acc = 0.0;
      for (int t = 0; t < r; ++t) acc += Psi(t, r) * Psi(t, s);
      const double v = -acc / Psi(r, r);
      Psi(r, s) = v;
      log_f -= 0.5 * v * v;
    }
  }
  double log_w = 0.0;
  for (int r = 0; r < q; ++r) {
    for (int s = r + 1; s < q; ++s) {
      if (G(r, s) != 1) continue;
      double acc = 0.0;
      for (int t = 0; t <= r; ++t) acc += Psi(t, r) * Psi(t, s);
      const double k_rs = acc / (2.0 * beta);
      log_w -= k_rs * k_rs / (2.0 * sig2);
    }
  }
  return log_f + log_w;
}

// One coupled (or independent) pair of weight averages and the implied log J.
// G_minus must not contain edge (i, j) (0-based, i < j).
static void akm_ratio_core(const arma::imat& G_minus, int i, int j,
                           double delta, double beta, double sigma, int N,
                           bool coupled,
                           double& mean_minus, double& mean_plus,
                           double& log_J,
                           arma::vec* draws_minus, arma::vec* draws_plus) {
  const int q = G_minus.n_rows;
  const double sig2 = sigma * sigma;
  const double b = 2.0 * delta + 2.0;
  arma::imat G_plus = G_minus;
  G_plus(i, j) = 1; G_plus(j, i) = 1;

  arma::vec nu(q, arma::fill::zeros);
  for (int l = 0; l < q; ++l)
    for (int k = l + 1; k < q; ++k) if (G_minus(l, k) == 1) nu(l) += 1.0;

  arma::mat Psi_m(q, q, arma::fill::zeros), Psi_p(q, q, arma::fill::zeros);
  double sum_m = -INFINITY, sum_p = -INFINITY;
  for (int t = 0; t < N; ++t) {
    Psi_m.zeros(); Psi_p.zeros();
    for (int l = 0; l < q; ++l) {
      const double c_low = R::rchisq(b + nu(l));
      Psi_m(l, l) = std::sqrt(c_low);
      if (coupled) {
        if (l == i) {
          const double z = R::norm_rand();
          Psi_p(l, l) = std::sqrt(c_low + z * z);
        } else {
          Psi_p(l, l) = Psi_m(l, l);
        }
      } else {
        Psi_p(l, l) = std::sqrt(R::rchisq(b + nu(l) + ((l == i) ? 1.0 : 0.0)));
      }
    }
    for (int r = 0; r < q; ++r) {
      for (int s = r + 1; s < q; ++s) {
        if (G_minus(r, s) == 1) {
          const double z = R::norm_rand();
          Psi_m(r, s) = z;
          Psi_p(r, s) = coupled ? z : R::norm_rand();
        } else if (r == i && s == j) {
          Psi_p(r, s) = R::norm_rand();
        }
      }
    }
    double lw_m = akm_log_wf(Psi_m, G_minus, beta, sig2);
    double lw_p = akm_log_wf(Psi_p, G_plus, beta, sig2);
    // The recursive completion can overflow for incomplete graphs at
    // moderate q (the completed entries compound row by row); such draws
    // carry weight zero. Accumulate on the log scale with a -Inf-safe
    // log-sum-exp so the ratio of log-means stays well defined whenever at
    // least one draw per arm has nonzero weight.
    if (std::isnan(lw_m)) lw_m = -INFINITY;
    if (std::isnan(lw_p)) lw_p = -INFINITY;
    if (lw_m != -INFINITY) {
      if (sum_m == -INFINITY) sum_m = lw_m;
      else if (sum_m > lw_m) sum_m += std::log1p(std::exp(lw_m - sum_m));
      else sum_m = lw_m + std::log1p(std::exp(sum_m - lw_m));
    }
    if (lw_p != -INFINITY) {
      if (sum_p == -INFINITY) sum_p = lw_p;
      else if (sum_p > lw_p) sum_p += std::log1p(std::exp(lw_p - sum_p));
      else sum_p = lw_p + std::log1p(std::exp(sum_p - lw_p));
    }
    if (draws_minus) (*draws_minus)(t) = lw_m;
    if (draws_plus) (*draws_plus)(t) = lw_p;
  }
  const double log_mean_m = sum_m - std::log(static_cast<double>(N));
  const double log_mean_p = sum_p - std::log(static_cast<double>(N));
  mean_minus = std::exp(log_mean_m);   // may underflow to 0; log_J does not
  mean_plus = std::exp(log_mean_p);
  const double a = delta + 1.0 + 0.5 * nu(i);
  log_J = 0.5 * std::log(2.0 * M_PI * sig2)
        + std::log(beta) - 0.5 * std::log(M_PI)
        + R::lgammafn(a) - R::lgammafn(a + 0.5)
        + log_mean_m - log_mean_p;
}

// [[Rcpp::export]]
Rcpp::List akm_ratio_cpp(const arma::imat& G_minus, int i, int j,
                         double delta, double beta, double sigma, int N,
                         bool coupled = true, bool return_draws = false) {
  // i, j are 1-based here (R convention); require i < j and (i, j) not in G-
  const int i0 = i - 1, j0 = j - 1;
  if (i0 < 0 || j0 <= i0 || j0 >= static_cast<int>(G_minus.n_rows))
    Rcpp::stop("need 1 <= i < j <= q");
  if (G_minus(i0, j0) != 0) Rcpp::stop("edge (i, j) must be absent from G_minus");
  if (N <= 0) Rcpp::stop("N must be positive");
  double mean_m, mean_p, log_J;
  arma::vec dm, dp;
  arma::vec *pm = nullptr, *pp = nullptr;
  if (return_draws) {
    dm.set_size(N); dp.set_size(N);
    pm = &dm; pp = &dp;
  }
  akm_ratio_core(G_minus, i0, j0, delta, beta, sigma, N, coupled,
                 mean_m, mean_p, log_J, pm, pp);
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("log_J") = log_J,
    Rcpp::Named("J") = std::exp(log_J),
    Rcpp::Named("mean_minus") = mean_m,
    Rcpp::Named("mean_plus") = mean_p);
  if (return_draws) {
    out["draws_minus"] = dm;
    out["draws_plus"] = dp;
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List sampler_chain_cpp(const arma::mat& S, double n_obs, int q,
                             double delta, double beta, double sigma,
                             double p_inc, int spec,
                             const Rcpp::NumericVector& log_z,
                             int n_iter, int n_burn, int refresh_every,
                             bool track_masks, int n_aux = 0,
                             int stat_thin = 0, bool time_blocks = false) {
  const double sig2 = sigma * sigma;
  using clk = std::chrono::steady_clock;
  double t_sweep = 0.0, t_toggle = 0.0, t_accept = 0.0;

  // edge list: upper triangle in row-major order (0,1),(0,2),...,(1,2),...,
  // which is the order of R's combn(q, 2); edge k is bit k of the graph mask
  std::vector<int> ei, ej;
  for (int a = 0; a < q - 1; ++a)
    for (int b = a + 1; b < q; ++b) { ei.push_back(a); ej.push_back(b); }
  const int M = static_cast<int>(ei.size());

  // input validation
  if (spec < 0 || spec > 2)
    Rcpp::stop("spec must be 0 (joint), 1 (hierarchical, exact log_z table) or 2 (hierarchical, coupled estimator)");
  if (static_cast<int>(S.n_rows) != q || static_cast<int>(S.n_cols) != q)
    Rcpp::stop("S must be a q x q matrix");
  if (spec == 2 && n_aux <= 0)
    Rcpp::stop("spec = 2 requires n_aux > 0 auxiliary draws per candidate");
  if ((spec == 1 || track_masks) && M > 62)
    Rcpp::stop("spec = 1 and track_masks index graphs by a 64-bit mask and require q(q-1)/2 <= 62");
  if (spec == 1 && static_cast<long long>(log_z.size()) != (1LL << M))
    Rcpp::stop("spec = 1 requires log_z of length 2^(q(q-1)/2), indexed by graph mask + 1");
  if (n_burn < 0 || n_burn >= n_iter)
    Rcpp::stop("need 0 <= n_burn < n_iter");
  if (!(p_inc > 0.0 && p_inc < 1.0))
    Rcpp::stop("p_inc must lie strictly between 0 and 1");

  arma::imat G(q, q, arma::fill::zeros);
  arma::mat K = arma::eye(q, q) * ((delta + 1.0) / beta);
  arma::mat Sigma = arma::inv_sympd(K);

  unsigned long long mask = 0ULL;
  std::vector<double> mask_visits;
  if (track_masks) mask_visits.assign(1ULL << M, 0.0);

  arma::mat edge_freq(q, q, arma::fill::zeros);
  arma::vec m_hist(M + 1, arma::fill::zeros);
  long long n_att = 0, n_acc = 0, n_disc = 0;
  double max_drift = 0.0;
  double sum_m = 0.0;
  long long n_keep = 0;
  // optional thinned scalar records (SBC): m, k11, k12, log|K|
  std::vector<double> rec_m, rec_k11, rec_k12, rec_logdet;

  arma::mat W2(2, 2), U(q, 2);

  for (int it = 0; it < n_iter; ++it) {
    // ---- row-block sweep ----
    clk::time_point tb0;
    if (time_blocks) tb0 = clk::now();
    for (int l = 0; l < q; ++l) {
      arma::uvec Nl = arma::find(G.col(l) == 1);
      const double Sll = S(l, l);
      const double psi = R::rgamma(delta + n_obs / 2.0 + 1.0,
                                   1.0 / (beta + 0.5 * Sll));
      arma::vec r_new(q, arma::fill::zeros);
      if (Nl.n_elem == 0) {
        r_new(l) = psi;
      } else {
        arma::vec sl = Sigma.col(l);
        arma::vec slN = sl.elem(Nl);
        arma::mat C = Sigma.submat(Nl, Nl) - slN * slN.t() / Sigma(l, l);
        C = 0.5 * (C + C.t());
        arma::mat Q = arma::eye(Nl.n_elem, Nl.n_elem) / sig2 + (2.0 * beta + Sll) * C;
        arma::mat L;
        bool okc = arma::chol(L, 0.5 * (Q + Q.t()), "lower");
        if (!okc) {
          // Sigma has drifted (near-singular state): hard refresh and retry.
          arma::mat Sfx;
          bool oki = arma::inv_sympd(Sfx, 0.5 * (K + K.t()));
          if (!oki) oki = arma::inv(Sfx, 0.5 * (K + K.t()));
          if (!oki) Rcpp::stop("Sigma refresh failed: K numerically singular");
          Sigma = 0.5 * (Sfx + Sfx.t());
          sl = Sigma.col(l);
          slN = sl.elem(Nl);
          C = Sigma.submat(Nl, Nl) - slN * slN.t() / Sigma(l, l);
          C = 0.5 * (C + C.t());
          Q = arma::eye(Nl.n_elem, Nl.n_elem) / sig2 + (2.0 * beta + Sll) * C;
          okc = arma::chol(L, 0.5 * (Q + Q.t()), "lower");
          if (!okc) Rcpp::stop("chol failed after Sigma refresh");
        }
        arma::vec sN(Nl.n_elem);
        for (arma::uword t = 0; t < Nl.n_elem; ++t) sN(t) = S(l, Nl(t));
        arma::vec mu = -arma::solve(arma::trimatu(L.t()),
                                    arma::solve(arma::trimatl(L), sN));
        arma::vec z(Nl.n_elem);
        for (arma::uword t = 0; t < Nl.n_elem; ++t) z(t) = R::norm_rand();
        arma::vec k = mu + arma::solve(arma::trimatu(L.t()), z);
        for (arma::uword t = 0; t < Nl.n_elem; ++t) r_new(Nl(t)) = k(t);
        r_new(l) = psi + arma::dot(k, C * k);
      }
      // symmetric rank-two Sigma update for the row replacement
      arma::vec r_old = K.col(l);
      arma::vec d = r_new - r_old;
      if (arma::norm(d, "inf") > 0.0) {
        arma::vec u = d;
        u(l) -= 0.5 * d(l);
        U.zeros(); U(l, 0) = 1.0; U.col(1) = u;
        W2 = {{0.0, 1.0}, {1.0, 0.0}};
        K.col(l) = r_new; K.row(l) = r_new.t();
        smw_rank2(Sigma, U, W2);
      }
      // guard: a near-zero slack makes K near-singular and the rank-two
      // arithmetic ill-conditioned; refresh Sigma exactly in that case
      if (psi < 1e-6) {
        arma::mat Sfx;
        bool oki = arma::inv_sympd(Sfx, 0.5 * (K + K.t()));
        if (!oki) oki = arma::inv(Sfx, 0.5 * (K + K.t()));
        if (oki) Sigma = 0.5 * (Sfx + Sfx.t());
      }
    }

    // ---- edge toggles, random scan over all pairs ----
    if (time_blocks) {
      const clk::time_point tb1 = clk::now();
      t_sweep += std::chrono::duration<double>(tb1 - tb0).count();
      tb0 = tb1;
    }
    arma::uvec ord = arma::randperm(M);
    for (int t = 0; t < M; ++t) {
      const int k = ord(t);
      const int i = ei[k], j = ej[k];
      const double Sii_ = Sigma(i, i), Sjj_ = Sigma(j, j), Sij_ = Sigma(i, j);
      const double D = Sii_ * Sjj_ - Sij_ * Sij_;
      const double c2 = std::sqrt(Sjj_ / D);
      const double phi = -Sij_ / std::sqrt(Sjj_ * D);
      const double c1 = K(i, j) - c2 * phi;
      const double c3 = K(j, j) - phi * phi;
      const double lam = 1.0 / sig2 + (2.0 * beta + S(j, j)) / (c2 * c2);
      const double lin = -S(i, j) + (2.0 * beta + S(j, j)) * c1 / (c2 * c2);
      const double mu_x = lin / lam;
      const double log_q0 = 0.5 * std::log(lam) - 0.5 * LOG2PI - 0.5 * lam * mu_x * mu_x;
      const double log_pslab0 = -0.5 * LOG2PI - std::log(sigma);
      double log_J = 0.0;
      if (spec == 1) {
        const unsigned long long bit = 1ULL << k;
        const unsigned long long m_without = (G(i, j) == 1) ? (mask - bit) : mask;
        log_J = log_z[static_cast<R_xlen_t>(m_without)]
              - log_z[static_cast<R_xlen_t>(m_without + bit)];
      } else if (spec == 2) {
        // hierarchical with the coupled Monte Carlo estimator (n_aux draws)
        arma::imat Gm = G;
        Gm(i, j) = 0; Gm(j, i) = 0;
        double mean_m, mean_p;
        akm_ratio_core(Gm, i, j, delta, beta, sigma, n_aux, true,
                       mean_m, mean_p, log_J, nullptr, nullptr);
      }
      const double log_A = std::log(p_inc / (1.0 - p_inc)) + log_pslab0 - log_q0 + log_J;
      ++n_att;
      const double log_u = std::log(R::unif_rand());
      // decision-discrepancy probe: with spec = 1 and n_aux > 0, also compute
      // the estimated-J decision with the SAME uniform along the exact-ratio
      // trajectory; the chain itself always follows the exact decision.
      if (spec == 1 && n_aux > 0) {
        arma::imat Gm = G;
        Gm(i, j) = 0; Gm(j, i) = 0;
        double mean_m, mean_p, lj_est;
        akm_ratio_core(Gm, i, j, delta, beta, sigma, n_aux, true,
                       mean_m, mean_p, lj_est, nullptr, nullptr);
        const double log_A_est = std::log(p_inc / (1.0 - p_inc))
                               + log_pslab0 - log_q0 + lj_est;
        const bool adding = (G(i, j) == 0);
        const bool d_exact = adding ? (log_u < log_A) : (log_u < -log_A);
        const bool d_est = adding ? (log_u < log_A_est) : (log_u < -log_A_est);
        if (d_exact != d_est) ++n_disc;
      }
      double a = 0.0, b = 0.0;
      bool accept = false;
      if (G(i, j) == 0) {
        if (log_u < log_A) {
          const double xs = mu_x + R::norm_rand() / std::sqrt(lam);
          a = xs;
          b = c3 + std::pow((xs - c1) / c2, 2) - K(j, j);
          G(i, j) = G(j, i) = 1;
          K(i, j) = K(j, i) = xs;
          K(j, j) += b;
          if (spec == 1 || track_masks) mask |= (1ULL << k);
          accept = true;
        }
      } else {
        if (log_u < -log_A) {
          a = -K(i, j);
          b = c3 + std::pow(c1 / c2, 2) - K(j, j);
          G(i, j) = G(j, i) = 0;
          K(i, j) = K(j, i) = 0.0;
          K(j, j) += b;
          if (spec == 1 || track_masks) mask &= ~(1ULL << k);
          accept = true;
        }
      }
      if (accept) {
        ++n_acc;
        clk::time_point ta0;
        if (time_blocks) ta0 = clk::now();
        // K' = K + U W U^t, U = [e_i, e_j], W = [[0, a], [a, b]]
        U.zeros(); U(i, 0) = 1.0; U(j, 1) = 1.0;
        W2 = {{0.0, a}, {a, b}};
        smw_rank2(Sigma, U, W2);
        if (time_blocks)
          t_accept += std::chrono::duration<double>(clk::now() - ta0).count();
      }
    }
    if (time_blocks)
      t_toggle += std::chrono::duration<double>(clk::now() - tb0).count();

    // ---- periodic full refresh of Sigma (drift control) ----
    if (refresh_every > 0 && (it + 1) % refresh_every == 0) {
      arma::mat Sig_exact;
      bool ok = arma::inv_sympd(Sig_exact, 0.5 * (K + K.t()));
      if (!ok) ok = arma::inv(Sig_exact, 0.5 * (K + K.t()));
      if (ok) {
        const double drift = arma::abs(Sigma - Sig_exact).max();
        if (drift > max_drift) max_drift = drift;
        Sigma = Sig_exact;
      }
    }

    // ---- accumulate ----
    if (it >= n_burn) {
      ++n_keep;
      int m_now = 0;
      for (int k = 0; k < M; ++k) if (G(ei[k], ej[k]) == 1) {
        ++m_now;
        edge_freq(ei[k], ej[k]) += 1.0;
      }
      m_hist(m_now) += 1.0;
      sum_m += m_now;
      if (track_masks) mask_visits[mask] += 1.0;
      if (stat_thin > 0 &&
          static_cast<long long>(it - n_burn + 1) % stat_thin == 0) {
        rec_m.push_back(static_cast<double>(m_now));
        rec_k11.push_back(K(0, 0));
        rec_k12.push_back(K(0, 1));
        double ld = 0.0, sign = 1.0;
        arma::log_det(ld, sign, 0.5 * (K + K.t()));
        rec_logdet.push_back(ld);
      }
    }
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("edge_freq") = edge_freq / static_cast<double>(n_keep),
    Rcpp::Named("m_hist") = m_hist,
    Rcpp::Named("E_m") = sum_m / static_cast<double>(n_keep),
    Rcpp::Named("acc_rate") = static_cast<double>(n_acc) / static_cast<double>(n_att),
    Rcpp::Named("max_sigma_drift") = max_drift);
  if (track_masks) out["mask_visits"] = Rcpp::NumericVector(mask_visits.begin(), mask_visits.end());
  if (stat_thin > 0) {
    out["rec_m"] = Rcpp::NumericVector(rec_m.begin(), rec_m.end());
    out["rec_k11"] = Rcpp::NumericVector(rec_k11.begin(), rec_k11.end());
    out["rec_k12"] = Rcpp::NumericVector(rec_k12.begin(), rec_k12.end());
    out["rec_logdet"] = Rcpp::NumericVector(rec_logdet.begin(), rec_logdet.end());
  }
  if (time_blocks) {
    out["timing"] = Rcpp::NumericVector::create(
      Rcpp::Named("sweep_s") = t_sweep,
      Rcpp::Named("toggle_s") = t_toggle,
      Rcpp::Named("accept_s") = t_accept,
      Rcpp::Named("n_att") = static_cast<double>(n_att),
      Rcpp::Named("n_acc") = static_cast<double>(n_acc));
  }
  if (spec == 1 && n_aux > 0)
    out["disc_rate"] = static_cast<double>(n_disc) / static_cast<double>(n_att);
  return out;
}
