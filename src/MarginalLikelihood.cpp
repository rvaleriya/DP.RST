// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double MarginalLikelihood_cpp(const arma::mat&  Y,
                              const Rcpp::IntegerVector& cluster_assign,
                              const Rcpp::IntegerVector& team_assign,
                              int                       k,
                              int                       j,
                              double                    sigmasq_mu,
                              const arma::mat&          Sigma)
{
  const int n = Y.n_rows;
  const int p = Y.n_cols;

  /* ---------- 1. Build one-hot matrices X (n×k) and Z (k×j) ---------- */
  arma::mat X(n, k, arma::fill::zeros);
  for (int i = 0; i < n; ++i)
    X(i, cluster_assign[i] - 1) = 1.0;       // R is 1-based

  arma::mat Z(k, j, arma::fill::zeros);
  for (int c = 0; c < k; ++c)
    Z(c, team_assign[c] - 1) = 1.0;

  /* ---------- 2. Ω diagonal, its inverse, and log-det ---------- */
  arma::rowvec omega_diag = arma::sum(Z, 0) + sigmasq_mu;   // 1×j
  arma::mat    omega_inv  = arma::diagmat(1.0 / omega_diag.t());  // j×j
  const double log_det_omega = arma::accu(arma::log(omega_diag));

  /* ---------- 3. Δ, its inverse, and log-det ---------- */
  arma::mat delta =
  (1.0 / sigmasq_mu) * arma::eye<arma::mat>(k, k)
    + X.t() * X
  - (1.0 / sigmasq_mu) * Z * omega_inv * Z.t();

  arma::mat chol_delta   = arma::chol(delta);
  const double log_det_delta = 2.0 * arma::accu(arma::log(chol_delta.diag()));
  arma::mat delta_inv    = arma::inv_sympd(delta);    // uses chol under the hood

  /* ---------- 4. Σ inverse and log-det via Cholesky ---------- */
  arma::mat chol_Sigma   = arma::chol(Sigma);
  const double log_det_Sigma = 2.0 * arma::accu(arma::log(chol_Sigma.diag()));
  arma::mat Sigma_inv    = arma::inv_sympd(Sigma);    // stable for SPD

  /* ---------- 5. Trace term for exponent ---------- */
  arma::mat resid_mat = arma::eye<arma::mat>(n, n) - X * delta_inv * X.t();
  double trace_term   = arma::trace(Sigma_inv * ( Y.t() * resid_mat * Y ));
  double exp_term     = -0.5 * trace_term;

  /* ---------- 6. Assemble final log marginal likelihood ---------- */
  const double marg_like =
  -0.5 * n * p * std::log(2.0 * M_PI)
  + 0.5 * p * (j - k) * std::log(sigmasq_mu)
  - 0.5 * n * log_det_Sigma
  - 0.5 * p * log_det_omega
  - 0.5 * p * log_det_delta
  + exp_term;

  return marg_like;
}
