#' Compute the Log Marginal Likelihood
#'
#' This function calculates the log marginal likelihood based on "group" and "team" assignments.
#'
#' @param Y A numeric matrix (n x p), representing observed data.
#' @param cluster_assign An integer vector of length `n`, assigning each observation to a "group".
#' @param team_assign An integer vector of length `n`, assigning each observation to a "team".
#' @param k Integer, the number of spatial clusters ("groups").
#' @param j Integer, the number of "teams".
#' @param sigmasq_mu Numeric, variance parameter controlling spatial clusters variability.
#' @param Sigma A numeric (p x p) covariance matrix for the observed data.
#'
#' @returns A numeric value representing the log marginal likelihood.
#'
#' @keywords internal
#'
#' @noRd
.MarginalLikelihood_R <- function(Y, cluster_assign, team_assign, k, j, sigmasq_mu, Sigma){
  ### Here we calculate the log marginal likelihood ###
  n = nrow(Y)
  p = ncol(Y)

  # Build one-hot matrices:
  # X: n x k, each row (observation) assigned to one cluster.
  X <- table(seq_along(cluster_assign), cluster_assign)

  # Z: k x j, each cluster assigned to one team.
  Z <- table(seq_along(team_assign), team_assign)

  # omega = crossprod(Z, Z) + sigmasq_mu * diag(1, nrow = j)
  # omega_inv = solve(omega)
  # log_det_omega <- sum(log(diag(omega)))

  omega_diag <- colSums(Z) + sigmasq_mu
  omega_inv <- diag(1 / omega_diag, nrow = length(omega_diag))
  log_det_omega <- sum(log(omega_diag))

  delta = (1/sigmasq_mu) * diag(1, nrow = k) + crossprod(X, X) - (1/sigmasq_mu) * Z %*% tcrossprod(omega_inv, Z)
  chol_delta <- chol(delta)
  log_det_delta <- 2 * sum(log(diag(chol_delta)))
  delta_inv <- chol2inv(chol_delta)

  chol_Sigma <- chol(Sigma)
  Sigma_inv = chol2inv(chol_Sigma)
  log_det_Sigma <- 2 * sum(log(diag(chol_Sigma)))

  exp_term = -(1/2) * sum(diag(tcrossprod(Sigma_inv, Y) %*% (diag(1, nrow = n) - X %*% tcrossprod(delta_inv,X)) %*% Y))

  marg_like <- -(n*p/2)*log(2*pi) + (p*(j-k)/2)*log(sigmasq_mu) - (n/2)*log_det_Sigma -
    (p/2)*log_det_omega - (p/2) * log_det_delta + exp_term

  return(marg_like)
}
#' Compute the Log Marginal Likelihood (fast C++ backend)
#' @keywords internal
#' @noRd
.MarginalLikelihood <- function(Y, cluster_assign, team_assign,
                                k, j, sigmasq_mu, Sigma)
{
  MarginalLikelihood_cpp(
    Y,
    as.integer(cluster_assign),
    as.integer(team_assign),
    k, j, sigmasq_mu, Sigma
  )
}
