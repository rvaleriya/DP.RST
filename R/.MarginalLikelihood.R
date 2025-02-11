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
.MarginalLikelihood <- function(Y, cluster_assign, team_assign, k, j, sigmasq_mu, Sigma){
  ### Here we calculate the log marginal likelihood ###
  n = nrow(Y)
  p = ncol(Y)

  # Create a binary matrix from the cluster assignments
  X <- table(sequence(length(cluster_assign)), cluster_assign)
  # Create a binary matrix from the teams assignments
  Z <- table(sequence(length(team_assign)), team_assign)

  # Calculate terms of the marginal log-likelihood
  omega = crossprod(Z, Z) + sigmasq_mu * diag(1, nrow = j)
  omega_inv = solve(omega)

  delta = (1/sigmasq_mu) * diag(1, nrow = k) + crossprod(X, X) - (1/sigmasq_mu) * Z %*% tcrossprod(omega_inv, Z)
  delta_inv = solve(delta)

  Sigma_inv = solve(Sigma)

  exp_term = -(1/2) * sum(diag(tcrossprod(Sigma_inv, Y) %*% (diag(1, nrow = n) - X %*% tcrossprod(delta_inv,X)) %*% Y))

  marg_like <- -(n*p/2)*log(2*pi) + (p*(j-k)/2)*log(sigmasq_mu) - (n/2)*log(det(Sigma)) -
    (p/2)*log(det(omega)) - (p/2) * log(det(delta)) + exp_term

  return(marg_like)
}
