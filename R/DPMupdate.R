#' Perform Dirichlet Process Mixture (DPM) Update
#'
#' @param mu_m Matrix of group means.
#' @param teams_m Vector of current team assignments.
#' @param sigmasq_y_m Covariance matrix for observations.
#' @param sigmasq_mu Scaling constant for the covariance matrix.
#' @param alpha Concentration parameter for the Dirichlet process.
#'
#' @returns A list containing updated team assignments (`teams_m`) and the number of teams (`j_teams_m`).
#' @keywords internal
#'
#' @noRd
.DPM_update <- function(mu_m, teams_m, sigmasq_y_m, sigmasq_mu, alpha) {

  p = ncol(sigmasq_y_m)

  k_m <- length(teams_m)
  j_teams_m <- max(teams_m)

  sigmasq_y_m_inv <- chol2inv(chol(sigmasq_y_m))
  group_var_inv <- (1/sigmasq_mu) * sigmasq_y_m_inv
  B_var_new <- (1 + sigmasq_mu) * sigmasq_y_m

  for (ind in 1:k_m) {
    # Remove the k-th group
    mu_subset <- mu_m[-ind, , drop = FALSE]
    mu_k <- matrix(mu_m[ind,], ncol = p)

    # Update team assignments
    team_assign_subset <- dense_rank(teams_m[-ind])
    Z_subset <- table(sequence(length(team_assign_subset)), team_assign_subset)
    J_current <- max(team_assign_subset)

    # Initialize variables
    n_j <- numeric(J_current)
    prob_existing_table <- numeric(J_current)

    # Compute probabilities for existing teams
    for (j in 1:J_current) {
      n_j[j] <- sum(team_assign_subset == j)
      mu_subset_j <- mu_subset[which(Z_subset[, j] == 1), , drop = FALSE]

      delta_star_inv <- sigmasq_y_m_inv + n_j[j] * group_var_inv
      inter_mat <- chol2inv(chol(group_var_inv + delta_star_inv))
      theta_star <- chol2inv(chol(delta_star_inv)) %*% (n_j[j] * group_var_inv %*% colMeans(mu_subset_j))
      B_var <- chol2inv(chol(group_var_inv - group_var_inv %*% inter_mat %*% group_var_inv))
      beta_means <- B_var %*% (group_var_inv %*% inter_mat %*% delta_star_inv %*% theta_star)

      prob_existing_table[j] <- log(n_j[j]) + dmvnorm(mu_k, mu = beta_means, sigma = B_var, logged = TRUE)
    }

    # Compute probability for a new team
    prob_new_table <- log(alpha) + dmvnorm(mu_k, mu = rep(0, p), sigma = B_var_new, logged = TRUE)

    # Normalize probabilities
    probability_unnormalized <- exp(c(prob_existing_table, prob_new_table) - max(c(prob_existing_table, prob_new_table)))
    probability <- probability_unnormalized / sum(probability_unnormalized)

    # Sample new team assignment
    teams_m <- append(team_assign_subset, rcat(n = 1, prob = probability), after = (ind - 1))
    teams_m <- dense_rank(teams_m)
    j_teams_m <- max(teams_m)
  }

  return(list(teams_m = teams_m, j_teams_m = j_teams_m))
}
