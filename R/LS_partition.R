#' Least Squares Partition from DP-RST Output
#'
#' This function selects the partition with the lowest Frobenius norm from the
#' posterior samples of the DP-RST model.
#'
#' @param DP.RST_output A list containing the output from `DP.RST()`,
#' including `cluster_out`, `teams_out`, and `j_teams_out`.
#'
#' @returns A list with the least-squares partitions:
#' \describe{
#'   \item{`groups_partition`}{Vector of group assignments for the optimal iteration.}
#'   \item{`teams_partition`}{Vector of refined team assignments for the optimal iteration.}
#' }
#'
#' @export
LS_partition <- function(DP.RST_output) {

  ##### CHOOSE THE ITERATION #####
  groups_assign = DP.RST_output$cluster_out
  teams_assign = DP.RST_output$teams_out
  teams_number = DP.RST_output$j_teams_out

  n = length(groups_assign)

  # List to save co-clustering matrices
  W_cum <- matrix(0, nrow = n, ncol = n)

  # Get the mode of the team
  tbl <- table(teams_number)
  max_freq <- max(tbl)
  modes <- as.numeric(names(tbl[tbl == max_freq]))

  # Get the indices where teams_number equals the mode
  mode_indices <- which(teams_number == modes)

  # Run the loop to create co-clustering matrix for over the selected indices
  for (i in 1:length(mode_indices)) {
    # Binary matrix for groups and teams
    X <- table(sequence(length(groups_assign[mode_indices[i], ])), groups_assign[mode_indices[i], ])
    Z <- table(sequence(length(teams_assign[[mode_indices[i]]])), teams_assign[[mode_indices[i]]])

    # Get the membership of observations in each team
    obs_in_teams <- X %*% Z

    # Compute W such that w_ij = 1 if observation i and observation j share same team
    W <- obs_in_teams %*% t(obs_in_teams)

    W_cum <- W_cum + W
  }

  mean_matrix <- W_cum/length(mode_indices)

  # Calculate the vector with Frobenius norms
  norm = c()

  for (i in 1:length(mode_indices)) {
    # Binary matrix for groups and teams
    X <- table(sequence(length(groups_assign[mode_indices[i], ])), groups_assign[mode_indices[i], ])
    Z <- table(sequence(length(teams_assign[[mode_indices[i]]])), teams_assign[[mode_indices[i]]])

    # Get the membership of observations in each team
    obs_in_teams <- X %*% Z
    # Compute W such that w_ij = 1 if observation i and observation j share same team
    mat <- obs_in_teams %*% t(obs_in_teams)

    diff_matrix <- mat - mean_matrix
    norm[i] <- norm(diff_matrix, type = "F")
  }

  # Fins the minimum Frobenius norm index
  min_index <- mode_indices[which.min(norm)]

  groups_assign_out <- groups_assign[min_index, ]
  teams_assign_out <- teams_assign[[min_index]]

  return(list(groups_partition = groups_assign_out, teams_partition = teams_assign_out))
}
