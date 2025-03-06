#' Select Optimal Partition from DP-RST Output
#'
#' This function selects an optimal partition from the posterior samples of a
#' Dirichlet Process Random Spanning Tree (DP-RST) model using various methods
#' to summarize the posterior distribution of partitions.
#'
#' @param DP.RST_output A list containing the output from `DP.RST()`.
#' @param method Character string specifying the method for selecting the optimal partition:
#'   \describe{
#'     \item{`"mode_based"`}{Selects a partition from iterations with the most common number
#'           of "teams", minimizing Frobenius distance to the mean within that mode. This is
#'           the default method.}
#'     \item{`"frobenius"`}{Selects the partition closest to the posterior similarity matrix
#'           across all iterations, as measured by Frobenius norm.}
#'     \item{`"binder"`}{Minimizes Binder's loss function against the posterior similarity
#'           matrix, which often better represents the central tendency of the posterior.}
#'     \item{`"posterior_mode"`}{Selects the iteration with the highest marginal likelihood
#'           value, if available.}
#'   }
#' @param batch_size Integer specifying the number of iterations to process in each batch.
#'        Larger values may increase speed but require more memory. Default is 100.
#'
#' @return A list containing the selected partition:
#'   \describe{
#'     \item{`groups_partition`}{Vector of group assignments from the selected iteration.}
#'     \item{`teams_partition`}{Vector of refined team assignments from the selected iteration.}
#'     \item{`selected_iteration`}{Index of the selected iteration.}
#'   }
#' @export
partition <- function(DP.RST_output, method = "mode_based", batch_size = 100) {
  # Extract relevant outputs
  groups_assign <- DP.RST_output$cluster_out
  teams_assign <- DP.RST_output$teams_out
  teams_number <- DP.RST_output$j_teams_out

  n <- ncol(groups_assign)
  total_iterations <- length(teams_number)
  all_iterations <- 1:total_iterations

  # Handle mode_based approach (original approach)
  if (method == "mode_based") {
    # Get the mode of the team numbers
    tbl <- table(teams_number)
    max_freq <- max(tbl)
    modes <- as.numeric(names(tbl[tbl == max_freq]))
    # Get the indices where teams_number equals the mode
    mode_indices <- which(teams_number == modes)

    # Create co-clustering matrix in batches
    W_cum <- matrix(0, nrow = n, ncol = n)
    batch_count <- ceiling(length(mode_indices) / batch_size)

    for (b in 1:batch_count) {
      batch_start <- (b - 1) * batch_size + 1
      batch_end <- min(b * batch_size, length(mode_indices))
      batch_indices <- mode_indices[batch_start:batch_end]

      for (idx in batch_indices) {
        X <- table(sequence(length(groups_assign[idx, ])), groups_assign[idx, ])
        Z <- table(sequence(length(teams_assign[[idx]])), teams_assign[[idx]])
        obs_in_teams <- X %*% Z
        W_cum <- W_cum + tcrossprod(obs_in_teams)

        # Clean up
        rm(X, Z, obs_in_teams)
      }
      gc()
    }

    mean_matrix <- W_cum / length(mode_indices)

    # Find minimum Frobenius norm in batches
    min_norm <- Inf
    best_idx <- mode_indices[1]

    for (b in 1:batch_count) {
      batch_start <- (b - 1) * batch_size + 1
      batch_end <- min(b * batch_size, length(mode_indices))
      batch_indices <- mode_indices[batch_start:batch_end]

      for (idx in batch_indices) {
        X <- table(sequence(length(groups_assign[idx, ])), groups_assign[idx, ])
        Z <- table(sequence(length(teams_assign[[idx]])), teams_assign[[idx]])
        obs_in_teams <- X %*% Z
        curr_W <- tcrossprod(obs_in_teams)
        current_norm <- sqrt(sum((curr_W - mean_matrix)^2))

        if (current_norm < min_norm) {
          min_norm <- current_norm
          best_idx <- idx
        }

        # Clean up
        rm(X, Z, obs_in_teams, curr_W)
      }
      gc()
    }
  } else { # Fixed: else must be on the same line as the closing brace
    if (method %in% c("frobenius", "binder")) {
      # Build PSM in batches using all iterations
      PSM <- matrix(0, nrow = n, ncol = n)
      batch_count <- ceiling(total_iterations / batch_size)

      cat("Building posterior similarity matrix in", batch_count, "batches...\n")

      for (b in 1:batch_count) {
        if (b %% 10 == 0 || b == 1) cat("Processing batch", b, "of", batch_count, "\n")

        batch_start <- (b - 1) * batch_size + 1
        batch_end <- min(b * batch_size, total_iterations)
        batch_iterations <- batch_start:batch_end

        for (idx in batch_iterations) {
          X <- table(sequence(length(groups_assign[idx, ])), groups_assign[idx, ])
          Z <- table(sequence(length(teams_assign[[idx]])), teams_assign[[idx]])
          obs_in_teams <- X %*% Z
          PSM <- PSM + tcrossprod(obs_in_teams)

          # Clean up
          rm(X, Z, obs_in_teams)
        }
        gc()
      }

      PSM <- PSM / total_iterations

      # Find the best partition among all iterations (in batches)
      min_loss <- Inf
      best_idx <- 1

      cat("Finding optimal partition in", batch_count, "batches...\n")

      for (b in 1:batch_count) {
        if (b %% 10 == 0 || b == 1) cat("Processing batch", b, "of", batch_count, "\n")

        batch_start <- (b - 1) * batch_size + 1
        batch_end <- min(b * batch_size, total_iterations)
        batch_iterations <- batch_start:batch_end

        for (idx in batch_iterations) {
          X <- table(sequence(length(groups_assign[idx, ])), groups_assign[idx, ])
          Z <- table(sequence(length(teams_assign[[idx]])), teams_assign[[idx]])
          obs_in_teams <- X %*% Z
          curr_W <- tcrossprod(obs_in_teams)

          # Different loss calculation based on method
          if (method == "frobenius") {
            # Frobenius norm
            current_loss <- sqrt(sum((curr_W - PSM)^2))
          } else if (method == "binder") { # Fixed: else must be on the same line
            # Binder loss function
            current_loss <- sum(PSM * (1 - curr_W)) + sum((1 - PSM) * curr_W)
          }

          if (current_loss < min_loss) {
            min_loss <- current_loss
            best_idx <- idx
          }

          # Clean up
          rm(X, Z, obs_in_teams, curr_W)
        }
        gc()
      }
    } else if (method == "posterior_mode") { # Fixed: else must be on the same line
      # Pick the iteration with highest marginal likelihood
      if (!is.null(DP.RST_output$marginal_likelihood_out)) {
        best_idx <- which.max(DP.RST_output$marginal_likelihood_out)
      } else { # Fixed: else must be on the same line
        warning("Marginal likelihood not available, falling back to random selection")
        best_idx <- sample(total_iterations, 1)
      }
    }
  }

  return(list(
    groups_partition = groups_assign[best_idx, ],
    teams_partition = teams_assign[[best_idx]],
    selected_iteration = best_idx
  ))
}


# # LS_partition_memory_efficient <- function(DP.RST_output, method = "mode_based", burn_in = 0, thin = 1, batch_size = 100) {
#   # Extract relevant outputs
#   groups_assign <- DP.RST_output$cluster_out
#   teams_assign <- DP.RST_output$teams_out
#   teams_number <- DP.RST_output$j_teams_out
#
#   n <- ncol(groups_assign)
#   total_iterations <- length(teams_number)
#
#   # Handle mode_based approach
#   if (method == "mode_based") {
#     # Get the mode of the team numbers
#     tbl <- table(teams_number)
#     max_freq <- max(tbl)
#     modes <- as.numeric(names(tbl[tbl == max_freq]))
#     # Get the indices where teams_number equals the mode
#     mode_indices <- which(teams_number == modes)
#
#     # Create co-clustering matrix for mode iterations
#     W_cum <- matrix(0, nrow = n, ncol = n)
#     for (i in 1:length(mode_indices)) {
#       X <- table(sequence(length(groups_assign[mode_indices[i], ])), groups_assign[mode_indices[i], ])
#       Z <- table(sequence(length(teams_assign[[mode_indices[i]]])), teams_assign[[mode_indices[i]]])
#       # Get the membership of observations in each team
#       obs_in_teams <- X %*% Z
#       # Compute W efficiently without storing intermediate W
#       W_cum <- W_cum + tcrossprod(obs_in_teams)
#
#       # Clean up memory
#       rm(X, Z, obs_in_teams)
#       if (i %% 10 == 0) gc()
#     }
#     mean_matrix <- W_cum/length(mode_indices)
#
#     # Calculate the best partition with minimum Frobenius norm
#     min_norm <- Inf
#     best_idx <- mode_indices[1]
#
#     for (i in 1:length(mode_indices)) {
#       X <- table(sequence(length(groups_assign[mode_indices[i], ])), groups_assign[mode_indices[i], ])
#       Z <- table(sequence(length(teams_assign[[mode_indices[i]]])), teams_assign[[mode_indices[i]]])
#       # Get the membership of observations in each team
#       obs_in_teams <- X %*% Z
#       # Compute W
#       curr_W <- tcrossprod(obs_in_teams)
#       # Calculate Frobenius norm directly without storing diff_matrix
#       current_norm <- sqrt(sum((curr_W - mean_matrix)^2))
#
#       if (current_norm < min_norm) {
#         min_norm <- current_norm
#         best_idx <- mode_indices[i]
#       }
#
#       # Clean up memory
#       rm(X, Z, obs_in_teams, curr_W)
#       if (i %% 10 == 0) gc()
#     }
#
#     return(list(
#       groups_partition = groups_assign[best_idx, ],
#       teams_partition = teams_assign[[best_idx]],
#       selected_iteration = best_idx
#     ))
#   }
#
#   # For other methods (memory-efficient versions)
#   else {
#     # Calculate valid iterations after burn-in and thinning
#     valid_iterations <- seq(burn_in + 1, total_iterations, by = thin)
#
#     if (method %in% c("frobenius_lite", "binder_lite")) {
#       # Memory-efficient version that uses sampling and batching
#       # Randomly sample a subset of iterations to build PSM if dataset is large
#       if (length(valid_iterations) > 1000) {
#         sample_size <- min(1000, length(valid_iterations))
#         sample_iterations <- sample(valid_iterations, sample_size)
#       } else {
#         sample_iterations <- valid_iterations
#       }
#
#       # Build PSM in batches
#       PSM <- matrix(0, nrow = n, ncol = n)
#       batch_count <- ceiling(length(sample_iterations) / batch_size)
#
#       for (b in 1:batch_count) {
#         batch_start <- (b - 1) * batch_size + 1
#         batch_end <- min(b * batch_size, length(sample_iterations))
#         batch_iterations <- sample_iterations[batch_start:batch_end]
#
#         for (idx in batch_iterations) {
#           X <- table(sequence(length(groups_assign[idx, ])), groups_assign[idx, ])
#           Z <- table(sequence(length(teams_assign[[idx]])), teams_assign[[idx]])
#           obs_in_teams <- X %*% Z
#           PSM <- PSM + tcrossprod(obs_in_teams)
#
#           # Clean up
#           rm(X, Z, obs_in_teams)
#         }
#         gc()
#       }
#
#       PSM <- PSM / length(sample_iterations)
#
#       # Find the best partition among all valid iterations
#       min_loss <- Inf
#       best_idx <- valid_iterations[1]
#
#       # Process in batches
#       for (b in 1:ceiling(length(valid_iterations) / batch_size)) {
#         batch_start <- (b - 1) * batch_size + 1
#         batch_end <- min(b * batch_size, length(valid_iterations))
#         batch_iterations <- valid_iterations[batch_start:batch_end]
#
#         for (idx in batch_iterations) {
#           X <- table(sequence(length(groups_assign[idx, ])), groups_assign[idx, ])
#           Z <- table(sequence(length(teams_assign[[idx]])), teams_assign[[idx]])
#           obs_in_teams <- X %*% Z
#           curr_W <- tcrossprod(obs_in_teams)
#
#           # Different loss calculation based on method
#           if (method == "frobenius_lite") {
#             # Frobenius norm
#             current_loss <- sqrt(sum((curr_W - PSM)^2))
#           } else if (method == "binder_lite") {
#             # Binder loss function
#             # For each pair i,j: Binder loss penalizes:
#             # - putting i,j in same cluster when they should be separate
#             # - putting i,j in different clusters when they should be together
#             current_loss <- sum(PSM * (1 - curr_W)) + sum((1 - PSM) * curr_W)
#           }
#
#           if (current_loss < min_loss) {
#             min_loss <- current_loss
#             best_idx <- idx
#           }
#
#           # Clean up
#           rm(X, Z, obs_in_teams, curr_W)
#         }
#         gc()
#       }
#     }
#
#     else if (method == "posterior_mode") {
#       # Simply pick the iteration with highest marginal likelihood
#       if (!is.null(DP.RST_output$marginal_likelihood_out)) {
#         ml_values <- DP.RST_output$marginal_likelihood_out[valid_iterations]
#         best_idx <- valid_iterations[which.max(ml_values)]
#       } else {
#         warning("Marginal likelihood not available, falling back to random selection.")
#         best_idx <- sample(valid_iterations, 1)
#       }
#     }
#
#     return(list(
#       groups_partition = groups_assign[best_idx, ],
#       teams_partition = teams_assign[[best_idx]],
#       selected_iteration = best_idx
#     ))
#   }
# }
