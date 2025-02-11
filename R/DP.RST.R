# Register internal functions to avoid R CMD check notes
utils::globalVariables(c(".getEdgeStatus", ".proposeMST", ".MarginalLikelihood",
                         ".mergeCluster", ".splitCluster", ".updateMerge",
                         ".updateSplit", "generate_pairs"))

#' DP-RST: Dirichlet Process Mixture of Random Spanning Trees
#'
#' Implements the DP-RST algorithm using MCMC sampling.
#'
#' @param Y Standardized matrix of Principal Components.
#' @param graph0 Spatial graph. Should be an `igraph` object containing `length(Y)` vertices.
#' @param init_val Named list of initial values containing:
#'        \describe{
#'          \item{`trees`}{List of initial spanning trees (M `igraph` objects).}
#'          \item{`cluster`}{Initial cluster membership matrix (`length(Y) * M`).}
#'          \item{`mu`}{List of initial cluster means.}
#'          \item{`teams`}{Initial refined partition matrix (`length(Y) * M`).}
#'          \item{`mu_teams`}{List of initial team means.}
#'          \item{`sigmasq_y`}{Initial noise variance.}
#'        }
#' @param hyperpar Named list of hyperparameters containing:
#'        \describe{
#'          \item{`M`}{Number of temperatures for Paralell Tempering. M=1 for one chain.}
#'          \item{`sigmasq_mu`}{Variance of Gaussian prior for `mu`.}
#'          \item{`lambda_s`}{Scale parameter of Inverse-Wishart prior.}
#'          \item{`nu`}{Degrees of freedom for Inverse-Wishart prior.}
#'          \item{`k_max`}{Maximum number of clusters in spatial partition.}
#'          \item{`j_max`}{Maximum number of clusters in refined partition.}
#'          \item{`temp`}{A vector of temperatures of size M for Parallel Tempering.}
#'        }
#' @param MCMC Number of MCMC iterations.
#' @param BURNIN Number of burn-in iterations.
#' @param THIN Thinning interval (retains samples every `THIN` iterations).
#' @param PT Logical; whether to use Parallel Tempering.
#' @param seed Random seed.
#'
#' @returns A list of MCMC samples:
#'         \describe{
#'          \item{`cluster_out`}{List of posterior cluster assignments.}
#'          \item{`teams_out`}{List of refined cluster assignments.}
#'          \item{`tree_out`}{List of spanning trees for each partition.}
#'          \item{`k_out`}{Vector with the number of clusters in spatial partition.}
#'          \item{`j_out`}{Vector with the number of clusters in refined partition.}
#'          \item{`mu_out`}{List of posterior samples of `mu`.}
#'          \item{`mu_teams_out`}{List of posterior samples of means for refined partition.}
#'          \item{`sigmasq_y_out`}{List of posterior samples of noise variance.}
#'          \item{`marginal_likelihood_out`}{Vector of marginal likelihood values.}
#'         }
#' @export
DP.RST <- function(Y, graph0, init_val, hyperpar,
                   MCMC, BURNIN, THIN,
                   PT = TRUE, seed = 1234) {

  set.seed(seed) # Set seed for reproducibility

  n = vcount(graph0)  # Number of vertices (observations)
  p = ncol(Y)  # Number of response variables (features)

  ### Extract hyperparameters and initial values
  sigmasq_mu <- lambda_s <- nu <- M <- temp <- k_max <- j_max <- alpha <- NULL  # Declare variables
  list2env(hyperpar, envir = environment())  # Load hyperparameters

  mstgraph_lst <- mu <- cluster <- teams <- mu_teams <- sigmasq_y <- NULL # Declare variables
  list2env(init_val, envir = environment()) # Load initial values

  # Delete vertex names if they exist
  if('name' %in% names(vertex_attr(graph0))) {
    graph0 = delete_vertex_attr(graph0, 'name')
  }
  inc_mat = get.edgelist(graph0, names = F) # Get the edge list
  adj_list = lapply(as_adj_list(graph0), FUN = function(x) {x$vid}) # Get adjacency list
  adj_edge_list = lapply(as_adj_edge_list(graph0), FUN = function(x) {x$eid}) # Get adjacency edge list

  k = as.numeric(apply(cluster, 2, max))  # Number of clusters for the crude spatial partition
  j_teams = rep(max(teams), times = M) # Number of clusters for the refined partition

  moves_track = c()  # Tracking moves during MCMC

  csize = list()  # Crude spatial cluster sizes
  subgraphs = list()  # Subgraphs for each cluster
  eid_btw_mst = list()  # Edges between clusters in the MST

  # Loop over temperatures (M)
  for(m in 1:M) {
    cluster_m = cluster[, m]
    teams_m = teams[, m]
    csize[[m]] = Table(cluster_m)
    mstgraph_m = mstgraph_lst[[m]]

    # Split graph into subgraphs for each cluster
    clust_vid_m = split(1:n, cluster_m)
    subgraphs[[m]] = lapply(clust_vid_m, function(vids, mstgraph) {
      induced_subgraph(mstgraph, vids)
    }, mstgraph_m)

    # Find edges between clusters in the MST
    inc_mat_mst = get.edgelist(mstgraph_m, names = F) # Get the edge list for the MST
    c1_m = cluster_m[inc_mat_mst[, 1]] # Clusters on one side of each edge
    c2_m = cluster_m[inc_mat_mst[, 2]] # Clusters on the other side
    idx_btw = which(c1_m != c2_m) # Find edges between different clusters
    eid_btw_mst[[m]] = (E(mstgraph_m)$eid)[idx_btw] # Store the edges between clusters
  }

  # Determine if edges are within a cluster or between clusters (n*M matrix)
  edge_status = apply(cluster, 2, FUN = .getEdgeStatus, inc_mat)

  ################# MCMC ####################

  ## Prepare storage for MCMC results
  cluster_out = array(0, dim = c((MCMC-BURNIN)/THIN, n)) # To store cluster assignments
  teams_out = list() # To store refined partition assignments
  tree_out = list() # To store spanning trees
  mu_out = list() # To store posterior samples of mu
  mu_teams_out = list() # To store posterior samples of mu for teams
  sigmasq_y_out = list() # To store posterior samples of noise variance
  marginal_likelihood_out = numeric((MCMC-BURNIN)/THIN)  # To store log marginal likelihoods
  k_out = numeric((MCMC-BURNIN)/THIN) # To store the number of clusters for spatial partition
  j_teams_out = numeric((MCMC-BURNIN)/THIN) # To store the number of clusters for refined partition
  moves_track_out = numeric((MCMC-BURNIN)/THIN) # To track accepted moves

  # Matrix for storing marginal likelihood during MCMC iterations
  marginal_likelihood = matrix(0, nrow = MCMC, ncol = M)

  ## Main MCMC iteration loop
  for(iter in 1:MCMC) {
    for(m in 1:M) {

      # Update the Minimum Spanning Tree (MST)
      mstgraph = .proposeMST(graph0, edge_status[, m], subgraphs[[m]])
      mstgraph_lst[[m]] = mstgraph$mstgraph # Update the MST for temperature m
      eid_btw_mst[[m]] = mstgraph$eid_btw_mst  # Update edges between clusters in the MST
      subgraphs[[m]] = mstgraph$subgraphs  # Update the subgraphs

      k_m = k[m] # Number of clusters for spatial partition at temperature m
      j_teams_m = j_teams[m] # Number of clusters for refined partition at temperature m
      mstgraph_m = mstgraph_lst[[m]] # Updated spanning tree for temperature m
      edge_status_m = edge_status[, m] # Status of edges (within or between clusters)
      cluster_m = cluster[, m] # Cluster assignments for temperature m
      teams_m = teams[, m] # Refined partition assignments for temperature m
      csize_m = csize[[m]] # Sizes of the clusters in the spatial partition
      subgraphs_m = subgraphs[[m]] # Subgraphs for clusters in spatial partition
      eid_btw_mst_m = eid_btw_mst[[m]] # Edges between clusters in the MST
      sigmasq_y_m = sigmasq_y[[m]] # Noise variance for temperature m

      # Calculate the marginal likelihood for the current state
      marginal_likelihood[iter, m] <- .MarginalLikelihood(Y = Y, cluster_assign = cluster_m,
                                                         k = k_m, j = j_teams_m, team_assign = teams_m,
                                                         sigmasq_mu = sigmasq_mu,
                                                         Sigma = sigmasq_y_m)

      # Set move to "change move" (split/merge) only
      move = 3

      if(move == 3) { ## change move
        # Perform the death move (merge two clusters) (c1, c2) -> c2
        merge_res = .mergeCluster(mstgraph_m, eid_btw_mst_m, subgraphs_m, csize_m,
                                 cluster_m, inc_mat, change = T)
        # Perform the birth move (split a cluster)
        split_res = .splitCluster(mstgraph_m, k_m-1, merge_res$subgraphs, merge_res$csize)

        # Update cluster structure after merging
        update_res_check_merge = .updateMerge(merge_res, subgraphs_m, csize_m, eid_btw_mst_m,
                                             cluster_m, edge_status_m, adj_list, adj_edge_list, mstgraph_m)

        # Update cluster structure after splitting
        update_res_check_split = .updateSplit(split_res, update_res_check_merge$subgraphs, k_m-1, update_res_check_merge$csize,
                                             update_res_check_merge$eid_btw_mst, update_res_check_merge$cluster,
                                             update_res_check_merge$estatus, adj_list, adj_edge_list)

        # Calculate new marginal log-likelihood
        marginal_likelihood_new <- .MarginalLikelihood(Y = Y, team_assign = teams_m,
                                                      cluster_assign = update_res_check_split$cluster,
                                                      k = k_m, j = j_teams_m,
                                                      sigmasq_mu = sigmasq_mu,
                                                      Sigma = sigmasq_y_m)

        # Find marginal log-likelihood ratio
        log_L = marginal_likelihood_new - marginal_likelihood[iter, m]

        if (PT == TRUE) {
          # Log-likelihood ratio adjusted for the current temperature
          log_L = temp[m] * log_L
        }

        moves_track[iter] <- "Move 3 Rejected"

        # Calculate acceptance probability
        acc_prob = min(0, log_L)
        acc_prob = exp(acc_prob)
        if(runif(1) < acc_prob){
          # Accept the move if condition met
          moves_track[iter] <- "Move 3 Accepted"

          update_res = update_res_check_split

          subgraphs[[m]] = update_res$subgraphs # Update subgraphs
          csize[[m]] = update_res$csize # Update cluster sizes
          eid_btw_mst[[m]] = update_res$eid_btw_mst # Update edges between clusters
          cluster[, m] = update_res$cluster # Update cluster assignments
          edge_status[, m] = update_res$estatus # Update edge status
          marginal_likelihood[iter, m] = marginal_likelihood_new # Update marginal likelihood
        }
      }


      ### Perform DPM algoirithm
      for(ind in 1:k[m]){
        # Remove the k-th group from the groups means (clusters of observations)
        mu.subset = mu[[m]][-ind,]

        if (!is.matrix(mu.subset)) {
          mu.subset = matrix(mu.subset, ncol = p)
        }

        # The mean of the removed group (p-dimensional vector)
        mu.k = matrix(mu[[m]][ind,], ncol = p)
        # Remove the k-th group assignment from the team indicator
        team_assign.subset = teams_m[-ind]
        # We need to rank the observations again to check if by deleting the i^th indivdual any team gets empty or not
        team_assign.subset = dense_rank(team_assign.subset)
        # Update Z matrix
        Z.subset <- table(sequence(length(team_assign.subset)), team_assign.subset)

        # Check what is the maximum number of occupied tables
        J.current = max(team_assign.subset)

        # Initialize the storing values
        n_j = 0; prob.existing.table = c()

        group_var_inv = (1/sigmasq_mu) * ginv(sigmasq_y_m)

        # Calculate team specific statistics
        for(j in 1:J.current){
          # Calculate number of groups in each team
          n_j[j] = sum(team_assign.subset == j)

          mu.subset.j = mu.subset[which(Z.subset[, j] == 1), ]
          # Make sure that mu.subset.j is matrix
          if (!is.matrix(mu.subset.j)) {
            mu.subset.j = matrix(mu.subset.j, ncol = p)
          }

          # Calculate the terms to find the density value
          delta_star <- ginv(ginv(sigmasq_y_m) + n_j[j] * group_var_inv)
          theta_star <- delta_star %*% (n_j[j] * group_var_inv %*% colMeans(mu.subset.j))
          B_var = ginv(group_var_inv - group_var_inv %*% ginv(group_var_inv + ginv(delta_star)) %*% group_var_inv)
          beta_means <- B_var %*% (group_var_inv %*% ginv(group_var_inv + ginv(delta_star)) %*% ginv(delta_star) %*% theta_star)

          # Log probability for joining existing teams
          prob.existing.table[j] = log(n_j[j]) + dmvnorm(mu.k, mean = beta_means, sigma = B_var, log = TRUE)
        }

        B_var_new = ginv(group_var_inv - group_var_inv %*% ginv(group_var_inv + ginv(sigmasq_y_m)) %*% group_var_inv)
        # Log probability for starting a new team
        prob.new.table = log(alpha) + dmvnorm(mu.k, mean = rep(0, p), sigma = B_var_new, log = TRUE)

        # Calculate the unnormalized probabilities of occupying the existing J.current tables and opening a new table
        probability.unnormalized = c(prob.existing.table, prob.new.table)
        # Apply the log-sum exponential trick to stablize numerical computation
        probability.unnormalized = exp(probability.unnormalized  - max(probability.unnormalized))
        # Calculate the probability
        probability = probability.unnormalized/sum(probability.unnormalized)
        # Draw a categorcal randm variable with probability proportional to that calculated
        # The i^th table in the original assignment of the table is re-assigned the drawn value
        teams_m_cat = rcat(n = 1, prob = probability)
        teams_m = append(team_assign.subset, teams_m_cat, after = (ind-1))
        # We need to rank the observations again to check if by deleting the i^th indivdual any table gets empty or not
        teams_m = dense_rank(teams_m)
        j_teams[m] = max(teams_m)
      } # end of DPM loop

      # Update variables
      teams[, m] = teams_m
      k_m = k[m]
      cluster_m = cluster[, m]
      csize_m = csize[[m]]
      mu_teams_m = mu_teams[[m]]
      j_teams_m = j_teams[m]

      ## Update mu_teams (team means for refined partition)
      Z <- table(sequence(length(teams[, m])), teams[, m]) # Refined partition binary assignments matrix
      Delta = t(Z) %*% Z + sigmasq_mu * diag(1, j_teams_m)
      mu_teams[[m]] <- rmatnorm(M = ginv(Delta) %*% t(Z) %*% mu[[m]],
                                              U = ginv(Delta),
                                              V = sigmasq_mu*sigmasq_y_m)

      ## Update mu (group means for spatial partition)
      X <- table(sequence(length(cluster_m)), cluster_m)
      omega = sigmasq_mu * t(X) %*% X + diag(1, k_m)
      mu[[m]] <- rmatnorm(M = ginv(omega) %*% (sigmasq_mu*t(X) %*% Y + Z %*% mu_teams[[m]]),
                                        U = ginv(omega),
                                        V = sigmasq_mu*sigmasq_y_m)


      ## Update Sigma
      Y_hat = X %*% mu[[m]]
      mu_hat = Z %*% mu_teams[[m]]
      sigmasq_y[[m]] = riwish(v = nu + n + k_m + j_teams_m,
                                        S = lambda_s + t((Y - Y_hat)) %*% (Y - Y_hat) +
                                          1/sigmasq_mu * t(mu[[m]] - mu_hat) %*% (mu[[m]] - mu_hat) +
                                          t(mu_teams[[m]]) %*% mu_teams[[m]])

    } # end of M loop


    # Track the algorithm progress
    if(iter %% 100 == 0)
      cat('Iteration', iter, 'done\n')

    if (PT == TRUE) { # Swap temperatures

      # Check if the current iteration is a multiple of 10
      if (iter %% 10 == 0) {

        # Generate pairs of temperatures
        pairs <- generate_pairs(temp, diff_threshold = 0.1)

        # Compare the posterior distribution of the pairs of chains
        for (j in 1:length(pairs)) {

          # Get the values of the jth pair from the list
          idx_val_1 = pairs[[j]][1]; idx_val_2 = pairs[[j]][2]
          # Get the index in temp vector of the jth pair from the list
          idx_1 = which(temp == idx_val_1); idx_2 = which(temp == idx_val_2)

          # Calculate new and old log-marginal likelihoods for swapping temperatures
          log_marg_lik_new <- temp[idx_1] * marginal_likelihood[iter, idx_2] +
            temp[idx_2] * marginal_likelihood[iter, idx_1]
          log_marg_lik_old <- temp[idx_2] * marginal_likelihood[iter, idx_2] +
            temp[idx_1] * marginal_likelihood[iter, idx_1]

          # Compute acceptance probability
          mh_ratio <- log_marg_lik_new - log_marg_lik_old
          acc_prob = min(0, mh_ratio)
          acc_prob = exp(acc_prob)

          # If the condition is met, perform the swapping
          if (runif(1) < acc_prob) {

            # Fix the current values from one temperature
            subgraphs_temp = subgraphs[[idx_1]]
            csize_temp = csize[[idx_1]]
            eid_btw_mst_temp = eid_btw_mst[[idx_1]]
            cluster_temp = cluster[, idx_1]
            k_temp = k[idx_1]
            edge_status_temp = edge_status[, idx_1]
            mstgraph_lst_temp = mstgraph_lst[[idx_1]]
            marginal_likelihood_temp = marginal_likelihood[iter, idx_1]
            sigmasq_y_temp = sigmasq_y[[idx_1]]
            mu_temp = mu[[idx_1]]
            j_teams_temp = j_teams[idx_1]
            mu_teams_temp = mu_teams[[idx_1]]
            teams_temp = teams[, idx_1]

            # Perform the first change
            subgraphs[[idx_1]] = subgraphs[[idx_2]]
            csize[[idx_1]] = csize[[idx_2]]
            eid_btw_mst[[idx_1]] = eid_btw_mst[[idx_2]]
            cluster[, idx_1] = cluster[, idx_2]
            k[idx_1] = k[idx_2]
            edge_status[, idx_1] = edge_status[, idx_2]
            mstgraph_lst[[idx_1]] = mstgraph_lst[[idx_2]]
            marginal_likelihood[iter, idx_1] = marginal_likelihood[iter, idx_2]
            sigmasq_y[[idx_1]] = sigmasq_y[[idx_2]]
            mu[[idx_1]] = mu[[idx_2]]
            j_teams[idx_1] = j_teams[idx_2]
            mu_teams[[idx_1]] = mu_teams[[idx_2]]
            teams[, idx_1] = teams[, idx_2]

            # Assign the fixed temperatures to the second temperature
            subgraphs[[idx_2]] = subgraphs_temp
            csize[[idx_2]] = csize_temp
            eid_btw_mst[[idx_2]] = eid_btw_mst_temp
            cluster[, idx_2] = cluster_temp
            k[idx_2] = k_temp
            edge_status[, idx_2] = edge_status_temp
            mstgraph_lst[[idx_2]] = mstgraph_lst_temp
            marginal_likelihood[iter, idx_2] = marginal_likelihood_temp
            sigmasq_y[[idx_2]] = sigmasq_y_temp
            mu[[idx_2]] = mu_temp
            j_teams[idx_2] = j_teams_temp
            mu_teams[[idx_2]] = mu_teams_temp
            teams[, idx_2] = teams_temp
          }
        }
      } else {temp = temp}
    }

    ## Save results
    if(iter > BURNIN & (iter - BURNIN) %% THIN == 0) {
      tree_out[[(iter-BURNIN)/THIN]] = mstgraph_lst[[length(temp)]]
      cluster_out[(iter-BURNIN)/THIN, ] = cluster[, length(temp)]
      teams_out[[(iter-BURNIN)/THIN]] = teams[, length(temp)]
      k_out[(iter-BURNIN)/THIN] = k[length(temp)]
      j_teams_out[(iter-BURNIN)/THIN] = j_teams[length(temp)]
      marginal_likelihood_out[(iter-BURNIN)/THIN] = marginal_likelihood[iter, length(temp)]
      mu_out[[(iter-BURNIN)/THIN]] = mu[[length(temp)]]
      mu_teams_out[[(iter-BURNIN)/THIN]] = mu_teams[[length(temp)]]
      sigmasq_y_out[[(iter-BURNIN)/THIN]] = sigmasq_y[[length(temp)]]
      #moves_track_out[(iter-BURNIN)/THIN] = moves_track[iter]
    }

  } # end of iteration loop


  mode(cluster_out) = 'integer'  # to save memory

  # Return the final MCMC results
  return(list('cluster_out' = cluster_out, # Posterior samples of cluster assignments
              'teams_out' = teams_out, # Posterior samples of refined partition (team) assignments
              'tree_out' = tree_out, # Spanning trees for each temperature for each partition
              'k_out' = k_out, # Number of clusters in the spatial partition for each posterior sample
              'j_teams_out' = j_teams_out, # Number of refined clusters (teams) for each posterior sample
              'marginal_likelihood_out' = marginal_likelihood_out, # Marginal likelihood values
              'mu_out' = mu_out,# Posterior samples of the means (mu) for clusters in the spatial partition
              'mu_teams_out' = mu_teams_out, # Posterior samples of the means (mu) for clusters in the refined partition
              'sigmasq_y_out' = sigmasq_y_out)) # Posterior samples of the noise variance))
}
