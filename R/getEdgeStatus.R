#' Determine Edge Status Based on Cluster Membership
#'
#' This function assigns a status to each edge in a graph based on whether it
#' connects nodes within the same cluster ('w') or between two clusters ('b').
#'
#' @param membership An integer vector where each element represents the cluster assignment of a node.
#' @param inc_mat A numeric matrix (E x 2), where each row represents an edge
#'        with two columns indicating the indices of the connected nodes.
#'
#' @returns A character vector of length `E`, where each element is either:
#' \item{'w'}{If the edge is **within** a cluster.}
#' \item{'b'}{If the edge is **between** two clusters.}
#'
#' @keywords internal
#'
#' @noRd
.getEdgeStatus <- function(membership, inc_mat) {
  membership_head = membership[inc_mat[, 1]]
  membership_tail = membership[inc_mat[, 2]]
  edge_status = rep('w', nrow(inc_mat))
  edge_status[membership_head != membership_tail] = 'b'
  return(edge_status)
}
