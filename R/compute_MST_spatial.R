#' Compute MST based on the initial clustering assignment
#'
#' #' This function computes a minimum spanning tree (MST) for a set of points with
#' an initial clustering assignment. It penalizes edges that connect points from different
#' clusters, computes the MST on the penalized graph, and then identifies spatial components
#' (i.e., disconnected subgraphs) within each cluster. The output is a list containing the MST
#' (as an igraph object) and the input data frame augmented with a new column, \code{spatial_cluster},
#' that gives the spatial component membership for each point.
#'
#' @param data A data frame containing at least coordinate columns and a clustering column.
#' @param coords_cols Character vector of length 2 specifying the coordinate column names (default: \code{c("x", "y")}).
#' @param cluster_col Character string specifying the name of the column with the initial cluster membership (default: \code{"cluster"}).
#' @param bnd A list with elements \code{x} and \code{y} representing the boundary. If \code{NULL} (default), the convex hull of the points is used.
#' @param threshold Numeric value for the maximum allowable edge length in \code{constrainedDentri} (default: 5000).
#' @param penalty Numeric penalty added to inter-cluster edges. If \code{NULL} (default), it is set to \code{max(edge weight) + 1}.
#'
#' @returns A list with two elements:
#' \item{mst}{The MST as an igraph object (with all points connected).}
#' \item{data}{The original data frame augmented with a new column \code{spatial_cluster} containing the spatial component membership for each point.}
#'
#' @export
compute_MST_spatial <- function(data,
                                coords_cols = c("x", "y"),
                                cluster_col = "cluster",
                                bnd = NULL,
                                threshold = 5000,
                                penalty = NULL) {

  # Convert coordinates to a matrix
  coords <- as.matrix(data[, coords_cols])

  # If no boundary is supplied, compute the convex hull and form a boundary
  # if (is.null(bnd)) {
  #   hull_idx <- chull(coords)
  #   # Ensure boundary is closed by repeating the first point at the end
  #   bnd <- list(x = coords[c(hull_idx, hull_idx[1]), 1],
  #               y = coords[c(hull_idx, hull_idx[1]), 2])
  # }

  # Generate a mesh (requires your gen2dMesh function)
  mesh <- gen2dMesh(coords, bnd)

  # Create the graph using constrainedDentri (assumes function available)
  graph0 <- constrainedDentri(n = nrow(coords), mesh = mesh, threshold = threshold)
  E(graph0)$eid = c(1:ecount(graph0))  # edge id
  V(graph0)$vid = c(1:vcount(graph0))  # vertex id

  # Assign the supplied clustering as a vertex attribute
  V(graph0)$cluster <- data[[cluster_col]]

  # Get the original edge weights (computed in constrainedDentri)
  orig_weights <- E(graph0)$weight

  # If penalty is not provided, set it to max weight + 1
  if (is.null(penalty)) {
    penalty <- max(orig_weights, na.rm = TRUE) + 1
  }

  # Extract edge endpoints and their cluster memberships
  edge_ends <- as.data.frame(get.edgelist(graph0))
  colnames(edge_ends) <- c("V1", "V2")
  edge_ends$V1 <- as.numeric(edge_ends$V1)
  edge_ends$V2 <- as.numeric(edge_ends$V2)

  clusters_v1 <- V(graph0)$cluster[edge_ends$V1]
  clusters_v2 <- V(graph0)$cluster[edge_ends$V2]

  # Identify edges connecting different clusters
  inter_cluster <- clusters_v1 != clusters_v2

  # Adjust edge weights: add penalty for inter-cluster edges
  new_weights <- orig_weights + penalty * inter_cluster
  E(graph0)$new_weight <- new_weights

  # Compute the MST using the penalized weights
  mst_graph <- mst(graph0, weights = E(graph0)$new_weight)

  # Identify spatial components:
  # Remove MST edges that connect points with different initial clusters.
  head_idx <- as.numeric(head_of(mst_graph, E(mst_graph)))
  tail_idx <- as.numeric(tail_of(mst_graph, E(mst_graph)))
  diff_clusters <- V(mst_graph)$cluster[head_idx] != V(mst_graph)$cluster[tail_idx]

  mst_same_label <- delete_edges(mst_graph, which(diff_clusters))

  # Compute connected components (each component is a "spatial cluster")
  comp <- igraph::components(mst_same_label)

  # Append the spatial cluster membership to the data frame
  data$spatial_cluster <- comp$membership

  # Return a list containing the MST and the updated data frame
  return(list(mst = mst_graph, data = data))
}
