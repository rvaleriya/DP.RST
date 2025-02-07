#' Generate a Constrained Delaunay Triangulation Graph
#'
#' This function constructs a graph representation of a constrained Delaunay triangulation from a given mesh.
#' It removes edges that connect boundary nodes and applies a length threshold to filter long edges.
#'
#' @param n An integer specifying the number of interior nodes to consider.
#' @param mesh A mesh object containing node coordinates and edge information.
#' @param threshold A numeric value specifying the maximum edge length allowed in the graph. Defaults to 5000.
#'
#' @returns An undirected graph (`igraph` object) where nodes represent mesh points and edges represent valid triangulation connections.
#' @export
#'
#' @examples
constrainedDentri <- function(n, mesh, threshold = 5000) {
  coords = mesh$nodes[1:n, ]

  # drop edges that connect boundary nodes
  rid_drop = mesh$bnd_edges
  edge_list = mesh$edges[!rid_drop, ]

  # compute edge length
  distance = sqrt( rowSums((coords[edge_list[, 1], ] - coords[edge_list[, 2], ]) ^ 2) )

  rid_drop = distance > threshold
  edge_list = edge_list[!rid_drop, ]
  distance = distance[!rid_drop]

  graph0 = graph_from_edgelist(edge_list, directed = F)
  E(graph0)$weight = distance

  return(graph0)
}
