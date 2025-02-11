#' Propose a New Minimum Spanning Tree (MST)
#'
#' This internal function generates a new MST by assigning weights to edges
#' based on their status ('w' for within cluster, 'b' for between clusters)
#' and then computing the MST of the updated graph.
#'
#' @param graph0 An `igraph` object representing the input graph.
#' @param edge_status A character vector indicating edge statuses ('w' for weak, 'b' for strong).
#' @param subgraphs A list of `igraph` subgraphs extracted from `graph0`.
#'
#' @returns A list with three elements:
#' \item{mstgraph}{The computed Minimum Spanning Tree (`igraph` object).}
#' \item{subgraphs}{A list of updated subgraphs derived from the new MST.}
#' \item{eid_btw_mst}{A vector of edge IDs in the MST with strong "between' clusters edges.}
#'
#' @keywords internal
#'
#' @noRd
.proposeMST <- function(graph0, edge_status, subgraphs) {
  nedge = length(edge_status)
  nb = sum(edge_status == 'b')
  nw = nedge - nb
  weight = numeric(nedge)
  weight[edge_status == 'w'] = runif(nw)
  weight[edge_status == 'b'] = runif(nb, 10, 20)
  E(graph0)$weight = weight
  mstgraph = mst(graph0, weights = E(graph0)$weight)

  # update subgraphs
  subgraphs_new = lapply(subgraphs, function(g, mstgraph) {induced_subgraph(mstgraph, V(g)$vid)},
                         mstgraph)
  # update eid_btw_mst
  eid_btw_mst = E(mstgraph)$eid[E(mstgraph)$weight >= 10]

  return(list(mstgraph = mstgraph, subgraphs = subgraphs_new, eid_btw_mst = eid_btw_mst))
}
