#' Split an Existing Cluster Based on a Minimum Spanning Tree (MST)
#'
#' This function selects a cluster, removes an edge from its MST, and splits
#' it into two connected components.
#'
#' @param mstgraph An `igraph` object representing the full Minimum Spanning Tree (MST).
#' @param k Integer, the total number of spatial clusters.
#' @param subgraphs A list of `igraph` subgraphs, each representing a cluster.
#' @param csize An integer vector where each element represents the size of a cluster.
#'
#' @returns A list containing:
#' \item{vid_old}{A vector of vertex IDs remaining in the original cluster.}
#' \item{vid_new}{A vector of vertex IDs forming the new cluster after the split.}
#' \item{eid_cutted}{The ID of the edge that was removed.}
#' \item{clust_old}{The ID of the cluster that was split.}
#' \item{idx_new}{A logical vector indicating which nodes belong to the new cluster.}
#'
#' @noRd
.splitCluster <- function(mstgraph, k, subgraphs, csize) {
  clust_split = sample.int(k, 1, prob = csize - 1)
  mst_subgraph = subgraphs[[clust_split]]

  while(length(E(mst_subgraph)) == 0) {
    clust_split = sample.int(k, 1, prob = csize - 1)
    mst_subgraph = subgraphs[[clust_split]]
    # print("No edges in subgraph")
  }

  intc = sample.int(csize[clust_split]-1, 1)
  while (intc > length(E(mst_subgraph))) {
    intc = sample.int(csize[clust_split]-1, 1)
    # print("More number of observations than edges")
  }

  edge_cutted = E(mst_subgraph)[intc]
  # edge_cutted = E(mst_subgraph)[sample.int(csize[clust_split]-1, 1)]
  eid_cutted = edge_cutted$eid
  mst_subgraph = delete.edges(mst_subgraph, edge_cutted)
  connect_comp = components(mst_subgraph)
  idx_new = (connect_comp$membership == 2)
  vid_new = V(mst_subgraph)$vid[idx_new]
  vid_old = V(mst_subgraph)$vid[!idx_new]

  return(list(vid_old = vid_old, vid_new = vid_new, eid_cutted = eid_cutted,
              clust_old = clust_split, idx_new = idx_new))
}
