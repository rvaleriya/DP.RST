#' Merge Two Existing Clusters in an MST
#'
#' This function selects an edge from `eid_btw_mst` to merge two clusters
#' and updates the cluster assignments accordingly.
#'
#' @param mstgraph An `igraph` object representing the minimum spanning tree (MST).
#' @param eid_btw_mst An integer vector of edge IDs that are between clusters in the MST.
#' @param subgraphs A list of `igraph` subgraphs representing current clusters.
#' @param csize An integer vector where each element represents the size of a cluster.
#' @param cluster An integer vector mapping each node to its cluster.
#' @param edge_list A matrix where each row represents an edge and its corresponding nodes.
#' @param change Logical (`TRUE` or `FALSE`), whether to update `subgraphs` and `csize` directly.
#'
#' @returns A list with:
#' \item{vid_old}{A vector of vertex IDs from the merged cluster (`c2`).}
#' \item{vid_new}{A vector of vertex IDs in the newly formed cluster (`c1`).}
#' \item{clust_old}{The ID of the merged (removed) cluster.}
#' \item{clust_new}{The ID of the new merged cluster.}
#' \item{edge_merge}{The index of the selected merging edge.}
#' \item{subgraphs}{Updated list of subgraphs (if `change = TRUE`).}
#' \item{csize}{Updated cluster sizes (if `change = TRUE`).}
#'
#' @noRd
.mergeCluster <- function(mstgraph, eid_btw_mst, subgraphs, csize, cluster,
                          edge_list, change = F) {
  # edge for merging
  edge_merge = sample.int(length(eid_btw_mst), 1)
  # update cluster information
  # clusters of endpoints of edge_merge
  eid_merge = eid_btw_mst[edge_merge]
  clusters_merge = cluster[edge_list[eid_merge, ]]
  clusters_merge = sort(clusters_merge)
  c1 = clusters_merge[1]; c2 = clusters_merge[2] # note c1 < c2
  # merge c2 to c1

  # vid of vertices in c2
  vid_old = V(subgraphs[[c2]])$vid
  # vid in merged cluster
  vid_new = c(V(subgraphs[[c1]])$vid, vid_old)

  csize_new = NULL; subgraphs_new = NULL
  if(change) {
    subgraphs_new = subgraphs
    subgraphs_new[[c1]] = induced_subgraph(mstgraph, vid_new)
    subgraphs_new[[c2]] = NULL

    csize_new = csize
    csize_new[c1] = length(vid_new)
    csize_new = csize_new[-c2]
  }

  # now drop c2
  return(list(vid_old = vid_old, vid_new = vid_new, clust_old = c2, clust_new = c1,
              edge_merge = edge_merge, subgraphs = subgraphs_new, csize = csize_new))
}
