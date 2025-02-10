#' Update Cluster Structure After a Merge Move
#'
#' This function updates the cluster structure, edge status, and adjacency information
#' after a merge move is accepted.
#'
#' @param res_merge A list containing information from `mergeCluster()`, including merged clusters and edges.
#' @param subgraphs A list of `igraph` subgraphs representing current clusters.
#' @param csize An integer vector where each element represents the size of a cluster.
#' @param eid_btw_mst An integer vector of edge IDs that were between clusters in the MST.
#' @param cluster An integer vector mapping each node to its cluster.
#' @param edge_status A character vector indicating edge statuses ('w' for within cluster, 'b' for between clusters).
#' @param adj_list A list where each element is a vector of adjacent vertex IDs.
#' @param adj_edge_list A list where each element is a vector of corresponding edge IDs.
#' @param mstgraph An `igraph` object representing the minimum spanning tree.
#'
#' @returns A list with:
#' \item{subgraphs}{Updated list of subgraphs after merging.}
#' \item{csize}{Updated cluster sizes.}
#' \item{cluster}{Updated cluster assignments for each node.}
#' \item{eid_btw_mst}{Updated edge list for edges between clusters.}
#' \item{estatus}{Updated edge status vector.}
#'
#' @noRd
.updateMerge <- function(res_merge, subgraphs, csize, eid_btw_mst, cluster,
                        edge_status, adj_list, adj_edge_list, mstgraph) {
  clust_old = res_merge$clust_old; clust_new = res_merge$clust_new
  vid_old = V(subgraphs[[clust_old]])$vid
  vid_new = c(V(subgraphs[[clust_new]])$vid, vid_old)
  subgraphs[[clust_new]] = induced_subgraph(mstgraph, vid_new)
  subgraphs[[clust_old]] = NULL

  csize[clust_new] = length(vid_new)
  csize = csize[-clust_old]

  cluster[vid_old] = clust_new
  idx = which(cluster > clust_old)
  cluster[idx] = cluster[idx] - 1

  eid_btw_mst = eid_btw_mst[-res_merge$edge_merge]

  # update edge status
  adj_vid_old = unlist(adj_list[vid_old])
  adj_eid_old = unlist(adj_edge_list[vid_old])
  clust_adj_old = cluster[adj_vid_old]
  idx_within = which(clust_adj_old == clust_new)
  eid_within = adj_eid_old[idx_within]
  edge_status[eid_within] = 'w'

  return(list(subgraphs = subgraphs, csize = csize, cluster = cluster,
              eid_btw_mst = eid_btw_mst, estatus = edge_status))
}
