#' Update Cluster Structure After a Split Move
#'
#' This function updates the cluster structure, edge status, and adjacency information
#' after a split move is accepted.
#'
#' @param split_res A list containing information from `splitCluster()`, including split clusters and edges.
#' @param subgraphs A list of `igraph` subgraphs representing current clusters.
#' @param k Integer, the total number of spatial clusters before the split.
#' @param csize An integer vector where each element represents the size of a cluster.
#' @param eid_btw_mst An integer vector of edge IDs that were between clusters in the MST.
#' @param cluster An integer vector mapping each node to its cluster.
#' @param edge_status A character vector indicating edge statuses ('w' for within cluster, 'b' for between clusters).
#' @param adj_list A list where each element is a vector of adjacent vertex IDs.
#' @param adj_edge_list A list where each element is a vector of corresponding edge IDs.
#'
#' @returns A list with:
#' \item{subgraphs}{Updated list of subgraphs after splitting.}
#' \item{csize}{Updated cluster sizes.}
#' \item{cluster}{Updated cluster assignments for each node.}
#' \item{eid_btw_mst}{Updated edge list for edges between clusters.}
#' \item{estatus}{Updated edge status vector.}
#'
#' @keywords internal
#'
#' @noRd
.updateSplit <- function(split_res, subgraphs, k, csize, eid_btw_mst, cluster,
                         edge_status, adj_list, adj_edge_list) {
  clust_split = split_res$clust_old
  vid_old = split_res$vid_old; vid_new = split_res$vid_new

  subgraph_split = subgraphs[[clust_split]]
  idx_new = split_res$idx_new
  subgraphs[[clust_split]] = induced_subgraph(subgraph_split, !idx_new)  # subgraph of old cluster
  subgraphs[[k+1]] = induced_subgraph(subgraph_split, idx_new) # subgraph of new cluster

  csize[clust_split] = length(vid_old)
  csize[k+1] = length(vid_new)

  cluster[vid_new] = k + 1
  eid_btw_mst = c(eid_btw_mst, split_res$eid_cutted)

  # update edge status
  adj_vid_old = unlist(adj_list[vid_old])
  adj_eid_old = unlist(adj_edge_list[vid_old])
  clust_adj_old = cluster[adj_vid_old]
  idx_btw = which(clust_adj_old != clust_split)
  eid_btw = adj_eid_old[idx_btw]
  edge_status[eid_btw] = 'b'

  return(list(subgraphs = subgraphs, csize = csize, cluster = cluster,
              eid_btw_mst = eid_btw_mst, estatus = edge_status))
}
