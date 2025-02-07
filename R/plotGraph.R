#' Plot a Graph Using ggplot2
#'
#' This function visualizes a graph by plotting its edges based on provided coordinates.
#'
#' @param coords A matrix or data frame with two columns representing node coordinates (x, y).
#' @param graph An `igraph` object representing the graph structure.
#' @param title A character string specifying the plot title (optional).
#'
#' @returns A `ggplot2` object visualizing the graph edges.
#' @export
#'
#' @examples
#' # Generate a mesh
#' coords <- matrix(runif(20), ncol = 2)
#' bnd <- list(x = c(0,1,1,0,0), y = c(0,0,1,1,0))
#' mesh <- gen2dMesh(coords, bnd)
#'
#' # Create a constrained Delaunay triangulation graph
#' graph <- constrainedDentri(n = mesh$n_int, mesh = mesh, threshold = 0.5)
#'
#' # Plot the graph
#' plotGraph(coords, graph, title = "Graph from Mesh")
#' # Plot the graph with the boundary
#' spatial_graph <- plotGraph(coords, graph, title = "Graph with Boundary from Mesh") +
#'                   geom_boundary(bnd)
#' spatial_graph
plotGraph <- function(coords, graph, title = NULL){
  edgelist = get.edgelist(graph)
  edgedata = data.frame(coords[edgelist[,1 ], ], coords[edgelist[, 2], ])
  colnames(edgedata) = c("x1", "y1", "x2", "y2")

  ggplot() + geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data = edgedata,
                          size = 0.5, colour = "grey") +
    labs(title = title, x = "lon", y = "lat")+
    theme(plot.title = element_text(hjust = 0.5))
}
