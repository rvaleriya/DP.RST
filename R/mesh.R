#' Generate a 2D Mesh
#'
#' This function generates a 2D triangular mesh from given interior and boundary coordinates.
#'
#' @param coords A numeric matrix with two columns (x, y) representing the coordinates of interior nodes.
#' @param bnd A list with elements `x` and `y`, representing the boundary node coordinates.
#' @param ... Additional arguments passed to `fdaPDE::create.mesh.2D`.
#'
#' @return A mesh object created using `fdaPDE::create.mesh.2D`, with two additional attributes:
#' \itemize{
#'   \item \code{n_int}: Number of interior nodes.
#'   \item \code{bnd_edges}: A logical vector indicating which edges connect boundary nodes.
#' }
#' @export
#'
#' @examples
#' coords <- matrix(runif(20), ncol=2)
#' bnd <- list(x = c(0,1,1,0,0), y = c(0,0,1,1,0))
#' mesh <- gen2dMesh(coords, bnd)
#' print(mesh)
gen2dMesh <- function(coords, bnd, ...) {
  # note the first and last boundary nodes are the same
  n = nrow(coords); n_bnd = length(bnd$x) - 1
  coords_all = rbind(coords, cbind(x = bnd$x, y = bnd$y)[1:n_bnd, ])

  # get boundary segments
  segments = cbind( (n+1):(n+n_bnd), c((n+2):(n+n_bnd), n+1) )

  mesh = create.mesh.2D(coords_all, segments = segments, ...)
  mesh$n_int = n  # number of interior nodes
  # edges that connect boundary nodes
  mesh$bnd_edges = apply(mesh$edges, 1, FUN = function(x) any(x > n))

  return(mesh)
}
