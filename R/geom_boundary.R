#' Create a ggplot2 segment layer for a complex boundary
#'
#' @param bnd A list with elements `x` and `y`, representing the boundary coordinates.
#' @param ... Additional arguments passed to `ggplot2::geom_segment()`.
#'
#' @returns A `ggplot2` `geom_segment` layer representing the boundary.
#' @export
geom_boundary <- function(bnd, ...) {
  n = length(bnd$x)
  segments = data.frame(
    x1 = bnd$x[-n], y1 = bnd$y[-n],
    x2 = bnd$x[-1], y2 = bnd$y[-1]
  )

  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = segments, ...)
}
