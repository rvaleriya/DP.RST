#' Plot Spatial Clusters
#'
#' Generates a scatter plot of spatial clusters, displaying the assigned groups
#' within the spatial boundary.
#'
#' @param coords A matrix or dataframe with two columns representing X and Y coordinates.
#' @param group_assign A vector of group assignments corresponding to `coords`.
#' @param bound A spatial boundary object to overlay on the plot.
#' @param point_size The size of the points in the plot (default: 2).
#'
#' @returns A ggplot object showing the spatial clusters.
#' @export
groups_plot <- function(coords, group_assign, bound, point_size = 2) {

  # Create a dataframe of coordinates and groups' assignment
  group_data <- as.data.frame(cbind(coords, group_assign))
  group_data$group_assign <- as.factor(group_data$group_assign)

  # colourCount = length(unique(group_assign))
  # getPalette = colorRampPalette(brewer.pal(8, "Accent"))
  # colors <- viridis::viridis_pal(option = "magma")(colourCount)

  # Make a plot
  group_pred <- ggplot() +
    geom_boundary(bound) +
    # geom_point(aes(x = group_data[,1], y = group_data[,2], colour = group_assign),
    #            data = group_data, size = point_size) +
    geom_point(aes(x = group_data[,1], y = group_data[,2], colour = as.factor(group_data[,3])),
               data = group_data, size = point_size) +
    # scale_colour_manual(values = colors, name = 'Groups', guide = guide_legend(ncol = 2)) +
    # scale_colour_manual(values = getPalette(colourCount), name = 'groups', guide = guide_legend(ncol = 2)) +
    labs(x = 'Coordinate X', y = 'Coordinate Y') +
    ggtitle('Groups of Observations')

  return(group_pred)
}
