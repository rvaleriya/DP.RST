#' Plot Refined Partition (Teams)
#'
#' Generates a scatter plot of the refined partition, displaying assigned teams
#' within the spatial boundary.
#'
#' @param coords A matrix or dataframe with two columns representing X and Y coordinates.
#' @param team_assign A vector of team assignments corresponding to `coords`.
#' @param group_assign A vector of group assignments corresponding to `coords`.
#' @param bound A spatial boundary object to overlay on the plot.
#' @param point_size The size of the points in the plot (default: 2).
#'
#' @returns A ggplot object showing the refined partition (teams).
#' @export
teams_plot <- function(coords, team_assign, group_assign, bound, point_size = 2) {

  X = table(sequence(length(group_assign)), group_assign)
  Z = table(sequence(length(team_assign)), team_assign)

  obs_in_teams = X %*% Z
  obs_teams_assign = obs_in_teams %*% sort(unique(team_assign))

  # Create a dataframe of coordinates and teams' assignment
  team_data <- as.data.frame(cbind(coords, obs_teams_assign))
  team_data$V3 <- as.factor(team_data$V3)

  # colourCount = length(unique(team_assign))
  # getPalette = colorRampPalette(brewer.pal(8, "Accent"))

  # Make a plot
  team_pred <- ggplot() +
    geom_boundary(bound) +
    geom_point(aes(x = team_data[,1], y = team_data[,2], colour = as.factor(team_data[,3])),
               data = team_data, size = point_size) +
    # geom_point(aes(x = team_data[,1], y = team_data[,2], colour = team_data[,3]),
    #            data = team_data, size = point_size) +
    # scale_colour_manual(values = getPalette(colourCount), name = 'teams') +
    labs(x = 'Coordinate X', y = 'Coordinate Y') +
    ggtitle('Teams of Observations')
  #+ annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, color = "red", fill = NA, size = 1.3)

  return(team_pred)
}
