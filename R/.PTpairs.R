#' Generate Pairs of Temperatures for Swapping on a Difference Threshold
#'
#' This function randomly selects pairs of values from a vector where
#' the absolute difference between paired values is within a given threshold.
#'
#' @param vec A numeric vector of values from which pairs will be generated.
#' @param diff_threshold A numeric value specifying the maximum absolute
#'        difference allowed between paired values. Default is `0.1`.
#'
#' @returns A list where each element is a numeric vector of length 2,
#'          representing a valid pair. If no valid pairs exist, returns an empty list.
#'
#' @keywords internal
#'
#' @noRd
generate_pairs <- function(vec, diff_threshold = 0.1) {
  # If there are less than 2 elements, no pairs can be made
  if (length(vec) < 2) stop("Vector length must be at least 2.")

  # Shuffle the vector
  vec <- sample(vec)

  # Initialize an empty list to store the pairs
  pairs <- list()

  while(length(vec) >= 2) {
    # Get the first element
    first_elem <- vec[1]

    # Find elements within the diff_threshold
    valid_elems <- vec[abs(vec - first_elem) <= diff_threshold]

    # If there are valid elements
    if(length(valid_elems) > 1) {
      # Choose a second element
      second_elem <- sample(valid_elems[valid_elems != first_elem], 1)  # Remove the first element before sampling
      pairs[[length(pairs) + 1]] <- c(first_elem, second_elem)

      # Remove these elements from vec
      vec <- setdiff(vec, c(first_elem, second_elem))
    } else {
      # If there are no valid elements, remove the current element from vec
      vec <- vec[-1]
    }
  }

  return(pairs)
}
