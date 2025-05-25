#' Simple internal assertion helper
#' @keywords internal
#' @noRd
.assert <- function(cond, msg) {
  if (!cond) stop(msg, call. = FALSE)
  invisible(TRUE)
}
