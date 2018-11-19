#' sd_nonzero
#'
#' @param x 
#'
#' @return a list of the standard deviation of non-zero values
#' @export
#'
#' @examples
sd_nonzero <- function(x) {
  sd(x[!x == 0])
}
