#' @title normalize_data
#'
#' @description Simple convenience function to library and log normalize a matrix
#'
#' @param A matrix to be normalized
#'
#' @return a normalized matrix
#' @export
#'
#' @examples
normalize_data <- function(A) {
  #  

  totalUMIPerCell <- rowSums(A)
  if (any(totalUMIPerCell == 0)) {
    toRemove <- which(totalUMIPerCell == 0)
    A <- A[-toRemove, ]
    totalUMIPerCell <- totalUMIPerCell[-toRemove]
    message(sprintf("Removed {length(toRemove)} cells which did not express any genes\n"))
  }

  A_norm <- sweep(A, 1, totalUMIPerCell, "/")
  A_norm <- A_norm * 10E3
  A_norm <- log(A_norm + 1)
}
