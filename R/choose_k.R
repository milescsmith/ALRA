#' @title choose_k
#'
#' @description Heuristic for choosing rank k for the low rank approximation based on
#  statistics of the spacings between consecutive singular values. Finds
#  the smallest singular value \sigma_i such that $\sigma_i - \sigma_{i-1}
#  is significantly different than spacings in the tail of the singular values.
#'
#' @param A_norm The log-transformed expression matrix of cells (rows) vs. genes (columns)
#' @param K Number of singular values to compute. Must be less than the smallest dimension of the matrix.
#' @param pval_thresh The threshold for ``significance''
#' @param noise_start Index for which all smaller singular values are considered noise
#' @param q Number of additional power iterations
#'
#' @importFrom rsvd rsvd
#'
#' @return A list with three items
#       1) Chosen k
#       2) P values of each possible k
#       3) Singular values of the matrix A_norm
#' @export
#'
#' @examples
choose_k <- function(A_norm, 
                     K = 100, 
                     pval_thresh = 1E-10, 
                     noise_start = 80, 
                     q = 2) {
  
  if (K > min(dim(A_norm))) {
    stop("For an m by n matrix, K must be smaller than the min(m,n).\n")
  }
  if (noise_start > K - 5) {
    stop("There need to be at least 5 singular values considered noise.\n")
  }
  noise_svals <- noise_start:K
  rsvd_out <- rsvd(A = A_norm, 
                   k = K, 
                   q = q)
  diffs <- diff(rsvd_out$d)
  pvals <- pnorm(diffs, 
                 mean(diffs[noise_svals - 1]), 
                 sd(diffs[noise_svals - 1]))
  k <- which(pvals < pval_thresh) %>% max()
  return(list(k = k, pvals = pvals, d = rsvd_out$d))
}
