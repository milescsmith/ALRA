#' @title alra
#' 
#' @description Computes the k-rank approximation to A_norm and adjusts 
#' it according to the error distribution learned from the negative values.
#'
#' @param A_norm The log-transformed expression matrix of cells (rows) vs. genes (columns)
#' @param k the rank of the rank-k approximation. Set to 0 for automated choice of k.
#' @param q the number of additional power iterations in randomized SVD
#'
#' @import magrittr
#' @importFrom future plan multiprocess
#' @importFrom future.apply future_apply
#' @importFrom glue glue
#' @importFrom rsvd rsvd
#' @importFrom Matrix crossprod tcrossprod
#'
#' @return  A list with three items
#       1) The rank k approximation of A_norm.
#       2) The rank k approximation of A_norm, adaptively thresholded
#       3) The rank k approximation of A_norm, adaptively thresholded and
#       with the first two moments of the non-zero values matched to the
#       first two moments of the non-zeros of A_norm. This is the completed
#       matrix most people will want to work with
#' @export 
#' @examples result.completed <- adjusted_svd(A_norm,15)
#'     A_norm_rank15 <- result.completed[[1]]     # The low rank approximation for reference purposes...not suggested for matrix completion
#'     A_norm_rank15_cor <- result.completed[[3]] # The actual adjusted, completed matrix
alra <- function(A_norm, k = 0, q = 10, do.par = FALSE) {

  if (do.par) {
    plan(multiprocess)
    apply_fun <- future_apply
  } else {
    apply_fun <- apply
  }
  
  message(glue("Read matrix with {nrow(A_norm)} cells and {ncol(A_norm)} genes\n"))
  if (class(A_norm) != "matrix") {
    stop(glue("A_norm is of class {class(A_norm)}, but it should be of class matrix. Did you forget to run as.matrix()?"))
  }

  if (k == 0) {
    k_choice <- choose_k(A_norm)
    k <- k_choice$k
    message(glue("Chose k={k}\n"))
  }

  message("Getting nonzeros\n")
  originally_nonzero <- A_norm > 0

  message("Randomized SVD\n")
  fastDecomp_noc <- rsvd(A_norm, k, q = q)
  A_norm_rank_k <- fastDecomp_noc$u[, 1:k] %*% diag(fastDecomp_noc$d[1:k]) %>% tcrossprod(fastDecomp_noc$v[, 1:k])
  
  message("Find mins\n")
  A_norm_rank_k_mins <- apply_fun(A_norm_rank_k, 2, min) %>% abs()
  message("Sweep\n")
  A_norm_rank_k_cor <- replace(A_norm_rank_k, 
                               A_norm_rank_k <= A_norm_rank_k_mins[col(A_norm_rank_k)], 
                               0)

  sigma_1 <- apply_fun(A_norm_rank_k_cor, 2, sd_nonzero)
  sigma_2 <- apply_fun(A_norm, 2, sd_nonzero)
  mu_1 <- colSums(A_norm_rank_k_cor) / colSums(!!A_norm_rank_k_cor)
  mu_2 <- colSums(A_norm) / colSums(!!A_norm)

  toscale <- !is.na(sigma_1) & !is.na(sigma_2)

  message(glue("Scaling all except for {sum(!toscale)} columns\n"))

  sigma_1_2 <- sigma_2 / sigma_1
  toadd <- -1 * mu_1 * sigma_2 / sigma_1 + mu_2

  A_norm_rank_k_temp <- A_norm_rank_k_cor[, toscale]
  A_norm_rank_k_temp <- sweep(A_norm_rank_k_temp, 2, sigma_1_2[toscale], FUN = "*")
  A_norm_rank_k_temp <- sweep(A_norm_rank_k_temp, 2, toadd[toscale], FUN = "+")

  A_norm_rank_k_cor_sc <- A_norm_rank_k_cor
  A_norm_rank_k_cor_sc[, toscale] <- A_norm_rank_k_temp
  A_norm_rank_k_cor_sc[A_norm_rank_k_cor == 0] <- 0

  lt0 <- A_norm_rank_k_cor_sc < 0
  A_norm_rank_k_cor_sc[lt0] <- 0
  num_of_zeros <- 100 * sum(lt0) / (nrow(A_norm) * ncol(A_norm))
  message(glue("{format(num_of_zeros, digits = 2)} of the values became negative in the scaling process and were set to zero\n"))

  A_norm_rank_k_cor_sc[originally_nonzero & A_norm_rank_k_cor_sc == 0] <- A_norm[originally_nonzero & A_norm_rank_k_cor_sc == 0]

  colnames(A_norm_rank_k_cor) <- colnames(A_norm)
  colnames(A_norm_rank_k_cor_sc) <- colnames(A_norm)
  colnames(A_norm_rank_k) <- colnames(A_norm)

  original_nz <- sum(A_norm > 0) / (nrow(A_norm) * ncol(A_norm))
  completed_nz <- sum(A_norm_rank_k_cor_sc > 0) / (nrow(A_norm) * ncol(A_norm))
  cat(glue("The matrix went from {format((100 * original_nz), digits = 2)} nonzero to {format((100 * completed_nz), digits = 2)} nonzero\n"))

  list(A_norm_rank_k = A_norm_rank_k, A_norm_rank_k_cor = A_norm_rank_k_cor, A_norm_rank_k_cor_sc = A_norm_rank_k_cor_sc)
}
