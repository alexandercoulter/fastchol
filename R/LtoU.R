#' Fast L to U Transpose
#'
#' @description
#' This function copies the lower-triangular part of square matrix `L`
#' (including the main diagonal) into the upper-triangular portion of output
#' matrix `U`, other values zero. For `L` that is only lower-triangular, this is
#' effectively a faster transpose operation.
#' 
#' @param L Any square matrix such that its lower-triangular part (including
#' main diagonal) will be transposed to create an upper-triangular matrix.
#'
#' @returns An upper-triangular square matrix whose non-zero entries are the
#' lower-triangular part of `L`.
#' @export
#' @import Rcpp, RcppArmadillo
#'
#' @examples
LtoU = function(L){
  
  U = matrix(0, nrow(L), ncol(L))
  LtoU_Rcpp(U, L)
  U
  
}