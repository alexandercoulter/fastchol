#' Fast Triangular Transpose
#'
#' @description
#' This function copies the lower- or upper-triangular part of square matrix `M`
#' (including the main diagonal) into the upper-triangular (or respectively,
#' lower-triangular part) portion of output
#' matrix `U`, other values zero. For `L` that is only lower-triangular, this is
#' effectively a faster transpose operation.
#' 
#' @param M Any square matrix such that its lower- (upper-)triangular part,
#' including main diagonal, will be transposed to create an upper- 
#' (lower-)triangular matrix.
#' @param lower_to_upper Boolean, dictating whether the lower-triangular part of
#' `M` will be transposed to the upper-triangular part (`TRUE`, default), or
#' whether the upper-triangular part of `M` will be transposed to the
#' lower-triangular part (`FALSE`)
#' 
#' @returns The transpose
#' @export
#'
#' @examples
#' # Generate positive definite matrix
#' p = 20
#' M = diag(p) + tcrossprod(matrix(rnorm(p * 2), p, 2))
#' 
#' # Calculate lower and upper Cholesky factors of M
#' U = chol(M)
#' L = t(U)
#' 
#' # Transpose and compare
#' max(abs(L - trit(U)))
#' max(abs(U - trit(L)))
trit = function(M, lower_to_upper = TRUE){
  
  R = matrix(0, nrow(M), ncol(M))
  if(lower_to_upper) LtoU_Rcpp(R, M) else UtoL_Rcpp(R, M)
  R
  
}