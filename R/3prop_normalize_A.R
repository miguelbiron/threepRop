#' Normalize an affinity matrix
#'
#' This function provides different modes of normalizing an affinity matrix.
#'
#' All methods for normalization start by removing possible diagonal elements.
#' The available methods are:
#' \describe{
#'     \item{asym}{asymmetric normalization, which yields a singly stochastic
#'     matrix, where rows add to 1.}
#'     \item{sym}{symmetric normalization, yielding a symmetric matrix.}
#'     \item{rm_diag}{only removes possible diagonal elements.}
#' }
#'
#' @param A a symmetric affinity matrix. Can be dense or sparse from package
#' \code{Matrix}.
#' @param M_type a character string with one of "asym", "sym", or "as_is" (see
#' details).
#'
#' @return A sparse normalized affinity matrix.
#'
#' @examples
#' # simulate_simple_SBM does not produce self connections
#' A = as.matrix(simulate_simple_SBM(25L, 0.2, 0.001, 0.25)$A)
#' all.equal(normalize_A(A, "rm_diag"), A)
#'
#' @references
#' Mostafavi, S., Goldenberg, A., & Morris, Q. (2012). Labeling nodes using
#' three degrees of propagation. \emph{PloS one, 7}(12), e51947.
#'
#' @export
normalize_A = function(A, M_type){

  if(!Matrix::isSymmetric(A)){stop("Matrix A is not symmetric.")}

  # copy A and remove diagonal elements if present
  M = A
  if(max(abs(Matrix::diag(M))) > 0){ Matrix::diag(M) = 0 }

  # normalize using chosen method
  switch(
    EXPR    = M_type,
    asym    = get_asym_M(M),
    sym     = get_sym_M(M),
    rm_diag = M,
    stop(sprintf("'%s' is an unrecognized type of normalization", M_type))
  )

}

#######################################
# symmetric variant (S)
#######################################

get_sym_M = function(A){
  d = Matrix::rowSums(A)
  d = ifelse(d == 0, 1, d) # do not alter rows/columns (A is sym) which sum to 0
  D_inv_sq = Matrix::Diagonal(x = 1/sqrt(d))
  return(D_inv_sq %*% A %*% D_inv_sq)
}

#######################################
# asymmetric variant (P)
#######################################

get_asym_M = function(A){
  d = Matrix::rowSums(A)
  d = ifelse(d == 0, 1, d) # do not alter rows which sum to 0
  return(Matrix::Diagonal(x = 1/d) %*% A)
}
