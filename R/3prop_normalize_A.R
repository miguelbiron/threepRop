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
#'     \item{rm_diag}{simply returns \code{A} without doing anything else than
#'     removing possible diagonal elements.}
#' }
#'
#' @param A a symmetric affinity matrix. Can be dense or sparse from package
#' \code{Matrix}.
#' @param M_type a character string with one of "asym", "sym", or "as_is" (see
#' details).
#'
#' @return A normalized affinity matrix.
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

  # remove diagonal elements
  M = A - Matrix::diag(Matrix::diag(A))

  # normalize using chosen method
  M = switch(
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
  d_inv_sqrt = 1/sqrt(d)

  # multiply each row
  M = sweep(
    x      = A,
    MARGIN = 1L,
    STATS  = d_inv_sqrt,
    FUN    = "*",
    check.margin = FALSE
  )

  # multiply each column
  M = sweep(
    x      = M,
    MARGIN = 2L,
    STATS  = d_inv_sqrt,
    FUN    = "*",
    check.margin = FALSE
  )

  return(M)
}

# # test
# n = 500L
# A = matrix(runif(n*n), nrow = n)
# D_inv_sqrt = diag(1/sqrt(Matrix::rowSums(A)))
# all.equal(get_sym_M(A), D_inv_sqrt %*% A %*% D_inv_sqrt)
# microbenchmark::microbenchmark(
#   sweep  = {get_sym_M(A)},
#   matrix = {
#     D_inv_sqrt = diag(1/sqrt(Matrix::rowSums(A)))
#     D_inv_sqrt %*% A %*% D_inv_sqrt
#   }
# ) # sweep version about 40 times faster for n=500L (for smaller n they are =)

#######################################
# asymmetric variant (P)
#######################################

get_asym_M = function(A){
  d = Matrix::rowSums(A)
  d = ifelse(d == 0, 1, d) # do not alter rows which sum to 0
  # normalize rows to sum to 1
  M = sweep(
    x      = A,
    MARGIN = 1L,
    STATS  = 1 / d,
    FUN    = "*",
    check.margin = FALSE
  )
  return(M)
}

# # test
# n = 500L
# A = matrix(runif(n*n), nrow = n)
# all.equal(c(1,1), range(Matrix::rowSums(get_asym_M(A)))) # check it's a Markov transition matrix
# all.equal(get_asym_M(A), diag(1/Matrix::rowSums(A)) %*% A)
# microbenchmark::microbenchmark(
#   sweep  = {get_asym_M(A)},
#   matrix = {diag(1/Matrix::rowSums(A)) %*% A}
# ) # sweep version about 25 times faster for n=500L (for smaller n they are =)
