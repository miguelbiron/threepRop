#' Run Generic Label Propagation (GLP) with a fixed parameter
#'
#' This function runs the iterative procedure necessary to compute GLP for a
#' given parameter \code{lambda}, exactly as in Equation 1 in Mostafavi et al.
#' (2012).
#'
#' @param M normalized affinity matrix, as returned by \code{normalize_A()}, of
#' size \code{N x N}. Must be of class \code{sparseMatrix}.
#' @param y vector of labels, of length \code{N}. Must be of class
#' \code{sparseVector}.
#' @param lambda GLP parameter (between 0 and 1).
#' @param tol tolerance for the stopping creiterion.
#' @param max_iter maximum number of iterations.
#' @param verbose print every \code{verbose} (0 means no output).
#'
#' @return A vector of length N with the estimated scores.
#'
#' @examples
#' sim_SBM = simulate_simple_SBM(N = 2500L, p_1 = 0.2, D = 0.04, R = 0.25)
#' M = normalize_A(sim_SBM$A, "asym")
#' glp(M=M, y=sim_SBM$y, lambda=0.1, max_iter = 10L, verbose = 1L)
#'
#' @references
#' Mostafavi, S., Goldenberg, A., & Morris, Q. (2012). Labeling nodes using
#' three degrees of propagation. \emph{PloS one, 7}(12), e51947.
#'
#' @export
glp = function(M, y, lambda, tol = 1e-4, max_iter = 100L, verbose = 0L){

  # we only work with sparse objects
  stopifnot(is(M, "sparseMatrix") && is(y, "sparseVector"))

  # sanity check
  stopifnot(lambda >= 0 && lambda <= 1)

  # start with random scores
  f = runif(length(y))

  # first iteration
  k = 1L
  f_new = lambda*(M %*% f) + (1-lambda)*y
  err = max(abs(f_new - f))

  # loop
  while(err > tol){
    k = k + 1L
    f_new = lambda*(M %*% f) + (1-lambda)*y
    err = max(abs(f_new - f))
    f = f_new
    if(k %% verbose == 0L) {
      cat(sprintf("Iteration %d: max_abs_change = %.4f.\n", k, err))
    }
  }

  return(as.vector(f))

}
