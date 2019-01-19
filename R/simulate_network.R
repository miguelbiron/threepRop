#' Simulate simplified Stochastic Block Model
#'
#' This function simulates a special case of the Stochastic Block Model (SBM).
#'
#' In this simplified version, there are only 2 distinct labels. Nodes with
#' different labels have an edge between them with probability \code{p_d},
#' while the probability is \code{p_s} if they have matching labels. In the
#' language of the SBM, the matrix of edge probabilities is completely specified
#' by one number associated with the diagonal (\code{p_s}), and another for the
#' off-diagonal elements (\code{p_d}).
#'
#' In spite of the above, the preferred way of specifying the model is by the
#' tuple \code{(N, p_1, D, R)}, where \code{D} is the expected density of the
#' network (between 0 and 1), and \code{R} (non-negative) controls its
#' assortativeness: for \code{R>1}, links are more common between pairs of nodes
#' with equal labels (assortative network). For \code{R<1}, the opposite is true
#' (disassortative network).
#'
#' The use can choose between supplying \code{(N, p_1, D, R)} or
#' \code{(N, p_1, p_s, p_d)}, but they cannot be mixed.
#'
#' Note that this implementation takes advantage of the structure of the SBM
#' to efficiently sample sparse random objects. One consequence of this is that
#' the elements of \code{y} are ordered, starting with all the positive nodes,
#' and then all the non-positives. The same happens with \code{A}, which is
#' partitioned into according blocks. If this is an issue, the user can always
#' permute the elements of \code{y} and \code{A} in a random but consistent
#' fashion.
#'
#' @param N number of nodes.
#' @param p_1 probability of a node being labelled as 1.
#' @param D expected density of the network, between 0 and 1.
#' @param R assortativeness parameter (non-negative). See Details.
#' @param p_d (alternative parameterization) probability of an edge existing
#' between nodes with different labels (see Details).
#' @param p_s (alternative parameterization) probability of an edge existing
#' between nodes with equal labels (see Details).
#' @param verb prints minor details about the simulation.
#'
#' @return a named list with \code{y} being a sparse vector of labels, and
#' \code{A} a sparse adjacency matrix.
#'
#' @examples
#' res = simulate_simple_SBM(N=25L, p_1=0.2, D = 0.1, R = 0.25)
#' str(res)
#' res2 = simulate_simple_SBM(N=25L, p_1=0.2, p_s = 0.001, p_d = 0.25)
#' str(res2)
#' \dontrun{
#'     # invalid combination of parameterizations
#'     simulate_simple_SBM(N=25L, p_1=0.2, D = 0.001, p_d = 0.25)
#' }
#'
#' @export
simulate_simple_SBM = function(N, p_1, D, R, p_d=NULL, p_s=NULL, verb = FALSE){

  # compute p_s and p_d if D and R are supplied
  if(is.null(p_s) && is.null(p_d)){
    p_s = (D*R)/((p_1^2 + (1-p_1)^2)*(1 + R))
    p_d = D/(2*p_1*(1-p_1)*(1 + R))
  }

  # sanity checks
  stopifnot(0 <= p_1 && p_1 <= 1 && 0 <= p_d && p_d <= 1 && 0 <= p_s && p_s <= 1)

  # sample labels
  n_pos = rbinom(n = 1L, size = N, prob = p_1) # sample number of positives
  n_neg = N - n_pos
  y = Matrix::sparseVector(i=seq_len(n_pos), length = N) # fill first n_pos elems

  # build A_pp
  if(verb) cat(sprintf("Building A_pp of size %d times %d\n", n_pos, n_pos))
  A_pp = build_A_xx(n_x = n_pos, p_edge = p_s)

  # build A_nn
  if(verb) cat(sprintf("Building A_nn of size %d times %d\n", n_neg, n_neg))
  A_nn = build_A_xx(n_x = n_neg, p_edge = p_s)

  # build A_pn, which is of size n_pos \times n_neg
  if(verb) cat(sprintf("Building A_pn of size %d times %d\n", n_pos, n_neg))
  n_edges = n_pos * n_neg
  n_on = rbinom(n = 1L, size = n_edges, prob = p_d) # number of edges on
  on = sample.int(n = n_edges, size = n_on) # sample id of edges
  on = get_coords(n_pos, n_neg, on) # coord pairs from index scanning rows
  A_pn = Matrix::sparseMatrix(
    i = on[, 1],
    j = on[, 2],
    dims = c(n_pos, n_neg)
  )

  # paste A components together and force A symmetric
  A = cbind(rbind(A_pp, Matrix::t(A_pn)), rbind(A_pn, A_nn))
  A = Matrix::forceSymmetric(A, uplo = "L")

  # convert to numerical (adds x slot equal to 1) and return in named list
  return(list(A=as(A, "dsCMatrix"), y=as(y, "dsparseVector")))

}

###############################################################################
# utils
###############################################################################

# builder for A_pp and A_nn
build_A_xx = function(n_x, p_edge){
  n_edges = choose(n_x, 2L) # number of possible edges
  n_on = rbinom(n = 1L, size = n_edges, prob = p_edge) # number of edges on
  if(n_on == 0L){
    A_xx = Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      dims = c(n_x, n_x),
      symmetric = TRUE
    )
  } else{
    on = sample.int(n = n_edges, size = n_on) # sample id of edges in lower.tri
    on = get_coords_sym(n=n_x, u=on)
    A_xx = Matrix::sparseMatrix(
      i = on[, 1],
      j = on[, 2],
      dims = c(n_x, n_x),
      symmetric = TRUE
    )
  }
  return(A_xx)
}


# get coordinates of matrix elements given indices of entries scanned row-wise
get_coords = function(n, m, u){
  cbind((u-1L) %/% m + 1L, (u-1L) %% m + 1L)
}

# get coordinates of elements in a symmetric matrix, given indices of entries
# in the lower triangle indexed by scanning columns
get_coords_sym = function(n, u){
  v = cumsum((n-1L):1L)
  res = matrix(0L, nrow = length(u), ncol = 2L)
  for(k in seq_along(u)){
    d = u[k] - v
    j = which.min(abs(d))
    if(d[j]>0L) j=j+1L
    res[k,] = c(n + d[j], j)
  }
  return(res)
}
