#' Estimate 3prop weights using double cross-validation
#'
#' This function estimates the weights for the 3prop algorithm using a double
#' cross-validation procedure.
#'
#' By double cross-validation (CV), we mean that there are two CV loops, one
#' nested within the other. The first loop deals with calculating the random
#' walk features, while the second loop deals with estimating the coefficients
#' alpha associated to those features.
#'
#' Additionally, the function returns computes the area under the ROC curve
#' (AUROC) for every fold in the outer CV loop.
#'
#' @param M normalized affinity matrix, as returned by \code{normalize_A()}, of
#' size \code{N x N}. Must be of class \code{sparseMatrix}.
#' @param y vector of labels, of length \code{N}. Must be of class
#' \code{sparseVector}.
#' @param R maximum length of the random walks to consider. The default 3 is the
#' usual for 3prop.
#' @param n_folds number of CV folds to use in both inner and outer CV loops.
#' @param reg regularization parameter for LDA.
#' @param method string for method to compute the coefficients. Can be "LDA"
#' (default) or "Ridge". Both use the parameter \code{reg}.
#'
#' @return A list with two elements: a matrix of size \code{R x n_folds},
#' containing the weights estimated for each (outer) CV iteration, and a vector
#' of length \code{n_folds} containing the AUROC for each such iteration.
#'
#' @examples
#' sim_SBM = simulate_simple_SBM(N = 2500L, p_1 = 0.2, D = 0.04, R = 0.25)
#' M = normalize_A(sim_SBM$A, "asym")
#' three_prop_cv(M=M, y=sim_SBM$y)
#' three_prop_cv(M=M, y=sim_SBM$y, method = "Ridge")
#'
#' @references
#' Mostafavi, S., Goldenberg, A., & Morris, Q. (2012). Labeling nodes using
#' three degrees of propagation. \emph{PloS one, 7}(12), e51947.
#'
#' @export
three_prop_cv = function(M, y, R = 3L, n_folds = 3L, reg = 1e-09, method = "LDA"){

  # we only work with sparse objects
  stopifnot(is(M, "sparseMatrix") && is(y, "sparseVector"))

  # useful quantities
  n = length(y)
  ind_all = seq_len(n)

  # sample a random partition of 1:n into n_folds subsets
  fold_id = sample.int(n_folds, size = n, replace = TRUE)

  # define storage for things to output
  alpha_mat = matrix(nrow = R, ncol = n_folds)
  AUROC_vec = numeric(n_folds)

  # loop folds
  for(j in seq_len(n_folds)){
    ind = which(fold_id == j) # which observations are in the current fold

    # copy y but set to zero labels for current fold
    y_cv = y; y_cv[ind] = 0

    # define relevant indices
    pos = y_cv@i # positives not in current test set
    neg = setdiff(ind_all, c(pos, ind)) # non-positives not in current test set

    # estimate coefficients
    alpha = est_weights_cv(
      M        = M,
      R        = R,
      y        = y_cv,
      pos      = pos,
      neg      = neg,
      n_folds  = n_folds,
      reg      = reg,
      method   = method
    )
    alpha_mat[, j] = alpha / sum(abs(alpha)) # normalize and store

    # compute X
    X = matrix(0, nrow = n, ncol = R)
    X[,1L] = as.vector(M %*% y_cv)
    if(R >= 2L){
      for(r in 2L:R){ X[,r] = as.vector(M %*% X[,r-1L]) }
    }

    # compute f_score and convert to vector
    f_scores = as.vector(X %*% alpha_mat[, j])

    # compute AUROC
    AUROC_vec[j] = AUROC(p=f_scores[ind], y=y[ind], plot_roc=FALSE, warn_inv=FALSE)

  }

  # return named list
  return(list(alpha=alpha_mat, AUROC=AUROC_vec))

}
