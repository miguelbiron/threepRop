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
#' size \code{N x N}. Can be dense or sparse from package \code{Matrix}.
#' @param y vector of labels, of length \code{N}. Can be dense or sparse from
#' package \code{Matrix}.
#' @param R maximum length of the random walks to consider. The default 3 is the
#' usual for 3prop.
#' @param n_folds number of CV folds to use in both inner and outer CV loops.
#'
#' @return A list with two elements: a matrix of size \code{R x n_folds},
#' containing the weights estimated for each (outer) CV iteration, and a vector
#' of length \code{n_folds} containing the AUROC for each such iteration.
#'
#' @examples
#' sim_SBM = simulate_simple_SBM(25L, 0.2, 0.001, 0.25)
#' M = normalize_A(sim_SBM$A, "asym")
#' three_prop_cv(M=M, y=sim_SBM$y)
#'
#' @references
#' Mostafavi, S., Goldenberg, A., & Morris, Q. (2012). Labeling nodes using
#' three degrees of propagation. \emph{PloS one, 7}(12), e51947.
#'
#' @export
three_prop_cv = function(M, y, R = 3L, n_folds = 3L){

  # useful quantities
  n = length(y)
  ind_all = seq_len(n)

  # compute list of M powers
  M_powers = vector(mode = "list", length = R)
  M_powers[[1L]] = M
  for(r in 2L:R){ M_powers[[r]] = M %*% M_powers[[r-1L]]}

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
    pos = Matrix::which(y_cv > 0) # positives not in current test set
    neg = setdiff(ind_all, c(pos, ind)) # non-positives not in current test set

    # estimate coefficients
    alpha = est_weights_lda_cv(
      M_powers = M_powers,
      y        = y_cv,
      pos      = pos,
      neg      = neg,
      n_folds  = n_folds
    )
    alpha_mat[, j] = alpha / sum(abs(alpha)) # normalize and store

    # # compute scores using loop to avoid creating X explicitly
    # # using a loop shouldn't be a problem, since R is very small (3 for 3prop)
    # f_scores = 0*y # create empty vector of the same type as y (dense or sparse)
    # for(r in seq_len(R)){
    #   f_scores = f_scores + alpha_mat[r, j]*(M_powers[[r]] %*% y_cv)
    # }

    # compute X
    X = do.call(cbind,
                lapply(M_powers, FUN = function(B){B %*% y_cv})
                )

    # compute f_score and convert to vector
    f_scores = as.vector(X %*% alpha_mat[, j])

    # compute AUROC
    AUROC_vec[j] = AUROC(p=f_scores[ind], y=y[ind], plot_roc=FALSE, warn_inv=FALSE)

  }

  # return named list
  return(list(alpha=alpha_mat, AUROC=AUROC_vec))

}
