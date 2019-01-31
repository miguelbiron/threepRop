est_weights_lda_cv = function(M, R, y, pos, neg, n_folds, reg){
  # M: sparseMatrix of size N times N
  # y : labels in {0,1}, with components of outer CV zeroed
  # pos: indices of positive instances in y which are not in the outer CV fold
  # neg: indices of non-positive instances which are not in the outer CV fold
  # n_folds: number of CV folds
  # reg: small positive number to ensure positive definitiness of cov matrix
  # OUTPUT: vector of alpha coefficients (length R)

  n = length(y)

  # sample a random partition of 1:n into n_folds subsets
  fold_id = sample.int(n_folds, size = n, replace = TRUE)

  # define matrix to store alpha vectors
  alpha_mat = matrix(nrow = R, ncol = n_folds)

  # loop folds
  for(k in seq_len(n_folds)){

    ind = which(fold_id == k) # which observations are in the current fold
    ind_test_neg = intersect(neg,ind) # is neg, not in outer fold, and in this fold
    ind_test_pos = intersect(pos,ind) # is pos, not in outer fold, and in this fold

    # set to zero labels for current fold
    # note that y already comes with components of outer CV zeroed
    y_cv = y; y_cv[ind] = 0

    # compute X
    X = matrix(0, nrow = n, ncol = R)
    X[,1L] = as.vector(M %*% y_cv)
    for(r in 2L:R){ X[,r] = as.vector(M %*% X[,r-1L]) }

    # compute cov
    # uses only nodes in this fold but not in the outer fold
    C = cov(X[c(ind_test_pos, ind_test_neg),])

    # compute averages for covariates in pos and neg sets
    # uses only nodes in this fold but not in the outer fold
    if(length(ind_test_pos)==1L){
      x_p = X[ind_test_pos,]
    } else{
      x_p = colMeans(X[ind_test_pos,])
    }
    if(length(ind_test_neg)==1L){
      x_n = X[ind_test_neg,]
    } else{
      x_n = colMeans(X[ind_test_neg,])
    }

    # compute alpha with regularized covariance matrix
    alpha_mat[,k] = solve(C+reg*diag(R)) %*% (x_p - x_n)

  }
  return(rowMeans(alpha_mat, na.rm = TRUE))
}

