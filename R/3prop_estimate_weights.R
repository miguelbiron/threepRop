est_weights_lda_cv = function(M_powers, y, pos, neg, n_folds = 3L){
  # M_powers: list of powers of M (length R)
  # y : labels in {0,1}, with components of outer CV zeroed
  # pos: indices of positive instances in y which are not in the outer CV fold
  # neg: indices of non-positive instances which are not in the outer CV fold
  # n_folds: number of CV folds
  
  n = length(y); R = length(M_powers)
  
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
    X = do.call(cbind,
                lapply(M_powers, function(B){B %*% y_cv})
    )
    
    # compute cov, allowing for sparse X
    # note: does not consider all nodes
    C = cov_sparse(X[c(ind_test_pos, ind_test_neg),])
    
    # compute alpha
    alpha_mat[,k]=solve(C)%*%(colMeans(X[ind_test_pos,])-colMeans(X[ind_test_neg,]))
    
  }
  return(rowMeans(alpha_mat))
}

