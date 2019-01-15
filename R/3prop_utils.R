#' Compute area under the ROC curve
#'
#' Function for computing the area under the ROC curve (AUROC)
#'
#' It is assumed that a higher score increases the probability of an object
#' having label 1.
#'
#' @param p vector of scores. Can be sparse (from package Matrix), but will be
#' coerced to dense.
#' @param y vector of labels (0 or 1). Can be sparse (from package Matrix), but
#' will be coerced to dense.
#' @param plot_roc if \code{TRUE}, plots the ROC curve.
#' @param warn_inv raises a warning if high values of \code{p} are associated
#' with the 0 label, which produces an AUROC lower than 0.5.
#'
#' @return The value of the AUROC.
#'
#' @export
AUROC = function(p, y, plot_roc = TRUE, warn_inv = TRUE){

  stopifnot(length(p) == length(y)) # sanity check
  N = length(y)

  # coerce to dense vectors
  if(is(y, "sparseVector")) y=as.vector(y)
  if(is(p, "sparseVector")) p=as.vector(p)


  # build and sort data.frame by decreasing p
  ind = order(p, decreasing = TRUE)
  dta = data.frame(y=y[ind],p=p[ind])

  # compute cumulative sums
  dta$sum_one = cumsum(dta$y > 0)
  dta$sum_zero = cumsum(dta$y < 1)

  # compute TPR and FPR
  dta$TPR = dta$sum_one / dta$sum_one[N]
  dta$FPR = dta$sum_zero / dta$sum_zero[N]

  # compute AUC
  dta$dFPR = c(dta$FPR[1L], diff(dta$FPR))
  auc = sum(dta$dFPR*dta$TPR)

  if(warn_inv && auc<0.49){ warning("AUC<0.5 because p means P(y==0) (inverted).")}

  # plot
  if(plot_roc){
    op=par(mar = c(5.1, 4.1, 2.1, 2.1))
    plot(dta$FPR, dta$TPR, type = "l", asp = 1,
         xlab = "False Positive Rate (1-Specificity)",
         ylab = "True Positive Rate (Sensitivity)")
    abline(a=0,b=1, lty="dotted")
    par(op)
  }

  return(auc)

}

# # test
# n = 100000L
# b = 1.8 # higher b means higher AUC
# x = sin(1:n)
# p = 1/(1+exp(-b*x))
# y = 1*(runif(n) < p)
# as.numeric(pROC:::auc.default(y, p)) == AUROC(p, y)
# microbenchmark::microbenchmark(
#   pROC = pROC:::auc.default(y, p),
#   mine = AUROC(p, y)
# ) # mine is about 22 times faster for n = 100000!



##############################################################################
# expand cov function (covariance) to accomodate sparse matrices
# returns dense matrix because X is assumed to be n \times p, p small
##############################################################################

cov_sparse = function(X){
  # if not sparse, use base method
  if(!is(X, 'sparseMatrix')) return(cov(X))

  n = nrow(X)

  # compute column means as dense vector
  mean_v = Matrix::colMeans(X, sparseResult=FALSE)

  # return dense covariance matrix
  return((1/n)*as.matrix(Matrix::crossprod(X)) - tcrossprod(mean_v))
}

# # test
# S = matrix(c(1, 0.5, 0.5, 1), nrow = 2L)
# X = MASS::mvrnorm(n=1000L, mu = c(12,-63), Sigma = S)
# X_s = as(X, "sparseMatrix")
# all.equal(cov_sparse(X), cov(X))
