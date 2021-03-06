% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glp.R
\name{glp}
\alias{glp}
\title{Run Generic Label Propagation (GLP) with a fixed parameter}
\usage{
glp(M, y, lambda, tol = 1e-04, max_iter = 100L, verbose = 0L)
}
\arguments{
\item{M}{normalized affinity matrix, as returned by \code{normalize_A()}, of
size \code{N x N}. Must be of class \code{sparseMatrix}.}

\item{y}{vector of labels, of length \code{N}. Must be of class
\code{sparseVector}.}

\item{lambda}{GLP parameter (between 0 and 1).}

\item{tol}{tolerance for the stopping creiterion.}

\item{max_iter}{maximum number of iterations.}

\item{verbose}{print every \code{verbose} (0 means no output).}
}
\value{
A vector of length N with the estimated scores.
}
\description{
This function runs the iterative procedure necessary to compute GLP for a
given parameter \code{lambda}, exactly as in Equation 1 in Mostafavi et al.
(2012).
}
\examples{
sim_SBM = simulate_simple_SBM(N = 2500L, p_1 = 0.2, D = 0.04, R = 0.25)
M = normalize_A(sim_SBM$A, "asym")
glp(M=M, y=sim_SBM$y, lambda=0.1, max_iter = 10L, verbose = 1L)

}
\references{
Mostafavi, S., Goldenberg, A., & Morris, Q. (2012). Labeling nodes using
three degrees of propagation. \emph{PloS one, 7}(12), e51947.
}
