# `threepRop`: R Implementation of the 3prop Label Propagation Algorithm

## Description

R implementation of the algorithm proposed by Mostafavi et al. (2012). This is a port of the MATLAB code provided by the authors, available [here](http://ai.stanford.edu/~saram/3propcode.zip). For testing purposes, the example dataset contained in the link has been incorporated into `threepRop` as a dataset called `example_dta`.

Additionally, the package provides a function to efficiently simulate Stochastic Block Models for two classes, whose output can be used to test 3prop (view examples in the help page of the main function `three_prop_cv`). Also, `threepRop` contains a utiliy for calculating the area under the ROC curve (AUROC) for performance assessment.

For completeness, we also provide a simple implementation of the Generic Label Propagation (GLP) algorithm in function `glp`, as described by Equation 1 in Mostafavi et al. (2012).

The implementation relies on sparse data classes of package `Matrix` for the linear algebra operations needed on sparse matrices and vectors. Both `three_prop_cv` and `glp` only accept inputs from these classes, and will throw errors when used with dense objects.

## Installation

``` r
if (!require(devtools)) {
    install.packages('devtools')
}
devtools::install_github('miguelbiron/threepRop')
```

## References

Mostafavi, S., Goldenberg, A., & Morris, Q. (2012). Labeling nodes using three degrees of propagation. _PloS one, 7_(12), e51947.
