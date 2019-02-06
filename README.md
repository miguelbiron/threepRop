# `threepRop`: R Implementation of the 3prop Label Propagation Algorithm

## Description

R implementation of the algorithm proposed by Mostafavi et al. (2012). This is a port of the MATLAB code provided by the authors, available [here](http://ai.stanford.edu/~saram/3propcode.zip).

Additionally, the package provides a function to efficiently simulate Stochastic Block Models for two classes, whose output can be used to test 3prop (view examples in the help page of the main function `three_prop_cv`). Additionally, it contains a utiliy for calculating the area under the ROC curve (AUROC) for performance assessment.

We also provide a simple implementation of the Generic Label Propagation (GLP) algorithm in function `glp`, as described by Equation 1 in Mostafavi et al. (2012).

## Installation

``` r
if (!require(devtools)) {
    install.packages('devtools')
}
devtools::install_github('miguelbiron/threepRop')
```

## References

Mostafavi, S., Goldenberg, A., & Morris, Q. (2012). Labeling nodes using three degrees of propagation. _PloS one, 7_(12), e51947.
