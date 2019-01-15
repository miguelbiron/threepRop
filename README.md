# `threepRop`: R Implementation of the 3prop Label Propagation Algorithm

## Description

R implementation of the algorithm proposed by Mostafavi et al. (2012). This is a port of the MATLAB code provided by the authors, available [here](http://ai.stanford.edu/~saram/3propcode.zip).

## Installation

``` r
if (!require(devtools)) {
    install.packages('devtools')
}
devtools::install_github('miguelbiron/threepRop')
```

## References

Mostafavi, S., Goldenberg, A., & Morris, Q. (2012). Labeling nodes using three degrees of propagation. _PloS one, 7_(12), e51947.
