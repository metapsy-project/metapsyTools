# Sampling variance-covariance matrix of the CHE working model

Constant sampling correlation `r` within clusters, and the average
sampling variance of a cluster for all of its effect sizes (Pustejovsky
& Tipton, 2022). Replaces
`clubSandwich::impute_covariance_matrix(..., smooth_vi = TRUE)`, which
is deprecated. Rows must be sorted by cluster.

## Usage

``` r
imputeCheV(vi, cluster, r)
```
