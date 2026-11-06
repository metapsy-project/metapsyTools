# Cluster-robust inference for CCREM models (CR1 estimator)

[`metafor::robust`](https://wviechtb.github.io/metafor/reference/robust.html)
evaluates its `cluster` argument non-standardly; this helper calls it
directly in a frame in which the cluster variable can be found.

## Usage

``` r
robustCcrem(m, study.var)
```
