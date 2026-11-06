# Pooled effect, confidence and prediction interval of multilevel models

Implements the same inference as shown in the summary table: model-based
if `use.rve` is `FALSE`, the CR2 estimator for three-level (CHE) models,
and the CR1 estimator, clustered on studies, for CCREMs (for which
clubSandwich cannot detect the clustering structure).

## Usage

``` r
multilevelInference(m, use.rve = TRUE)
```
