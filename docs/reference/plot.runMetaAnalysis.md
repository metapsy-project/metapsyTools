# Plot method for objects of class 'runMetaAnalysis'

Plot S3 method for objects of class `runMetaAnalysis`.

## Usage

``` r
# S3 method for class 'runMetaAnalysis'
plot(x, which = NULL, eb = FALSE, eb.labels = FALSE, ...)
```

## Arguments

- x:

  An object of class `runMetaAnalysis`.

- which:

  Model to be plotted. Can be one of `"overall"`, `"combined"`,
  `"lowest.highest"`, `"outliers"`, `"influence"`, `"threelevel"`,
  `"threelevel.che"`, `"ccrem"`, `"ccrem.che"`, `"baujat"`, `"loo-es"`,
  `"loo-i2"`, `"trimfill"`, `"limitmeta"` or `"selection"`.

- eb:

  Prints a forest plot with empirical Bayes point estimates and
  study-specific prediction intervals as proposed by van Aert (2021).
  Defaults to `FALSE`.

- eb.labels:

  If `eb` is `TRUE`, should the empirical Bayes estimates and prediction
  intervals for each study be printed on the right side of the forest
  plot? Defaults to `FALSE`.

- ...:

  Additional arguments.

## References

van Aert, R. C., Schmid, C. H., Svensson, D., & Jackson, D. (2021).
Study specific prediction intervals for random-effects meta-analysis: A
tutorial. *Research Synthesis Methods, 12*(4), 429-447.

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Paula Kuper
<paula.r.kuper@gmail.com>, Pim Cuijpers <p.cuijpers@vu.nl>
