# Best Linear Unbiased Predictions (BLUPs) for 'runMetaAnalysis' models.

Generates empirical Bayes (EB) estimates, also known as best linear
unbiased predictions (BLUPs), by merging the fitted values obtained from
fixed effects and estimated contributions of random effects. These
estimates represent the study-specific true effect sizes or outcomes and
are accompanied by standard errors and prediction interval bounds. Uses
the
[`metafor::blup.rma.uni()`](https://wviechtb.github.io/metafor/reference/blup.html)
function internally.

## Usage

``` r
# S3 method for class 'runMetaAnalysis'
eb(x, which = NULL, ...)
```

## Arguments

- x:

  An object of class `runMetaAnalysis`.

- which:

  Model for which estimates should be printed. Can be one of
  `"overall"`, `"combined"`, `"lowest.highest"`, `"outliers"`,
  `"influence"`, `"threelevel"`, `"threelevel.che"`, `"ccrem"`, or
  `"ccrem.che"`.

- ...:

  Additional arguments.
