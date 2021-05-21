# Profile Likelihood Plots for 'runMetaAnalysis' models.

Profiles the restricted log-likelihood of `threelevel`,
`threelevel.che`, `ccrem`, or `ccrem.che` models, using the
[`metafor::profile.rma()`](https://wviechtb.github.io/metafor/reference/profile.rma.html)
function. This functionality can be used to check if the heterogeneity
variances (\\\tau^2\\ within and between studies) were identifiable and
correctly estimated.

## Usage

``` r
# S3 method for class 'runMetaAnalysis'
profile(fitted, which = NULL, ...)
```

## Arguments

- fitted:

  An object of class `runMetaAnalysis`.

- which:

  Model for which estimates should be printed. Can be one of
  `"threelevel"` `"threelevel.che"`, `"ccrem"`, or `"ccrem.che"`.

- ...:

  Additional arguments.
