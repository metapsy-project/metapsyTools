# blup: Empirical Bayes estimates

Best Linear Unbiased Predictions (BLUPs) for 'runMetaAnalysis' models.

## Usage

``` r
blup(x, ...)
```

## Arguments

- x:

  Model

- ...:

  Other arguments

## Details

Generates empirical Bayes (EB) estimates, also known as best linear
unbiased predictions (BLUPs), by merging the fitted values obtained from
fixed effects and estimated contributions of random effects. These
estimates represent the study-specific true effect sizes or outcomes and
are accompanied by standard errors and prediction interval bounds.
