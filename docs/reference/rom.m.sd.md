# Calculate ratio of means (log scale) using means and standard deviations

Only meant to be used as part of
[`calculateEffectSizes`](calculateEffectSizes.md). Returns
log(mean_arm1/mean_arm2) and its SE (delta method).

## Usage

``` r
rom.m.sd(x, ...)
```

## Arguments

- x:

  data

- ...:

  columns `mean_arm1`, `mean_arm2`, `sd_arm1`, `sd_arm2`, `n_arm1`,
  `n_arm2`.

## Value

data.frame with `es` (log ROM) and `se`.
