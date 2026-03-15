# Calculate ratio of means (log scale) from change-score means and SDs

Only meant to be used as part of
[`calculateEffectSizes`](calculateEffectSizes.md).

## Usage

``` r
rom.change.m.sd(x, ...)
```

## Arguments

- x:

  data

- ...:

  columns `mean_change_arm1`, `mean_change_arm2`, `sd_change_arm1`,
  `sd_change_arm2`, `n_change_arm1`, `n_change_arm2`.

## Value

data.frame with `es` (log ROM) and `se`.
