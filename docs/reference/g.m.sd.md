# Calculate Hedges' g using means and standard deviations

Calculate Hedges' g using Mean and Standard Deviation. Only meant to be
used as part of [`calculateEffectSizes`](calculateEffectSizes.md).

## Usage

``` r
g.m.sd(x, ...)
```

## Arguments

- x:

  data

- ...:

  Effect size data. Data frame must include columns `mean_arm1`,
  `mean_arm2`, `sd_arm1`, `sd_arm2`, `n_arm1`, `n_arm2`. See the
  [Metapsy data
  standard](https://docs.metapsy.org/data-preparation/format/).
