# Calculate Hedges' g using within-group change data

Calculate Hedges' g based on change data. Only meant to be used as part
of [`calculateEffectSizes`](calculateEffectSizes.md).

## Usage

``` r
g.change.m.sd(x, ...)
```

## Arguments

- x:

  data

- ...:

  Change score effect size data. Data frame must include columns
  `mean_change_arm1`, `mean_change_arm2`, `sd_change_arm1`,
  `sd_change_arm2`, `n_change_arm1`, `n_change_arm2`. See the [Metapsy
  data standard](https://docs.metapsy.org/data-preparation/format/).
