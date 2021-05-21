# Forward pre-calculated log-risk ratios

Forwards pre-calculated log-risk ratios. Only meant to be used as part
of [`calculateEffectSizes`](calculateEffectSizes.md).

## Usage

``` r
rr.precalc(x, ...)
```

## Arguments

- x:

  data

- ...:

  Pre-calculated effect size data. Data frame must include columns
  `precalc_log_rr` and `precalc_log_rr_se`. See the [Metapsy data
  standard](https://docs.metapsy.org/data-preparation/format/).
