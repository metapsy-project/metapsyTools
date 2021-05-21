# Forward pre-calculated values of Hedges' g

Forwards pre-calculated values of Hedges' g. Only meant to be used as
part of [`calculateEffectSizes`](calculateEffectSizes.md).

## Usage

``` r
g.precalc(x, ...)
```

## Arguments

- x:

  data

- ...:

  Pre-calculated effect size data. Data frame must include columns
  `precalc_g` and `precalc_g_se`. See the [Metapsy data
  standard](https://docs.metapsy.org/data-preparation/format/).
