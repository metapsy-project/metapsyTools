# Calculate Hedges' g using binary outcome data

Calculates Hedges' g from binary outcome data. Only meant to be used as
part of [`calculateEffectSizes`](calculateEffectSizes.md).

## Usage

``` r
g.binary(x, cc = TRUE, ...)
```

## Arguments

- x:

  data

- cc:

  Should a continuity correction for zero cells be applied? Either
  `TRUE` (default; applies the Sweeting et al. (2004) treatment-arm
  continuity correction, TACC), `FALSE` (no correction), or a single
  numeric increment to be added to all four cells (e.g., `0.5`).

- ...:

  Binary effect size data. Data frame must include columns `event_arm1`,
  `event_arm2`, `totaln_arm1`, `totaln_arm2`. See the [Metapsy data
  standard](https://docs.metapsy.org/data-preparation/format/).
