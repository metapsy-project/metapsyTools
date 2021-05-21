# Calculate the log-risk ratio using binary outcome data

Calculate the log risk ratio using binary outcome data. Only meant to be
used as part of [`calculateEffectSizes`](calculateEffectSizes.md).

## Usage

``` r
rr.binary(x, cc = 0.5, ...)
```

## Arguments

- x:

  data

- cc:

  Should a continuity correction for zero cells be applied? Either
  `FALSE` or the increment to be added. Default is 0.5.

- ...:

  Binary effect size data. Data frame must include columns `event_arm1`,
  `event_arm2`, `totaln_arm1`, `totaln_arm2`. See the [Metapsy data
  standard](https://docs.metapsy.org/data-preparation/format/).
