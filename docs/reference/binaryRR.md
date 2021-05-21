# Calculate the log-risk ratio using binary outcome data.

Calculate the log risk ratio using binary outcome data. Only meant to be
used as part of [`calculateEffectSizes`](calculateEffectSizes.md).

## Usage

``` r
binaryRR(x, ...)
```

## Arguments

- x:

  data

- ...:

  Binary effect size data. Must be `Improved_N_trt1`, `Improved_N_trt2`,
  `Rand_N_trt1`, `Rand_N_trt2` (all `numeric`).
