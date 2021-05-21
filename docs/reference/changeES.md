# Calculate Hedges' g based on change data.

Calculate Hedges' g based on change data. Only meant to be used as part
of [`calculateEffectSizes`](calculateEffectSizes.md).

## Usage

``` r
changeES(x, ...)
```

## Arguments

- x:

  data

- ...:

  Change score effect size data. Must be `Change_m_trt1`,
  `Change_m_trt2`, `Change_SD_trt1`, `Change_SD_trt2`, `Change_N_trt1`,
  `Change_N_trt2` (all `numeric`).
