# Calculate the log-risk ratio using binary outcome data

Calculate the log risk ratio using binary outcome data. Only meant to be
used as part of [`calculateEffectSizes`](calculateEffectSizes.md).

## Usage

``` r
rr.binary(x, cc = TRUE, ...)
```

## Arguments

- x:

  data

- cc:

  Should a continuity correction for zero cells be applied? Either
  `TRUE` (default, applies TACC) or `FALSE`.

- ...:

  Binary effect size data. Data frame must include columns `event_arm1`,
  `event_arm2`, `totaln_arm1`, `totaln_arm2`. See the [Metapsy data
  standard](https://docs.metapsy.org/data-preparation/format/).

## Details

When `cc = TRUE`, a treatment-arm continuity correction (TACC) is
applied when any cell in the 2x2 table has zero count. Following
[Sweeting et al. (2004)](https://doi.org/10.1002/sim.1761), a correction
of `k = n_opposite / (n_arm1 + n_arm2)` is added to both cells of each
arm (i.e., `totaln_arm2 / (totaln_arm1 + totaln_arm2)` to the arm1 cells
and `totaln_arm1 / (totaln_arm1 + totaln_arm2)` to the arm2 cells)
before calculating `es` and `se`. The correction is not applied to the
stored cell counts. When arms are balanced this reduces to adding 0.5 to
all cells, but unlike a fixed 0.5 correction it remains unbiased under
unequal arm sizes.
