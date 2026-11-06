# Include information of rows with switched reference group

Adds effect size data and study information of rows with switched
reference arms.

## Usage

``` r
includeSwitchedArms(dat, ...)
```

## Arguments

- dat:

  Data set created by `calculateEffectSizes` in the wider format, which
  includes calculated effect sizes and standard errors in columns such
  as `.g` and `.g_se`. Only data sets created by `calculateEffectSizes`
  with `.trt.indicator` set to `"arm"` can be used. Signs of `.g`,
  `.log_rr` and (if available) `.log_rom` are reversed in the added
  rows.

- ...:

  Further arguments (not used).
