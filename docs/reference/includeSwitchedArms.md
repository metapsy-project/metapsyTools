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
  includes calculated effect sizes and standard errors in columns `es`
  and `se`, respectively. Only data sets created by
  `calculateEffectSizes` with `trt.indicator` set to `"trt"` can be
  used.

- ...:

  Further arguments (not used).
