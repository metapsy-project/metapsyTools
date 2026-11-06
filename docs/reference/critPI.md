# Critical value for prediction intervals of multilevel models

Models fitted with `test = "z"` have no degrees of freedom (`ddf` is
`NA`), so the standard normal quantile has to be used in that case.

## Usage

``` r
critPI(m, level = 0.95)
```
