# Evaluate a risk of bias filtering statement within a data set

Evaluates statements such as `"rob > 2"`, `"rob>2"`,
`"rob %in% c(3,4)"`, or `"rob > 2 & ac == 1"` within `data`. Rows for
which the statement cannot be evaluated because a variable is missing
are not selected.

## Usage

``` r
evalRobFilter(filter, data, purpose = "data")
```
