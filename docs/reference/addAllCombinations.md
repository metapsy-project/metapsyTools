# Add all trial arm combinations for multiarm trials in NMA.

Add all trial arm combinations for multiarm trials in NMA.

## Usage

``` r
addAllCombinations(
  data,
  vars.for.id = c("study", "outcome_type", "instrument", "time", "time_weeks", "rating"),
  vars.for.es = c("mean", "sd", "n", "mean_change", "sd_change", "n_change", "event",
    "totaln"),
  condition = "condition",
  condition.specification = "multi",
  groups.column.indicator = c("_arm1", "_arm2")
)
```
