# Test for baseline imbalances

This function returns a SMD for baseline imbalance tests, along with 99%
confidence intervals. If the `cluster` argument is specified, *p*-values
can be adjusted using methods available in
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

## Usage

``` r
testBaselineImbalance(
  m_arm1,
  m_arm2,
  sd_arm1,
  sd_arm2,
  n_arm1,
  n_arm2,
  cluster,
  study,
  data = NULL,
  p.adj = "holm",
  ...
)
```

## Arguments

- m_arm1:

  A vector with the baseline mean(s) in the first arm.

- m_arm2:

  A vector with the baseline mean(s) in the second arm.

- sd_arm1:

  A vector with the baseline standard deviation(s) in the first arm.

- sd_arm2:

  A vector with the baseline standard deviation(s) in the second arm.

- n_arm1:

  A vector with the baseline sample size(s) in the first arm.

- n_arm2:

  A vector with the baseline sample size(s) in the second arm.

- cluster:

  An optional vector of the same length as `m_arm1`, `m_arm2`, etc.,
  encoding to which larger cluster/study a comparison belongs to. If
  specified, *p* values of the comparisons will be adjusted according to
  the method provided in the `p.adj` argument.

- study:

  An optional vector with study labels.

- data:

  An optional `data.frame` that includes the effect size data.

- p.adj:

  A character string, specifying the type of *p*-value adjustment. For
  available options, see the "Details" section in
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Defaults
  to `"holm"`.

- ...:

  Additional arguments.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example data that follows the Metapsy data standard
data("depressionPsyCtr")

# This is an unexported function, so we need the metapsyTools::: prefix
# Test for differences without p-value adjustment
metapsyTools:::testBaselineImbalance(
  bl_mean_arm1, bl_mean_arm2, 
  bl_sd_arm1, bl_sd_arm2, 
  bl_n_arm1, bl_n_arm2, 
  data = depressionPsyCtr,
  study = study)

# Provide cluster variable: p-values are adjusted within clusters
metapsyTools:::testBaselineImbalance(
  bl_mean_arm1, bl_mean_arm2, 
  bl_sd_arm1, bl_sd_arm2, 
  bl_n_arm1, bl_n_arm2, 
  data = depressionPsyCtr,
  cluster = study,
  study = study)

# Change adjustment method
metapsyTools:::testBaselineImbalance(
  bl_mean_arm1, bl_mean_arm2, 
  bl_sd_arm1, bl_sd_arm2, 
  bl_n_arm1, bl_n_arm2, 
  data = depressionPsyCtr,
  cluster = study,
  study = study, p.adj = "bonferroni")
} # }
```
