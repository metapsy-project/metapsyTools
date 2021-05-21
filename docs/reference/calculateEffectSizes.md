# Calculate effect sizes

This is a function to calculate effect sizes of meta-analysis data
prepared in the [Metapsy data
format](https://docs.metapsy.org/data-preparation/format/).

## Usage

``` r
calculateEffectSizes(data,
                     funcs.g = list(g.m.sd = g.m.sd,
                                    g.change.m.sd = g.change.m.sd,
                                    g.binary = g.binary,
                                    g.precalc = g.precalc),
                     funcs.rr = list(rr.binary = rr.binary,
                                     rr.precalc = rr.precalc),
                     include.switched.arms = FALSE,
                     change.sign = NULL,
                     impute.response = FALSE,
                     vars.for.id = c("study", "outcome_type",
                                     "instrument", "time",
                                     "time_weeks",
                                     "rating"),
                     .condition = "condition",
                     .condition.specification = "multi",
                     .groups.column.indicator = c("_arm1", "_arm2"),
                     .trt.indicator = "arm",
                     .impute.response.vars = c(m.trt.pre = "baseline_m_arm1", 
                                               m.trt.post = "mean_arm1", 
                                               sd.trt.post = "sd_arm1", 
                                               n.trt = "n_arm1",
                                               m.ctr.pre = "baseline_m_arm2", 
                                               m.ctr.post = "mean_arm2", 
                                               sd.ctr.post = "sd_arm2", 
                                               n.ctr = "n_arm2"))
```

## Arguments

- data:

  Meta-analysis data set formatted using the [Metapsy
  guidelines](https://tools.metapsy.org/articles/metapsyTools.html#required-data-structure).

- funcs.g:

  `list` of functions. These functions will be used to calculate the
  effect sizes (Hedges' *g*) based on the raw data (see Details).

- funcs.rr:

  `list` of functions. These functions will be used to calculate risk
  ratios based on the raw event data (see Details).

- include.switched.arms:

  `logical`. Should all unique arm *comparisons* (in lieu of unique arm
  *combinations*) be calculated? Default is `FALSE`.

- change.sign:

  `character`. Name of a `logical` column in `data`, encoding if the
  sign of a calculated effect size should be reversed (`TRUE`) or not
  (`FALSE`). Set to `NULL` (default) if no changes should be made.

- impute.response:

  `logical`. When calculating the (log)-risk ratios, should response
  rates be computed using the [`imputeResponse`](imputeResponse.md)
  function? `FALSE` by default. If defined, the column specified in
  `change.sign` will also be considered when calculating the response.

- vars.for.id:

  `character` vector, containing column names of all variables used to
  construct unique comparison IDs.

- .condition:

  `character`. The prefix of the two variables in `data` in which the
  conditions (e.g. "guided iCBT", "waitlist") of the trial arm
  comparison are stored.

- .condition.specification:

  `character`. The prefix of the two variables in the dataset which
  provide a "specification" of the trial arm condition in multiarm
  trials.

- .groups.column.indicator:

  `character`. A character vector with two elements, representing the
  suffix used to differentiate between the first and second arm in a
  comparison.

- .trt.indicator:

  `character`. A character specifying the name used to indicate the
  treatment arm.

- .impute.response.vars:

  `list`. Named list with the names of columns in `data` and the
  specific argument in [`imputeResponse`](imputeResponse.md) they should
  be used for.

## Value

`calculateEffectSizes` returns the meta-analysis data set as class
`data.frame` in wide format (if results are saved to a variable). It
also generates the following columns, wich are added to the data:

- `.id`: Unique identifier for a trial arm comparison/row.

- `.g`: Calculated effect size (Hedges' *g*).

- `.g_se`: Standard error of Hedges' *g*.

- `.log_rr`: Calculated effect size (*logRR*).

- `.log_rr_se`: Standard error of *logRR*.

- `.event_arm1`: Number of events (responders, remission, deterioration
  cases) in the first trial arm.

- `.event_arm2`: Number of events (responders, remission, deterioration
  cases) in the second trial arm.

- `.totaln_arm1`: Total sample size in the first trial arm.

- `.totaln_arm2`: Total sample size in the second trial arm.

## Details

By default, `calculateEffectSizes` calculates the small-sample bias
corrected standardized mean difference (Hedges' *g*) and log-risk
ratios, as well their respective standard errors, if adequate raw effect
size data is available for a comparison.

It is essential that the data set in `data` contains a few required
columns for this to work. An overview of the required data format is
provided on the ["Get
Started"](https://tools.metapsy.org/articles/metapsytools) page of the
`metapsyTools` documentation.

**Standardized mean differences (Hedges' *g*)** can be calculated from
the following column types:

- \(1\) Continuous Outcome Data

  - **`mean_arm1`**: Mean of the outcome in the first arm at the
    measured time point.

  - **`mean_arm2`**: Mean of the outcome in the second arm at the
    measured time point.

  - **`sd_arm1`**: Standard deviation of the outcome in the first arm at
    the measured time point.

  - **`sd_arm2`**: Standard deviation of the outcome in the second arm
    at the measured time point.

  - **`n_arm1`**: Sample size in the first trial arm.

  - **`n_arm2`**: Sample size in the second trial arm.

- \(2\) Change Score Data

  - **`mean_change_arm1`**: Mean score change between baseline and the
    measured time point in the first arm.

  - **`mean_change_arm2`**: Mean score change between baseline and the
    measured time point in the second arm.

  - **`sd_change_arm1`**: Standard deviation of the mean change in the
    first arm.

  - **`sd_change_arm2`**: Standard deviation of the mean change in the
    second arm.

  - **`n_change_arm1`**: Sample size in the first trial arm.

  - **`n_change_arm2`**: Sample size in the second trial arm.

- \(3\) Dichotomous Outcome Data

  - **`event_arm1`**: Number of events (responders, remission,
    deterioration cases) in the first trial arm.

  - **`event_arm2`**: Number of events (responders, remission,
    deterioration cases) in the second trial arm.

  - **`totaln_arm1`**: Sample size in the first trial arm.

  - **`totaln_arm2`**: Sample size in the second trial arm.

- \(4\) Pre-calculated Hedges' *g*

  - **`precalc_g`**: The pre-calculated value of Hedges' *g*
    (small-sample bias corrected standardized mean difference; [Hedges,
    1981](https://journals.sagepub.com/doi/10.3102/10769986006002107)).

  - **`precalc_g_se`**: Standard error of *g*.

The **log-risk ratio** and its standard error can be calculated from the
followin column types:

- \(1\) Dichotomous Outcome Data

  - **`event_arm1`**: Number of events (responders, remission,
    deterioration cases) in the first trial arm.

  - **`event_arm2`**: Number of events (responders, remission,
    deterioration cases) in the second trial arm.

  - **`totaln_arm1`**: Sample size in the first trial arm.

  - **`totaln_arm2`**: Sample size in the second trial arm.

- \(2\) Pre-calculated log-risk ratio

  - **`precalc_log_rr`**: The pre-calculated value of the log-risk ratio
    logRR, comparing events in the first arm to events in the second
    arm.

  - **`precalc_log_rr_se`**: The standard error of the log-risk ratio
    logRR, comparing events in the first arm to events in the second
    arm.

Other functions can be added to the list provided to `funcs.g` and
`funcs.rr`. However, results of the function must result in a
`data.frame` that contains the following columns:

- `.id`: Unique identifier for a trial arm comparison/row.

- `.g`: Calculated effect size (Hedges' *g*).

- `.g_se`: Standard error of Hedges' *g*.

- `.log_rr`: Calculated effect size (logRR).

- `.log_rr_se`: Standard error of logRR.

- `.event_arm1`: Number of events (responders, remission, deterioration
  cases) in the first trial arm.

- `.event_arm2`: Number of events (responders, remission, deterioration
  cases) in the second trial arm.

- `.totaln_arm1`: Total sample size in the first trial arm.

- `.totaln_arm2`: Total sample size in the second trial arm.

It is possible to set one or several of these column entries to `NA`;
but the columns themselves must be included.

For more details see the [Get
Started](https://tools.metapsy.org/articles/metapsytools) vignette.

## See also

[`checkDataFormat`](checkDataFormat.md),
[`imputeResponse`](imputeResponse.md)

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Paula Kuper
<paula.r.kuper@gmail.com>, Pim Cuijpers <p.cuijpers@vu.nl>

## Examples

``` r
if (FALSE) { # \dontrun{
data("depressionPsyCtr")
depressionPsyCtr %>%
    checkDataFormat() %>%
    checkConflicts() %>% 
    calculateEffectSizes()
} # }
```
