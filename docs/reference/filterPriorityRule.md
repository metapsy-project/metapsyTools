# Filter data based on a priority rule

This function filters rows of a dataset based on a priority rule for
specific variables defined by the user.

## Usage

``` r
filterPriorityRule(.data, ..., .study.indicator = "study")
```

## Arguments

- .data:

  A `data.frame` containing the calculated effect sizes, as created by
  the [`calculateEffectSizes`](calculateEffectSizes.md) function.

- ...:

  \<[dplyr_data_masking](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>.
  A number of prioritized filtering rules for variables. Should follow
  the form `variable = c("prio1", "prio2", ...)`. To apply multiple
  priority filters, simply separate them using a comma. For each study,
  rows are then selected based on the specified hierarchy for a
  variable. The priorities are provided as a concatenated vector,
  representing the variable levels. The level to appear first in this
  vector has the highest priority, the second one the second-largest
  priority, and so on. If a study contains none of the variable levels
  specified in the function call, the study is omitted entirely.

- .study.indicator:

  `character`. Name of the variable in which the study IDs are stored.

## Value

`filterPriorityRule` returns the filtered data set as class
`data.frame`. The filtered data set should then be ready for
meta-analytic pooling, for example using
[metagen](https://rdrr.io/pkg/meta/man/metagen.html). Further filters
can be applied using [`filterPoolingData`](filterPoolingData.md).

## See also

[`filterPoolingData`](filterPoolingData.md)

For more details see the [Get
Started](https://tools.metapsy.org/articles/metapsytools) vignette.

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Paula Kuper
<paula.r.kuper@gmail.com>, Pim Cuijpers <p.cuijpers@vu.nl>

## Examples

``` r
if (FALSE) { # \dontrun{
# Load data and calculate effect size
data("depressionPsyCtr")
depressionPsyCtr %>%
  checkDataFormat() %>%
  checkConflicts() %>%
  calculateEffectSizes() -> data

# Filter using four priority rules
filterPriorityRule(data,
                   condition_arm1 = c("cbt", "pst"),
                   condition_arm2 = c("cau", "wl", "cbt"),
                   instrument = c("cesd", "phq-9", "scl", "hdrs"),
                   time = c("post", "fu")) -> res
} # }
```
