# Add information that varies between trial arms as extra columns to your meta-analysis dataset

Creates two additional columns for each selected variable in which
information is stored separately for the intervention and control group.
This is typically useful when trial-level variables (i.e. variables that
differ between trial arms) are to be included in the final meta-analysis
dataset.

## Usage

``` r
addTrialArmInfo(.data, ...,
            .group.indicator = "condition",
            .name.intervention.group = "ig",
            .name.control.group = "cg",
            .vars.for.id = c("study", "primary",
                             "Outc_measure",
                             "Time", "Time_weeks"))
```

## Arguments

- .data:

  A `data.frame` containing unique intervention-control group
  comparisons, as created by the
  [`expandMultiarmTrials`](expandMultiarmTrials.md) function.

- ...:

  \<[dplyr_data_masking](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>.
  The name of several columns (included in `.data`) that are trial-level
  variables to be added as columns to `.data`. To add multiple
  variables, simply separate them using a comma.

- .group.indicator:

  `character`. Name of the column in `.data` which encodes the
  intervention/control group rows.

- .name.intervention.group:

  `character`. Name used in the `.group.indicator` variable to identify
  the intervention group rows.

- .name.control.group:

  `character`. Name used in the `.group.indicator` variable to identify
  the control group rows.

- .vars.for.id:

  `character` vector, containing column names of all variables used to
  construct unique comparison IDs.

## Value

`addTrialArmInfo` returns a dataset as class `data.frame`. This dataset
contains all the information previously stored in `.data`, plus two
columns for each selected trial arm variable (one for the intervention
and one for the control group).

## Details

Before running the meta-analysis, it is necessary to select only the
rows containing calculated effect sizes ('`es`'). This results in an
information loss when data differs between trial arms within one study
(e.g. the sample size *n* is often not identical in both arms of a
study); only the row of the "active"/intervention arm is selected, and
the information of the control group arm is discarded.

`addTrialArmInfo` is a convenience function which allows to avoid this
information loss by adding trial-specific information as extra columns
in the dataset. Two columns are created for each feature: one containing
the value of the intervention arm, and another containing the
information in the control arm.

The function is only applicable to datasets with expanded multiarm
trial; that is, the output of
[`expandMultiarmTrials`](expandMultiarmTrials.md) (or
[`expandMultiarmTrials`](expandMultiarmTrials.md), followed by
[`calculateEffectSizes`](calculateEffectSizes.md)).

For more details see the help vignette:
[`vignette("metapsyTools")`](../articles/metapsyTools.md).

## See also

[`expandMultiarmTrials`](expandMultiarmTrials.md)

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Paula Kuper
<paula.r.kuper@gmail.com>, Pim Cuijpers <p.cuijpers@vu.nl>

## Examples

``` r
if (FALSE) { # \dontrun{

# Example 1: calculate effect sizes
# then add "Post_N" as trial arm variable
data("inpatients")
inpatients %>%
  checkDataFormat() %>%
  expandMultiarmTrials() %>%
  calculateEffectSizes() %>%
  addTrialArmInfo(Post_N) %>%
  filterPoolingData(primary == 1)

# Example 2: add several trial arm variables simultaneously
inpatients %>%
  checkDataFormat() %>%
  expandMultiarmTrials() %>%
  calculateEffectSizes() %>%
  addTrialArmInfo(Post_N, Rand_N, Cond_spec) %>%
  filterPoolingData(primary == 1)
} # }
```
