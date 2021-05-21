# Create study table

This function creates an overview table containing selected study
information.

## Usage

``` r
createStudyTable(.data, ...,
                 .round.by.digits = NULL,
                 .column.names = NULL,
                 .na.replace = "nr",
                 .html = TRUE)
```

## Arguments

- .data:

  Meta-analysis study data; typically data created by the
  [`expandMultiarmTrials`](expandMultiarmTrials.md) or
  [`calculateEffectSizes`](calculateEffectSizes.md) function. Trial
  arm-specific information (e.g. sample size in each group) can be added
  via [`addTrialArmInfo`](addTrialArmInfo.md). See 'Details'.

- ...:

  \<[dplyr_data_masking](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>.
  The name of several columns (included in `.data`) that should be added
  to the study table. Also allows to alter individual values/factor
  labels within a variable. See 'Details'.

- .round.by.digits:

  named `list`. Should contain the number of digits by which to round a
  numeric column in `.data`. The name of the column must be specified in
  the list element's name. Set to `NULL` if no rounding should be
  performed (default).

- .column.names:

  named `list`. If variable names should be renamed when producing the
  study table, the new name should be included in this list. The
  original column name must be specified as the name of the list
  element. Set to `NULL` if no renaming should be performed (default).

- .na.replace:

  `character` to replace `NA` values with; `"nr"` by default.

- .html:

  `logical`. Should an HTML table be produced? `TRUE` by default. See
  'Details'.

## Value

Returns a `data.frame` with all the selected variables. If `.html` is
`TRUE`, an HTML table will also be produced.

## Details

**General Purpose**: This function allows to select variables to be
included in a study table. Such study tables are typically part of a
meta-analysis report/article. Variables are included by adding their
names to the function call and separating them with commas. The columns
will appear in the exact same order as specified in the function.

**Trial-Arm Variables**: Before producing the final table,
`createStudyTable` will filter out all redundant rows based on the
selected variables. If you want to include information that differs
between the (two or more) trial arms (e.g. the sample size of each
group) as separate columns, you have to use
[`addTrialArmInfo`](addTrialArmInfo.md) first. This ensures that
individual columns are created for both the intervention and control
group (e.g. `N_ig` and `N_cg`), which can then be included in the call
to `createStudyTable`.

**Changing Values**: The function also allows to change specified values
within a variable. Factor levels encoded as numbers, for example (e.g.
`country = 1` for European studies, and so forth) can be changed by
adding a concatenated ([`c`](https://rdrr.io/r/base/c.html)) vector to
the name of the variable. This vector should contain the *new* value as
a `character` on the *left* side, and the *old* value on the *right*
side, separated by '`=`' (e.g. `country = c("Europe" = "1")`). The
values will then be recoded before producing the table.

**HTML Table**: By default, `createStudyTable` produces an HTML table
using [`kable`](https://rdrr.io/pkg/knitr/man/kable.html). The HTML
table makes copy & paste easier, particularly when working with MS Word,
since the table formatting is kept.

For more details see the help vignette:
[`vignette("metapsyTools")`](../articles/metapsyTools.md).

## See also

[`addTrialArmInfo`](addTrialArmInfo.md),
[`kable_styling`](https://rdrr.io/pkg/kableExtra/man/kable_styling.html)

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Paula Kuper
<paula.r.kuper@gmail.com>, Pim Cuijpers <p.cuijpers@vu.nl>

## Examples

``` r
if (FALSE) { # \dontrun{
# Filter out all primary outcomes, check data,
# calculate effect sizes, then produce study table 
# using selected information.
data("depressionPsyCtr")

depressionPsyCtr %>%
 filterPriorityRule(
   instrument = c("phq-9", "bdi-1", 
                  "hdrs", "ces-d")) %>%
 checkDataFormat() %>%
 checkConflicts() %>%
 calculateEffectSizes() %>%
 createStudyTable(
   study,
   diagnosis = c("Cutoff" = "3", "Mood" = "2",
                 "MDD" = "1"),
   age_group, instrument,
   mean_age, percent_women,
   condition_arm1 = c("CBT" = "cbt", "PST" = "pst",
                      "BA" = "bat", "LR" = "lrt",
                      "PDT" = "dyn", "IPT" = "ipt"),
   condition_arm2, n_arm1, n_arm2,
   country = c("Canada" = "4", "Europe" = "3", "USA" = "1",
               "Middle East" = "7"),
   sg, ac, ba, itt,
   .round.by.digits = list(mean_age = 0, n_arm1 = 0, 
                           n_arm2 = 0),
   .column.names = list(age_group = "age group",
                        n_arm1 = "N (arm1)",
                        n_arm2 = "N (arm2)",
                        percent_women = "% female")) -> table
} # }
```
