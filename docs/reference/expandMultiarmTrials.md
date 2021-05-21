# Expand multiarm trials

This is a function to expand the format of meta-analysis data for its
applicability to the `calculateEffectSizes` function of the
`metapsyTools`-package.

## Usage

``` r
expandMultiarmTrials(data,
                     vars.for.id = c("study", "primary",
                                     "Outc_measure",
                                     "Time", "Time_weeks",
                                     "sr_clinician"),
                     study.indicator = "study",
                     multiarm.indicator = "is.multiarm",
                     no.arms.indicator = "no.arms",
                     group.indicator = "condition",
                     condition.specification = "Cond_spec",
                     groups.column.indicator = c("_trt1", "_trt2"),
                     group.names = list("ig" = "ig",
                                        "cg" = "cg"),
                     data.format = NULL)
```

## Arguments

- data:

  Meta-analysis data stored as a `data.frame`.

- vars.for.id:

  `character` vector, containing column names of all variables used to
  construct unique comparison IDs.

- study.indicator:

  `character`, signifying the name of the variable containing the study
  name.

- multiarm.indicator:

  `numeric`, signifying if a row is part of a multiarm study (1) or not
  (0).

- no.arms.indicator:

  `character`, signifying the name of the variable containing the number
  of arms included in a study (typically 2).

- group.indicator:

  `character`, column name of the variable storing the study name.

- condition.specification:

  `character`, column name of the variable storing the trial condition
  name.

- groups.column.indicator:

  `character`. If the dataset is in wide format: a character vector with
  two elements, representing the suffix used to differentiate between
  the first and second treatment in a comparison.

- group.names:

  `list`, storing the name of the value corresponding to the
  intervention group (`"ig"`) and control group (`"cg"`).

- data.format:

  `character`. Either `"long"` or `"wide"`, depending on the format of
  the dataset in `data`. `NULL` by default, which lets the user define
  the format after the function has been called.

## Value

`expandMultiarmTrials` returns the meta-analysis data set as class
`data.frame` (if results are saved to a variable). The rows of multiarm
studies are expanded so that each intervention group has an
unambiguously assigned control group. It also generates the following
columns:

- `id` a *comparison*-specific ID variable.

- `study.id` a *study*-specific ID variable.

- `study` a study-specific variable containing the study name. For
  multiarm studies, this variable also specifies the active treatment
  indicated by `multiarm.group.indicator` behind the name of the study
  (e.g. `"Hauksson, 2017 -grp"`).

## Details

This function expands multiarm studies in a meta-analysis data set,
thereby ensuring that each comparison (intervention group vs. control
group in `condition`) is unique for a specific outcome, and thus has two
rows. For this purpose, it duplicates the corresponding row of the
control group condition if required. A specific study indicator variable
is created that enables further use, e.g. in 3-level models. For more
details see the help vignette:
[`vignette("metapsyTools")`](../articles/metapsyTools.md).

## See also

[`calculateEffectSizes`](calculateEffectSizes.md)

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Paula Kuper
<paula.r.kuper@gmail.com>, Pim Cuijpers <p.cuijpers@vu.nl>

## Examples

``` r
if (FALSE) { # \dontrun{
data("inpatients")
expandMultiarmTrials(inpatients)
} # }
```
