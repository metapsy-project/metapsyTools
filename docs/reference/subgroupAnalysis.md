# Run subgroup analyses

This function allows to simultaneously conduct different subgroup
analyses using `runMetaAnalysis` objects.

## Usage

``` r
subgroupAnalysis(.model,
                 ...,
                 .which.run = .model$which.run[1],
                 .round.digits = 2,
                 .nnt.cer = NULL,
                 .tau.common = FALSE,
                 .html = TRUE)
```

## Arguments

- .model:

  An object of class `"runMetaAnalysis"`, created by
  [`runMetaAnalysis`](runMetaAnalysis.md).

- ...:

  \<[dplyr_data_masking](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>.
  A number of subgroup variables included in the original dataset
  provided to [`runMetaAnalysis`](runMetaAnalysis.md), separated by
  commas.

- .which.run:

  The model in `.model` that should be used for the subgroup analyses.
  Uses the default analysis in `.model` if no value is specified by the
  user.

- .round.digits:

  `numeric`. Number of digits to round the (presented) results by.
  Default is `2`.

- .nnt.cer:

  `numeric`. Value between 0 and 1, indicating the assumed control group
  event rate to be used for calculating NNTs via the Furukawa-Leucht
  method. If set to `NULL` (default), the value saved in `.model` is
  (re-)used.

- .tau.common:

  `logical`. Should a common (`TRUE`) or subgroup-specific (`FALSE`)
  estimate of the between-study heterogeneity be calculated when
  analyzing the subgroups? `FALSE` by default. Note that subgroup
  analyses based on "multilevel" models automatically assume common
  heterogeneity estimates.

- .html:

  `logical`. Should an HTML table be created for the results? Default is
  `TRUE`.

## Value

Returns an object of class `"subgroupAnalysis"`. This object includes,
among other things, a `data.frame` with the name `summary`, in which all
subgroup analysis results are summarized. Other objects are the "raw"
subgroup analysis model objects returned. This allows to conduct further
operations on some subgroup analysis specifically.

## Details

For more details see the [Get
Started](https://tools.metapsy.org/articles/metapsytools) vignette.

## See also

[`runMetaAnalysis`](runMetaAnalysis.md)

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
  calculateEffectSizes() %>% 
  filterPoolingData(condition_arm2 %in% 
                      c("wl", "other ctr")) -> data

# Run the meta-analyses
runMetaAnalysis(data) -> res

# Subgroup analysis
subgroupAnalysis(res, condition_arm2, country,
                 .which.run = "combined",
                 .tau.common = TRUE) -> sg
plot(sg, "condition_arm2")
plot(sg, "country")
} # }
```
