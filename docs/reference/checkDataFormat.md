# Check data format

This function checks if a `data.frame` object conforms with the [Metapsy
data standard](https://docs.metapsy.org/data-preparation/format/).

## Usage

``` r
checkDataFormat(data,
                must.contain = c("study", "condition_arm1",
                                 "condition_arm2", 
                                 "multi_arm1", 
                                 "multi_arm2",
                                 "outcome_type", "instrument",
                                 "time", "time_weeks",
                                 "rating", "mean_arm1", "mean_arm2",
                                 "sd_arm1", "sd_arm2",
                                 "n_arm1", "n_arm2",
                                 "event_arm1", "event_arm2",
                                 "totaln_arm1", "totaln_arm2"),
                variable.class = list("study" = "character", 
                                      "condition_arm1" = "character",
                                      "condition_arm2" = "character", 
                                      "multi_arm1" = "character", 
                                      "multi_arm2" = "character",
                                      "outcome_type" = "character", 
                                      "instrument" = "character",
                                      "time" = "character", 
                                      "time_weeks" = "numeric",
                                      "rating" = "character", 
                                      "mean_arm1" = "numeric", 
                                      "mean_arm2" = "numeric",
                                      "sd_arm1" = "numeric", 
                                      "sd_arm2" = "numeric",
                                      "n_arm1" = "numeric", 
                                      "n_arm2" = "numeric",
                                      "event_arm1" = "numeric", 
                                      "event_arm2" = "numeric",
                                      "totaln_arm1" = "numeric", 
                                      "totaln_arm2" = "numeric"))
```

## Arguments

- data:

  A `data.frame` containint meta-analysis data.

- must.contain:

  `character` vector, containing all the variable names the data set
  should contain. Defaults correspond with the Metapsy data standard.

- variable.class:

  `list`, defining the required class for some or all variables. If the
  class differs in `data`, the function will try to convert the variable
  to the desired class. Defaults correspond with the Metapsy data
  standard.

## Value

`checkDataFormat` returns messages that specify if input variables,
values and classes of the variables are as defined. The output should
then be passed to [`checkConflicts`](checkConflicts.md) and then to
[`runMetaAnalysis`](runMetaAnalysis.md).

If default settings are used, `checkDataFormat()` can be used in
combination with [`checkConflicts`](checkConflicts.md) to determine if a
dataset follows the [Metapsy data
standard](https://docs.metapsy.org/data-preparation/format/). Datasets
that are formatted using this standard can be directly used in the
[analysis
module](https://tools.metapsy.org/articles/metapsytools#the-analysis-module)
in `metapsyTools`; for example the
[`runMetaAnalysis`](runMetaAnalysis.md) function.

## Details

The function checks if:

- the data set contains all relevant variables and

- variables have the desired class (if not, it tries to convert).

## See also

[`expandMultiarmTrials`](expandMultiarmTrials.md)

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Paula Kuper
<paula.r.kuper@gmail.com>, Pim Cuijpers <p.cuijpers@vu.nl>

## Examples

``` r
if (FALSE) { # \dontrun{
data("depressionPsyCtr")

# Example 1: Check with default arguments
checkDataFormat(depressionPsyCtr)

#Example 2: Check for non-default arguments
checkDataFormat(depressionPsyCtr,
                must.contain = c("study", "condition",
                                 "primary", "year"),
                variable.class = list(study = "character",
                                      no.arms = "numeric"))
} # }

```
