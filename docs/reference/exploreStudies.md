# Explore included treatments and comparisons

This function allows to summarize included treatments and treatment
comparisons in a data set.

## Usage

``` r
exploreStudies(data,
               which = c("treatments", "comparisons"),
               
               # Metapsy standard variables
               .study.var = "study",
               .condition = "condition",
               .condition.specification = "multi",
               .groups.column.indicator = c("_arm1", "_arm2"),
               .trt.indicator = "arm",
               .n.vars = c("n", "n_change", "totaln", "N"),
               
               # Output
               html = TRUE)
```

## Arguments

- data:

  `data.frame`. Effect size data in the wide format, as created by
  [`calculateEffectSizes`](calculateEffectSizes.md). For the other
  default settings to be applicable, the data set should follow the
  [Metapsy data
  standard](https://docs.metapsy.org/data-preparation/format/).
  Alternatively, one can also provide an `metapsyDatabase` object as
  returned by
  [`metapsyData::getData()`](https://rdrr.io/pkg/metapsyData/man/getData.html),
  or a meta-analysis object returned by
  [`runMetaAnalysis`](runMetaAnalysis.md).

- which:

  Should the data set be summarized with respect to the included
  treatments (`"treatments"`) or treatment comparisons
  (`"comparisons"`)? Defaults to `"treatments"`.

- .study.var:

  `character`. The name of the variable in the data set in which the
  study labels are stored.

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

- .n.vars:

  `character`. A character vector which includes the names of all
  variables in the data set in which sample size information is stored.
  Only the prefix is needed, where `.groups.column.indicator` provides
  the suffixes.

- html:

  `logical`. Should an HTML table be created for the results? Default is
  `TRUE`.

## Value

Returns an object of class `"exploreStudies"`. This object includes a
list object called `summary` in which the counts for distinct treatments
(`conditions`) and comparisons (`comparisons`) are summarized, as well
as a `data.frame` data. This data frame includes the initially provided
data set collapsed by study (so that each row represents one study). To
this data set, variables are added that encode how many arms with a
specific condition are included in the trial (e.g. if `cbt=2`, this
means that two CBT groups are included in the trial), as well as the
number of distinct comparisons, and the sample size of both (these
columns all start with `n.`). This can be helpful to perform further
descriptive analyses.

## Details

Using the variables provided in the `.n.vars` argument, `exploreStudies`
calculates the arm- and study-specific sample sizes. If no adequate
information is provided, sample sizes cannot be calculated for a study.
If this is the case, a warning is printed, pointing to the studies with
missing sample size information.

## See also

[`createStudyTable`](createStudyTable.md),
[`calculateEffectSizes`](calculateEffectSizes.md),
[`subgroupAnalysis`](subgroupAnalysis.md),
[`correctPublicationBias`](correctPublicationBias.md),
[`metaRegression`](metaRegression.md),
[`runMetaAnalysis`](runMetaAnalysis.md).

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Paula Kuper
<paula.r.kuper@gmail.com>, Pim Cuijpers <p.cuijpers@vu.nl>

## Examples

``` r
if (FALSE) { # \dontrun{
# Explore studies in built-in dataset
data("depressionPsyCtr")
exploreStudies(depressionPsyCtr, "treatments") 
exploreStudies(depressionPsyCtr, "comparisons") 

# - Extract metapsy database using metapsyData
# - Filter CBT and PST studies
# - Run a meta-analysis and explore synthesize studies
library(metapsyData)
getData("depression-psyctr", version="22.0.2") %>% 
  filterPoolingData(condition_arm1 %in% c("cbt", "pst")) %>% 
  runMetaAnalysis(which.run = c("combined")) -> res

exploreStudies(res)
exploreStudies(res, "comparisons")
} # }
```
