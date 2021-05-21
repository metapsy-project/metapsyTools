# Convert Metapsy database into network meta-analysis format

This function converts a database following the [Metapsy data
standard](https://docs.metapsy.org/data-preparation/format/) into format
that is suitable for network meta-analysis (e.g. using
[`netmeta::netmeta()`](https://rdrr.io/pkg/netmeta/man/netmeta.html)).

## Usage

``` r
convertNetmeta(
               
               # Continuous outcomes (endpoint scores)
               mean_arm1, mean_arm2, sd_arm1, sd_arm2, 
               n_arm1, n_arm2,
               
               # Continuous outcomes (change scores)
               mean_change_arm1, mean_change_arm2, 
               sd_change_arm1, sd_change_arm2,
               n_change_arm1, n_change_arm2, 
               
               # Response (event counts)
               event_arm1, event_arm2, totaln_arm1, 
               totaln_arm2, 
               
               # Study characteristics
               condition_arm1, condition_arm2, study, 
               ..., data = NULL)
```

## Arguments

- mean_arm1:

  Mean score in the first trial arm.

- mean_arm2:

  Mean score in the second trial arm.

- sd_arm1:

  Standard deviation in the first trial arm.

- sd_arm2:

  Standard deviation in the second trial arm.

- n_arm1:

  Sample size in the first trial arm.

- n_arm2:

  Sample size in the second trial arm.

- mean_change_arm1:

  Mean change scores in the first trial arm.

- mean_change_arm2:

  Mean change scores in the second trial arm.

- sd_change_arm1:

  Standard deviation of change scores in the first trial arm.

- sd_change_arm2:

  Standard deviation of change scores in the second trial arm.

- n_change_arm1:

  Sample size of change scores in the first trial arm.

- n_change_arm2:

  Sample size of change scores in the second trial arm.

- event_arm1:

  Number of responders in the first trial arm.

- event_arm2:

  Number of responders in the second trial arm.

- totaln_arm1:

  Total number of participants in the first trial arm.

- totaln_arm2:

  Total number of participants in the second trial arm.

- condition_arm1:

  Treatment or format in the first trial arm.

- condition_arm2:

  Treatment or format in the second trial arm.

- study:

  Study labels for each comparison.

- ...:

  Additional arguments. Can be used to specify additional columns to be
  included in the output (see Details).

- data:

  Dataset following the [Metapsy data
  standard](https://docs.metapsy.org/data-preparation/format/)
  (optional).

## Value

Returns a `data.frame` in wide-format, containing the calculated effect
sizes (standardized mean differences, SMDs) for each required
comparison. The following columns will be included in all outputs:

- `studlab`: Study label for each comparison. If response counts were
  used as outcome, "`(response)`" is appended to the study name.

- `treat1`: Condition or format used in the first trial arm.

- `treat2`: Condition or format used in the second trial arm.

- `TE`: The calculated effect size (SMD).

- `seTE`: Standard error of the calculated effect size.

Depending on whether continuous or binary outcomes (or both) were used,
the dataset will also include further columns containing the raw data
used to obtain the effect size:

- `n1`: Sample size in the first trial arm.

- `n2`: Sample size in the second trial arm.

- `mean1`: Mean (change) scores in the first trial arm.

- `mean2`: Mean (change) scores in the second trial arm.

- `sd1`: Standard deviation in the first trial arm.

- `sd2`: Standard deviation in the second trial arm.

- `event1`: Responders in the first trial arm.

- `event2`: Responders in the second trial arm.

Studies for which effect sizes could not be calculated will be saved as
a character vector in the `removed.studies` attribute. They can be
extracted using `attr(res, "removed.studies")`, where `res` is the
returned data frame.

## Details

This function converts a Metapsy database into a "wider" format dataset
that can be used to run network meta-analyses. Returned objects are
optimized for
[`netmeta::netmeta()`](https://rdrr.io/pkg/netmeta/man/netmeta.html) and
can be used "out-of-the box" in this package.

The function will perform an expansion of multi-arm trials, which is
required for most network meta-analysis implementations. Thus, the
function will calculate all three unique comparisons for three-arm
trials, all six comparisons for four-arm trials, etc.

Two additional formatting requirements must be met to conduct the
conversion:

- Each comparison in a trial is only allowed to provide exactly one
  effect size/contrast. This may be resolved by filtering the dataset
  beforehand using [filterPoolingData](filterPoolingData.md) or
  [filterPriorityRule](filterPriorityRule.md). The function will return
  an informative error message if non-unique comparisons are found.

- The function can only use raw continuous outcome and binary response
  data to calculate SMDs for each comparison. Rows with other effect
  size information (e.g. pre-calculated effects based on *t* or
  *F*-tests) will not be included. If comparisons had to be removed,
  affected studies will be printed into the console and should be added
  manually.

It is also possible to add additional columns (e.g., columns included in
the dataset provided in `data`) to the final data frame. These columns
have to be specified as additional arguments in the function call, where
the argument name will be used as the column name (see Examples).

## See also

[`netmeta::netmeta()`](https://rdrr.io/pkg/netmeta/man/netmeta.html),
[`checkDataFormat()`](checkDataFormat.md)

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Paula Kuper
<paula.r.kuper@gmail.com>, Pim Cuijpers <p.cuijpers@vu.nl>

## Examples

``` r
if (FALSE) {

# Filter database so that only unique comparisons remain for each study.
data <- depressionPsyCtr %>% 
  filterPriorityRule(
    instrument = c("ces-d", "phq-9", "scl", 
                   "hdrs", "bdi-2", "scid")) %>% 
  filterPoolingData(year >= 1985, study != "Barrett, 2001") 

# Convert endpoint, change score, and response data.
dat.netmeta <- convertNetmeta(
  
  # Continuous outcome data
  mean_arm1, mean_arm2, sd_arm1, sd_arm2, n_arm1, n_arm2,
  mean_change_arm1, mean_change_arm2, sd_change_arm1, 
  sd_change_arm2, n_change_arm1, n_change_arm2, 
  
  # Response data
  event_arm1 = event_arm1, event_arm2 = event_arm2, 
  totaln_arm1 = totaln_arm1, totaln_arm2 = totaln_arm2, 
  
  # Treatments to be used in NMA
  condition_arm1 = condition_arm1, 
  condition_arm2 = condition_arm2, 
  
  # Additional column to be added
  scale = instrument,
  
  # Study label and data
  study = study, data = data)

# Load netmeta and perform NMA
library(netmeta)
netmeta(TE, seTE, treat1, treat2, studlab, 
        data = dat.netmeta, reference.group = "wl")
        
        
# Example using metapsyData database
library(metapsyData)
d <- getData("depression-psyctr", version = "22.0.2")

d$data %>% 
  filterPriorityRule(instrument = c("phq-9", "ces-d", "hdrs")) %>% 
  filterPoolingData(!study %in%
                      c('Baumgartner, 2021', 'Brown, 1984', 'Fann, 2015', 
                        'Fledderus, 2012', 'Floyd, 2004', 'Kleiboer, 2015', 
                        'Lemma, 2013', 'Mohr, 2013', 'Nezu, 1989', 
                        'NystrÃ¶m, 2017', 'Pecheur, 1984', 'Propst, 1992', 
                        'Rehm, 1981', 'Rohan, 2007', 'Scogin, 1989', 'Selmi, 1990', 
                        'Smith, 2017a', 'Titov, 2010', 'Tomasino, 2017', 
                        'Watt, 2000', 'Westerhof, 2019', 'Araya, 2021', 
                        'Choi, 2014')) %>% 
  convertNetmeta(mean_arm1, mean_arm2, sd_arm1, sd_arm2, n_arm1, n_arm2,
                 event_arm1 = event_arm1, event_arm2 = event_arm2, 
                 totaln_arm1 = totaln_arm1, totaln_arm2 = totaln_arm2,
                 condition_arm1 = condition_arm1, condition_arm2 = condition_arm2,
                 study = study, format = format, data = .) -> dat.netmeta

# Extract studies for which no effect sizes could be calculated
attr(dat.netmeta, "removed.studies")

# Run network meta-analysis
netmeta(TE, seTE, treat1, treat2, studlab, data = dat.netmeta)


# Multi-arm expansion with single trials:
# - Using continuous outcome
convertNetmeta(mean_arm1 = c(4.12, 5.74),
               mean_arm2 = c(5.74, 6.41),
               sd_arm1 = c(4.22, 5.15),
               sd_arm2 = c(5.15, 2.79),
               n_arm1 = c(50, 50),
               n_arm2 = c(50, 50),
               condition_arm1 = c("cbt", "dyn"),
               condition_arm2 = c("dyn", "wl"),
               study = c("Doe, 1999", "Doe, 1999"))

# - using response outcome
convertNetmeta(event_arm1 = c(22, 12),
               event_arm2 = c(12, 5),
               totaln_arm1 = c(87, 89),
               totaln_arm2 = c(89, 92),
               condition_arm1 = c("cbt", "dyn"),
               condition_arm2 = c("dyn", "wl"),
               study = c("Doe, 1999", "Doe, 1999"))
               
# - using study format instead of conditions
format.data <- data.frame(study = c("Doe, 1999", "Doe, 1999", "Miller, 2000",
                                    "Willms, 2017", "Willms, 2017"),
                          format_arm1 = c("gsh", "ush", "ush", "gsh", "ush"),
                          format_arm2 = c("ush", "wl",  "cau", "ush", "cau"),
                          mean_arm1 = c(4.12, 5.74, 3.21, 4.99, 6.23),
                          mean_arm2 = c(5.74, 6.41, 6.29, 6.23, 6.41),
                          sd_arm1 = c(4.22, 5.15, 4.21, 4.00, 5.92),
                          sd_arm2 = c(5.15, 2.79, 4.52, 5.92, 3.12),
                          n_arm1 = c(50, 50, 76, 30, 30),
                          n_arm2 = c(50, 50, 75, 30, 30))

convertNetmeta(mean_arm1, mean_arm2, sd_arm1, sd_arm2, n_arm1, n_arm2,
               condition_arm1 = format_arm1, condition_arm2 = format_arm2,
               study = study, data = format.data)
}
```
