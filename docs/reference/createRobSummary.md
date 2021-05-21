# Create a summary risk of bias plot

If the `rob.data` argument has been specified, this function allows to
create a summary risk of bias plot for results of the
[`runMetaAnalysis`](runMetaAnalysis.md) function.

## Usage

``` r
createRobSummary(model, 
                 name.low, 
                 name.high, 
                 name.unclear, 
                 which.run = model$which.run[1])
```

## Arguments

- model:

  An object of class `runMetaAnalysis`, created by the
  [`runMetaAnalysis`](runMetaAnalysis.md) function.

- name.low:

  A `character` vector, specifying which code(s) have been used in the
  original data for studies with a low risk of bias.

- name.high:

  A `character` vector, specifying which code(s) have been used in the
  original data for studies with a high risk of bias.

- name.unclear:

  A `character` vector, specifying which code(s) have been used in the
  original data for studies with unclear risk of bias.

- which.run:

  The model in `model` that should be used for the summary risk of bias
  plot. Uses the default analysis in `model` if no value is specified by
  the user. Possible values are `"overall"`, `"combined"`, `"lowest"`,
  `"highest"`, `"outliers"`, `"influence"` and `"rob"`.

## Value

Creates a RevMan-type risk of bias summary plot.

## See also

[dmetar::rob.summary](http://dmetar.protectlab.org/reference/rob.summary.md),
[robvis::rob_summary](https://rdrr.io/pkg/robvis/man/rob_summary.html),
[meta::rob](https://rdrr.io/pkg/meta/man/rob.html)

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Paula Kuper
<paula.r.kuper@gmail.com>, Pim Cuijpers <p.cuijpers@vu.nl>

## Examples

``` r
if (FALSE) { # \dontrun{


# Define ROB data to be added to the models
robData = list(
  # Names of ROB variables included in 'data'
  domains = c("sg", "ac", "ba", "itt"),
  # Long-format labels for each ROB domain
  domain.names = c("Sequence Generation", 
                   "Allocation Concealment", 
                   "Blinding of Assessors", 
                   "ITT Analyses"),
  # Codes used to rate the risk of bias (sr=self-report)
  categories = c("0", "1", "sr"),
  # Symbols that should be used for these codes in forest plots
  symbols = c("-", "+", "s"),
  # Colors to be used in forest plots for each of these codes
  colors = c("red", "green", "yellow"))

# Run meta-analyses with ROB data
res <- depressionPsyCtr %>% 
  filterPoolingData(condition_arm1 %in% c("cbt", "pst", "3rd")) %>% 
  runMetaAnalysis(rob.data = robData)

# Create a summary plot
createRobSummary(res, 
                 name.low = "1", 
                 name.high = "0", 
                 name.unclear = "sr")

# Create a summary plot for the "combined" model
# - Recode 'sr' (self-report) as low risk of bias
createRobSummary(res, 
                 name.low = c("1", "sr"), 
                 name.high = "0", 
                 name.unclear = NULL,
                 which.run = "combined")
                 
                 
} # }
```
