# Check for inconsistencies between two RoB extraction sheets

Based on prepared extraction sheets by two independent raters, this
function allows you to automatically check for differing studies and/or
inconsistencies in the ratings. Cohen's kappa is calculated
automatically, differentiating between "Yes" vs. "No"/"No Information"
ratings to obtain the base rate of agreement.

## Usage

``` r
checkRobDiscrepancies(data.1, data.2)
```

## Arguments

- data.1:

  An RoB extraction sheet. Columns of this file must have the same name
  as in the [RoB extraction sheet
  template](https://www.metapsy.org/assets/files/rob-template.xlsx)
  provided by the Metapsy initative. If the Metapsy template has been
  used, make sure to delete the top rows before importing, so that only
  the metapsyTools variables remain as the column names. Required
  columns are: `study`, `d1_1`, `d1_2`, `d1_3`, `d1_4`, `d1_notes`,
  `d2_5`, `d2_6`, `d2_7`, `d2_8`, `d2_9`, `d2_notes`, `d3_10`, `d3_11`,
  `d3_12`, `d3_13`, `d3_14`, `d3_notes`, `d4_15`, `d4_16`, `d4_17`,
  `d4_18`, `d4_notes`, `d5_19`, `d5_20`, `d5_21`, `d5_22`, `d5_23`,
  `d5_24`, `d5_notes`.

- data.2:

  The same extraction sheet from another (i.e., second) rater.

## Value

If discrepancies are found, the function will return a data frame with
the respective study/studies, along with the diverging ratings
(`discrepancies` element). If different studies are included in both
sheets, they will be saved under `diff.studies`.

## See also

[`createRobRatings`](createRobRatings.md)

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Clara Miguel Sanz
<clara.miguelsanz@vu.nl>, Pim Cuijpers <p.cuijpers@vu.nl>

## Examples

``` r
if (FALSE) { # \dontrun{
library(readxl)

# Get example extraction sheet from metapsy.org/assets/files/rob_data.xlsx
rob_data <- read_excel("rob_data.xlsx")

# Create second sheet with partly different ratings
rob_data -> rob_data_2
rob_data_2[-1,] -> rob_data_2
rob_data_2[1,"d2_5"] = "NI"
rob_data_2[1,"d4_15"] = "No/PN"

# Check for discrepancies
tmp <- metapsyTools:::checkRobDiscrepancies(rob_data, rob_data_2)
tmp
} # }
```
