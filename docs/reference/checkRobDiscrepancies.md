# Check for inconsistencies between two RoB extraction sheets

Based on prepared extraction sheets by two independent raters, this
function automatically checks for differing studies and/or
inconsistencies in the ratings. Discrepancies are reported at the
original 3-category level (`"Yes/PY"`, `"No/PN"`, `"NI"`). Inter-rater
agreement is quantified using Cohen's \\\kappa\\, computed on a
binarised version of the ratings (`"Yes/PY"` vs. not-`"Yes/PY"`) over
pairwise complete cells.

## Usage

``` r
checkRobDiscrepancies(data.1, data.2)
```

## Arguments

- data.1:

  An RoB extraction sheet. Columns of this file must have the same name
  as in the [RoB extraction sheet
  template](https://www.metapsy.org/assets/files/rob-template.xlsx)
  provided by the Metapsy initiative. If the Metapsy template has been
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

A list containing (depending on what is found):

- `discrepancies`: a data frame with studies that differ between raters
  and the diverging ratings (returned only if discrepancies exist);

- `diff.studies`: studies included in only one of the two sheets
  (returned only if such studies exist);

- `kappa`: Cohen's \\\kappa\\ on the binarised ratings (`"Yes/PY"` vs.
  not) over pairwise complete cells.

## Details

**Discrepancy detection.** Ratings are compared cell-by-cell at the
original 3-category level. A discrepancy is flagged whenever the two
raters assigned different categories, or one rater provided a valid
rating and the other did not.

**Inter-rater agreement.** Cohen's \\\kappa\\ is computed after
binarising the ratings into `"Yes/PY"` vs. not-`"Yes/PY"` (the latter
combining `"No/PN"` and `"NI"`). Only cells in which both raters
provided a valid rating in {`"Yes/PY"`, `"No/PN"`, `"NI"`} are included
(pairwise complete cases). Let \\n\\ denote the number of such cells and
\\y_k^{(r)} = 1\\ if rater \\r\\ rated cell \\k\\ as `"Yes/PY"`, 0
otherwise. The observed agreement is \$\$p_0 = \frac{1}{n}
\sum\_{k=1}^{n} \mathbb{I}\\\left\[y_k^{(1)} = y_k^{(2)}\right\],\$\$
the marginal proportions of `"Yes/PY"` ratings are \$\$\pi_r =
\frac{1}{n} \sum\_{k=1}^{n} y_k^{(r)}, \quad r \in \\1, 2\\,\$\$ and the
chance-expected agreement is \$\$p_e = \pi_1 \pi_2 + (1 - \pi_1)(1 -
\pi_2).\$\$ Cohen's \\\kappa\\ is then \$\$\kappa = \frac{p_0 - p_e}{1 -
p_e}.\$\$

## See also

`createRobRatings`

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
