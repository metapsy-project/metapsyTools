# Generate RoB rating for a Metapsy database

Based on a prepared extraction sheet, this function allows you to
automatically generate risk of bias judgments for an entire Metapsy
database, using the RoB-2 Metapsy extension.

## Usage

``` r
createRobRatings(rob.data, database = NULL)
```

## Arguments

- rob.data:

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

- database:

  A meta-analytic database conforming with the Metapsy data standard. In
  particular, the variables `rand_ratio`, `rand_arm1`/`rand_arm2`,
  `attr_arm1`/`attr_arm2` and `rating` must be populated according to
  the Metapsy [data
  standard](https://docs.metapsy.org/data-preparation/format/). This
  argument is `NULL` by default, which means that only domain-specific
  and overall RoB ratings are generated from the RoB extraction sheet
  provided to `rob.data`.

## Value

Returns an object of class `"createRobRatings"`. This includes:

- `rob.data`. The original extraction sheet, with signalling question
  responses based on the database information appended. If a study has
  more than one comparison in the database, ratings are displayed for
  each comparison in a study.

- `database`. The database with the RoB ratings appended on a domain
  level: `rob_d1`, `rob_d2`, etc. This dataset also includes the
  comparison-level (`rob`) and study-level (`rob_study_lvl`) overall RoB
  judgments. `NULL` if no database was provided.

- `miss.studies`. A vector of studies that appear in the database, but
  were not found in the extraction sheet. `NULL` if no database was
  provided.

- `miss.studies.rob`. A vector of studies that appear in the extraction
  sheet, but were not found in the database. `NULL` if no database was
  provided.

## See also

[`checkRobDiscrepancies`](checkRobDiscrepancies.md)

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Clara Miguel Sanz
<clara.miguelsanz@vu.nl>, Pim Cuijpers <p.cuijpers@vu.nl>

## Examples

``` r
if (FALSE) {
library(readxl)
library(xlsx)

# Get example database from metapsy.org/assets/files/data.xlsx
data <- read_excel("data.xlsx")

# Get example extraction sheet from metapsy.org/assets/files/rob_data.xlsx
rob_data <- read_excel("rob_data.xlsx")

# Create ratings
tmp <- metapsyTools:::createRobRatings(rob.data = rob_data,
                                       database = data)
tmp

# Show database
tmp$database

# Show extraction sheet
tmp$rob.data

# Save both files into a new MS Excel file
xlsx::write.xlsx(tmp$database, "data_rated.xlsx", 
                 sheetName = "database", showNA = FALSE)
xlsx::write.xlsx(tmp$rob.data, "data_rated.xlsx", 
                 sheetName = "rob", showNA = FALSE,
                 append = TRUE)
}
```
