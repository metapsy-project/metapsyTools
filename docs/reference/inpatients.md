# The 'inpatients' dataset with 32 clinical trials

An example dataset containing data of 32 clinical trials. Its format
equals the format of data in the "metapsy" database that is similar to
the long format in *R*. The dataset also contains columns with study
characteristics that are important for effect size calculation and more.

## Usage

``` r
data("inpatients")
```

## Format

A `data.frame` with 179 rows and 46 variables:

- study:

  `character`, The study label containing the author(s) and year of the
  study.

- condition:

  `character`, The condition of the groups, either "intervention group"
  or "control group".

- Cond_spec:

  `character`, The specific intervention in conditions.

- is.multiarm:

  `numeric`, The dichotomized indication if a study has multiple arms or
  not.

- no.arms:

  `numeric`, The number of arms of a study.

- multiple.arms:

  `character`, The specification of arms in multiarm studies.

- Outc_type:

  `character`, The type of outcome.

- primary:

  `numeric`, The indication if a outcome is primary or not.

- Outc_measure:

  `character`, The outcome measure used.

- Time:

  `character`, The dichotomized time of assessment, either "post" or
  "FU".

- Time_weeks:

  `character`, The assessment time of FU in weeks.

- year:

  `numeric`, The year of the study.

- country:

  `character`, The country of the study.

- ...:

## Author

Mathias Harrer, Paula Kuper, Pim Cuijpers
