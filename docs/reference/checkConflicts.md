# Check for (potential) data format conflicts

This function checks for potential data formatting conflicts that may
produce errors or incorrect results when applying the
[`calculateEffectSizes`](calculateEffectSizes.md) or
[`runMetaAnalysis`](runMetaAnalysis.md) function later on.

## Usage

``` r
checkConflicts(.data,
               vars.for.id = c("study", "outcome_type",
                               "instrument", "time",
                               "time_weeks",
                               "rating"),
               .condition = "condition",
               .condition.specification = "multi",
               .groups.column.indicator = c("_arm1", "_arm2"))
```

## Arguments

- .data:

  Meta-analysis data stored as a `data.frame`, to be checked by the
  function.

- vars.for.id:

  `character` vector, containing column names of all variables used to
  construct unique comparison IDs.

- .condition:

  `character`. The prefix of the two variables in `data` in which the
  conditions (e.g. "guided iCBT", "waitlist") of the trial arm
  comparison are stored.

- .condition.specification:

  `character`, name of the column containing the specific condition in
  each trial arm. For multiarm trials, these conditions *must* be
  distinct (e.g. `"cbt-guided"` and `"cbt-unguided"`).

- .groups.column.indicator:

  `character`. A character vector with two elements, representing the
  suffix used to differentiate between the first and second arm in a
  comparison.

## Value

The type of data returned by `checkConflicts` depends on the outcome of
the evaluation. When no problems have been detected, the function simply
returns the data set provided in `.data`.

When (potential) formatting formatting issues have been detected, the
function throws a message and returns the affected studies/`data.frame`
columns. In particular, results are provided within a `list` of three
objects:

- `allConflicts`, a `data.frame` containing all affected rows,
  regardless of conflict type.

- `idConflicts`, a `data.frame` containing rows with ID/number of arms
  conflicts.

- `cgConflicts`, a `data.frame` containing rows with reference arm
  conflicts (there must be a unique control/reference group for each
  comparison).

The returned list has class `checkConflicts`.

## See also

[`checkDataFormat`](checkDataFormat.md)

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Paula Kuper
<paula.r.kuper@gmail.com>, Pim Cuijpers <p.cuijpers@vu.nl>

## Examples

``` r
if (FALSE) { # \dontrun{
data("depressionPsyCtr")

# Example 1: Use defaults and simply run checks
depressionPsyCtr %>%
  checkDataFormat() %>%
  checkConflicts() -> res

# Example 2: Overrule defaults; this will produce a conflict
depressionPsyCtr %>%
  checkDataFormat() %>%
  checkConflicts(vars.for.id = "study") -> res

} # }
```
