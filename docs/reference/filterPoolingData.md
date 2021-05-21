# Filter data to be pooled for meta-analysis

This convenience function allows to create a filtered data set, which is
then ready to be used for meta-analytic pooling.

## Usage

``` r
filterPoolingData(.data, ...,
                  .filter.missing.rows = FALSE,
                  .es.column = es)
```

## Arguments

- .data:

  A `data.frame` containing the calculated effect sizes, as created by
  the [`calculateEffectSizes`](calculateEffectSizes.md) function.

- ...:

  \<[dplyr_data_masking](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>.
  A number of filtering statements (using variables in `.data`) that
  return a logical value. To apply multiple filters, simply separate
  them using a comma. "OR" statements can be provided using the `|`
  operator. Multiple filter statements separated using commas will be
  combined using the AND (`&`) operator. See "Details".

- .filter.missing.rows:

  `logical`. Should rows with no effect sizes be filtered out? Default
  is `FALSE`.

- .es.column:

  Name of the column in `.data` to be used for filtering out rows with
  no effect sizes. Default is `es`.

## Value

`filterPoolingData` returns the filtered data set as class `data.frame`.
The filtered data set should then be ready for meta-analytic pooling,
for example using [`runMetaAnalysis`](runMetaAnalysis.md).

## Details

The `filterPoolingData` function allows to apply several filters to your
meta-analysis data set all at once. When used in a pipe (`%>%`), you
only need to supply several filtering statements separated by commas,
using the same column names as they appear in the data set (e.g.
`primary == 1, meanage == 58`). The filtering statements are then
connected using "AND" (`&`).

If you want to apply an "OR" filter, simply use `|` instead of a comma
(e.g. `type == "cbt" | format == 6`). To select all rows that contain
one of several values in a variable, use `%in%`; e.g.
`study %in% c("Bailey, 2017", "Barth 2005")`.

The [`Detect`](Detect.md) function can be used within the function call
to search for variable elements which *contain* one or several selected
words (separated by `|`). To include all rows which contain the word
"cbt" or "wl" or "cau" in the "Cond_spec_trt2" variable, we can use
`Detect(Cond_spect_trt2, "cbt|wl|cau")`. This will also filter out
elements like `"cbt (online)"`, because "cbt" is included.

For more details see the [Get
Started](https://tools.metapsy.org/articles/metapsytools) vignette.

## See also

[`filterPriorityRule`](filterPriorityRule.md)

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Paula Kuper
<paula.r.kuper@gmail.com>, Pim Cuijpers <p.cuijpers@vu.nl>

## Examples

``` r
if (FALSE) { # \dontrun{

# Example 1: calculate effect sizes and then use multiple AND filters.
data("depressionPsyCtr")
depressionPsyCtr %>%
  calculateEffectSizes() %>%
  filterPoolingData(time == "post", instrument == "hdrs")

# Example 2: use OR filter
data("depressionPsyCtr")
depressionPsyCtr %>%
  calculateEffectSizes() %>%
  filterPoolingData(time == "post" | instrument == "hdrs")

# Example 3: use %in% operator
data("depressionPsyCtr")
depressionPsyCtr %>%
  calculateEffectSizes() %>%
  filterPoolingData(instrument %in% c("hdrs", "phq-9"))

# Example 4: Search for studies using "fuzzy-ish" matching
data("depressionPsyCtr")
depressionPsyCtr %>%
  calculateEffectSizes() %>%
  filterPoolingData(Detect(instrument, "bdi"))
} # }
```
