# Display available SD-reference instruments

Prints the catalogue of patient or clinician-rated outcome measurement
instruments for which a reference standard deviation is available within
`metapsyTools`. The shorthand in the `instr` column is the lookup key
used by [`calculateEffectSizes`](calculateEffectSizes.md) when
`sd.reference = "fill"` or `"override"`: by default, this code is
matched (case-insensitively) against the `instrument` column of the
data.

## Usage

``` r
listSDReferences(pattern = NULL, print = TRUE)
```

## Arguments

- pattern:

  `character`. Optional regular-expression filter applied to both the
  `instr` shorthand and the `full` instrument name. Matching is
  case-insensitive. `NULL` (default) returns the entire catalogue.

- print:

  `logical`. If `TRUE` (default), the catalogue is pretty-printed to the
  console. The (possibly filtered) table is always returned invisibly so
  it can be captured into a variable.

## Value

Invisibly, a `data.frame` with one row per instrument and the columns
`instr`, `full`, `k`, `k_es`, `sd`, `lo`, `hi`, `tau2`, `i2`. The
`"version"`, `"version_date"`, `"source_url"` and `"citation"`
attributes of the bundled catalogue are preserved.

## Details

The bundled catalogue is derived from the Metapsy SD-reference database
and keeps only the headline random-effects pooled SD per instrument (no
subgroup, no sensitivity analysis). For subgroup- or
sensitivity-specific estimates, browse the full reference database at
[metapsy.org/database/sd-reference.html](https://metapsy.org/database/sd-reference.html).

## See also

[`calculateEffectSizes`](calculateEffectSizes.md)

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# Full catalogue
listSDReferences()

# Filter
listSDReferences("depression")
listSDReferences("^BDI")

# Capture without printing
tbl <- listSDReferences(print = FALSE)
} # }
```
