# The 'sdReference' dataset

Catalogue of patient or clinician-rated outcome measurement instruments
for which a reference standard deviation is available within
`metapsyTools`. The shorthand in the `instr` column is the lookup key
used by [`calculateEffectSizes`](calculateEffectSizes.md) when
`sd.reference = "fill"` or `"override"`: by default, this code is
matched (case-insensitively) against the `instrument` column of the
data.

## Usage

``` r
data("sdReference")
```

## Format

A `data.frame`.

## Author

Mathias Harrer, Bernát Fábián
