# Meta-Regression method for objects of class 'runMetaAnalysis'

Serves as a wrapper for `metareg` or `update.rma`, depending on the
class of the fitted model.

## Usage

``` r
metaRegression(x, ...)
```

## Arguments

- x:

  A model extracted from an object of class `runMetaAnalysis`.

- ...:

  Additional arguments.

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Paula Kuper
<paula.r.kuper@gmail.com>, Pim Cuijpers <p.cuijpers@vu.nl>

## Examples

``` r
if (FALSE) { # \dontrun{
metaRegression(res$model.combined, ~ rob + scale(year))
} # }
```
