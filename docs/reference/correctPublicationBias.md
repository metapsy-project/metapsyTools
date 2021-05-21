# Correct the effect size for publication bias/small-study effects

This function allows to add effect sizes estimates corrected for
publication bias/ small-study effects to results of the
`runMetaAnalysis` function.

## Usage

``` r
correctPublicationBias(model, 
                       which.run = model$which.run[1],
                       lower.is.better = TRUE,
                       selmodel.steps = 0.05,
                       ...)
```

## Arguments

- model:

  An object of class `runMetaAnalysis`, created by the `runMetaAnalysis`
  function.

- which.run:

  The model in `model` that should be used for the publication bias
  analyses. Uses the default analysis in `model` if no value is
  specified by the user. Possible values are `"overall"`, `"combined"`,
  `"lowest"`, `"highest"`, `"outliers"`, `"influence"` and `"rob"`.

- lower.is.better:

  Do lower values indicate better outcomes (i.e. higher effects)?
  Default is `TRUE`.

- selmodel.steps:

  Thresholds to be assumed for the step function in the selection model.
  Must be a vector of numbers referring to the cut-points in the
  selection models. If two-sided testing is assumed for the included
  studies, the cut-point must be doubled to obtain the assumed *p*-value
  (e.g. `selmodel.steps = c(0.03, 0.05)` means that *p*=0.06 and
  *p*=0.10 are assumed as selection thresholds). The default is `0.05`.

- ...:

  Additional arguments. See
  [trimfill.default](https://rdrr.io/pkg/meta/man/trimfill.html) and
  [limitmeta](https://rdrr.io/pkg/metasens/man/limitmeta.html).

## Value

Returns an object of class `"runMetaAnalysis"` and
`"correctPublicationBias"`. This object includes all original objects
included in `model`, but adds a `list` object with the name
`correctPublicationBias`. This list object includes all three fitted
publication bias analysis models, as well as the generated results.

## Details

The `correctPublicationBias` function is a wrapper running three
meta-analytic methods to control the pooled effect size for publication
bias and/or small-study effects:

- `"trimfill"`. Applies Duval and Tweedie's (2000a, 2000b) trim-and-fill
  algorithm, using the
  [trimfill](https://wviechtb.github.io/metafor/reference/trimfill.html)
  method in the `meta` package.

- `"limitmeta"`. Runs a limit meta-analysis as described in Rücker et
  al. (2011), using the implementation in the `limitmeta` package.

- `"selection"`. Runs a step function selection model using the
  [selmodel](https://wviechtb.github.io/metafor/reference/selmodel.html)
  function in `metafor`. For details see e.g. Vevea and Hedges (1995).

## References

Duval S & Tweedie R (2000a): A nonparametric "Trim and Fill" method of
accounting for publication bias in meta-analysis. *Journal of the
American Statistical Association, 95*, 89–98

Duval S & Tweedie R (2000b): Trim and Fill: A simple funnel-plot-based
method of testing and adjusting for publication bias in meta-analysis.
*Biometrics, 56*, 455–63

Rücker G, Schwarzer G, Carpenter JR, Binder H, Schumacher M (2011):
Treatment-effect estimates adjusted for small-study effects via a limit
meta-analysis. *Biostatistics, 12*, 122–42

Vevea, J. L., & Hedges, L. V. (1995). A general linear model for
estimating effect size in the presence of publication bias.
*Psychometrika, 60*(3), 419–435. ⁠https://doi.org/10.1007/BF02294384⁠

## See also

[`runMetaAnalysis`](runMetaAnalysis.md)

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Paula Kuper
<paula.r.kuper@gmail.com>, Pim Cuijpers <p.cuijpers@vu.nl>

## Examples

``` r
if (FALSE) { # \dontrun{
data("depressionPsyCtr")
depressionPsyCtr %>%
  checkDataFormat() %>%
  checkConflicts() %>%
  calculateEffectSizes() %>%
  filterPoolingData(condition_arm1 %in% c("cbt", "pst")) %>%
  runMetaAnalysis() -> res

# Correct for small-study-effects/publication bias
res %>% correctPublicationBias()
# Use additional arguments to control settings of the trim-and-fill
# and limit meta-analysis
correctPublicationBias(res,
                       which.run = "combined",
                       type = "R",
                       method.adjust = "mulim") 
 
# Generate plots
correctPublicationBias(res) %>% plot("trimfill")
correctPublicationBias(res) %>% plot("limitmeta")
correctPublicationBias(res) %>% plot("selection")

# Returned object is of class "runMetaAnalysis"; therefore,
# all S3 methods are available:
res <- res %>% correctPublicationBias()
metaRegression(res$model.overall, ~country)
} # }
```
