# Calculate the proportion of true effect sizes above a meaningful threshold

Based on results of the [`runMetaAnalysis`](runMetaAnalysis.md)
function, this function allows to estimate the proportion of true effect
sizes that exceed a user-defined meaningful (e.g. clinically relevant)
threshold.

## Usage

``` r
proportionMID(model, 
              mid = NULL, 
              which = "all", 
              test = "smaller", 
              plot = FALSE)
```

## Arguments

- model:

  A class `runMetaAnalysis` object, created by the
  [`runMetaAnalysis`](runMetaAnalysis.md) function.

- mid:

  A `numeric` value, indicating a clinically relevant effect threshold
  (e.g. a minimally important difference; \\MID\\; Cuijpers et al.,
  [2014](https://onlinelibrary.wiley.com/doi/full/10.1002/da.22249))
  that should be used to estimate the proportion of true effect sizes
  that exceed this cut-off. If the outcome measure used in `model` is
  Hedges' \\g\\, the provided value should also be a standardized mean
  difference. If the outcome measure of `model` is a risk ratio, the
  treshold should also be provided as an (untransformed) risk ratio.

- which:

  The model in `model` that should be used to estimate the proportions.
  Defaults to `"all"`, which means that proportions are calculated for
  all models. Alternatively, possible values are `"overall"`,
  `"combined"`, `"lowest"`, `"highest"`, `"outliers"`, `"influence"` and
  `"rob"`, if these models are available in the `model` object. If
  [`correctPublicationBias`](correctPublicationBias.md) has been run,
  `"trimfill"`, `limitmeta` and `selection` are also possible options.
  It is also possible to concatenate model names, meaning that
  proportions are calculated for all the supplied models.

- test:

  By default, the function estimates the proportion of true effects
  *below* the provided threshold in `mid` (`test="smaller"`).
  Alternatively, one can specify `test="bigger"`. This will calculate
  the proportion of true effects *above* the treshold.

- plot:

  Should a density plot illustrating the proportions be returned?
  Defaults to `FALSE`. Please note that an S3 `plot` method is available
  for outputs of this function even when `plot=FALSE` (see "Details").

## Value

Returns an object of class `"proportionMID"`. An S3 `plot` method is
defined for this object class, which allows to create a density plot
illustrating the estimated proportions, using the model-based estimate
of the pooled effect size and between-study heterogeneity \\\tau^2\\.

## Details

The `proportionMID` function implements an approach to estimate the
proportion of true effect sizes exceeding a (scientifically or
clinically) relevant threshold, as proposed by Mathur & VanderWeele
([2019](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8057)).
These estimated proportions have been suggested as a useful metric to
determine the impact that between-study heterogeneity in a meta-analysis
has on the "real-life" interpretation of results.

If, for example, a pooled effect is significant, high between-study
heterogeneity can still mean that a substantial proportion of true
effects in the studies population are practicially irrelevant, or even
negative. Conversely overall non-significant effects, in the face of
large heterogeneity, can still mean that a substantial proportion of
studies have non-negligible *true* effects.

As recommended by Mathur & VanderWeele
([2019](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8057)), the
`proportionMID` function also automatically calculates the proportion of
true effects exceeding the *"inverse"* of the user-defined effect (e.g.,
if `mid=-0.24`, by default, the function also estimates the proportion
of true effects that is larger than \\g\\=0.24; note the changed sign).
This can be used to check for e.g. clinically relevant negative effects.

When the `plot` method is used, or when `plot` is set to `TRUE` in the
function, a plot showing the assumed distribution of true effects based
on the estimated meta-analytic model is created. Notably, is is assumed
that the random-effects distribution of true effect sizes is
approximately normal. This simplifying assumption is required for this
(and many other meta-analytic methods; Jackson & White,
[2018](https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201800071))
to hold.

Confidence intervals provided by the functions are calculated using the
asymptotic closed-form solution derived using the Delta method in Mathur
& VanderWeele
([2019](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8057)).
Following their recommendations, a warning is printed when \\p\\\<0.15
or \\p\\\>0.85, since in this case the asymptotic CIs should be
interpreted cautiously; CIs based on boostrapping would be preferable in
this scenario and can be calculated using the
[`MetaUtility::prop_stronger()`](https://rdrr.io/pkg/MetaUtility/man/prop_stronger.html)
function.

## References

Cuijpers, P., Turner, E. H., Koole, S. L., Van Dijke, A., & Smit, F.
(2014). What is the threshold for a clinically relevant effect? The case
of major depressive disorders. *Depression and Anxiety, 31*(5), 374-378.

Jackson, D., & White, I. R. (2018). When should meta-analysis avoid
making hidden normality assumptions?. *Biometrical Journal, 60*(6),
1040-1058.

Mathur, M. B., & VanderWeele, T. J. (2019). New metrics for
meta-analyses of heterogeneous effects. *Statistics in Medicine, 38*(8),
1336-1342.

## See also

[`runMetaAnalysis`](runMetaAnalysis.md),
[`correctPublicationBias`](correctPublicationBias.md)

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Paula Kuper
<paula.r.kuper@gmail.com>, Pim Cuijpers <p.cuijpers@vu.nl>

## Examples

``` r
if (FALSE) { # \dontrun{
# Run meta-analysis; then estimate the proportion
depressionPsyCtr %>% 
  filterPoolingData(condition_arm1 == "cbt") %>% 
  runMetaAnalysis()  -> x

proportionMID(x, mid = -0.24)
proportionMID(x, mid = -0.24, "outliers") %>% plot()


# If bootstrap CIs are requested in runMetaAnalysis,
# calculation of CIs around p is possible for all available
# models.
depressionPsyCtr %>% 
  filterPoolingData(condition_arm1 == "cbt") %>% 
  runMetaAnalysis(i2.ci.boot = TRUE, nsim.boot = 1000) %>% 
  correctPublicationBias() -> x

proportionMID(x, mid = -0.33)

# Run meta-analysis based on RRs; then estimate proportion
# of true effects bigger than the defined threshold
depressionPsyCtr %>% 
  filterPoolingData(
    condition_arm1 %in% c("cbt", "pst", 
                          "dyn", "3rd wave")) %>% 
  runMetaAnalysis(es.measure = "RR") %>% 
  proportionMID(mid = 1.13, test = "bigger")
} # }
```
