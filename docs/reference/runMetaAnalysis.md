# Run different types of meta-analyses

This wrapper function allows to simultaneously pool effect sizes using
different meta-analytic approaches.

## Usage

``` r
runMetaAnalysis(data,

                # Models to run
                which.run = c("overall", "combined",
                              "lowest.highest", "outliers",
                              "influence", "rob", "threelevel",
                              "threelevel.che"),
                              
                # Effect size measure
                es.measure = c("g", "RR", "EER", "CER"),
                es.type = c("precalculated", "raw"),
                es.var = ifelse(identical(es.measure[1], "RR"), 
                                          ".log_rr", ".g"),
                se.var = ifelse(identical(es.measure[1], "RR"), 
                                          ".log_rr_se", ".g_se"),
                es.binary.raw.vars = 
                  c(".event_arm1", ".event_arm2",
                    ".totaln_arm1", ".totaln_arm2"),
                    
                # Estimator of the heterogeneity variance
                method.tau = "REML",
                method.tau.ci = "Q-Profile",
                i2.ci.boot = FALSE,
                nsim.boot = 5e3,
                hakn = TRUE,
                
                # Data specifications
                study.var = "study",
                arm.var.1 = "condition_arm1",
                arm.var.2 = "condition_arm2",
                measure.var = "instrument",
                low.rob.filter = "rob > 2",
                round.digits = 2,
                rob.data = NULL,
                
                # Model specifications
                which.combine = c("arms", "studies"),
                which.combine.var = "multi_arm1",
                which.outliers = c("overall", "combined"),
                which.influence = c("overall", "combined"),
                which.rob = c("overall", "combined"),
                which.waap.wls = c("overall", "combined"),
                nnt.cer = 0.2,
                rho.within.study = 0.6,
                phi.within.study = 0.9,
                power.within.study = 0.8,
                w1.var = ifelse(identical(es.measure[1], "g"), 
                                "n_arm1", "totaln_arm1"),
                w2.var = ifelse(identical(es.measure[1], "g"), 
                                "n_arm2", "totaln_arm2"),
                time.var = "time_weeks",
                vcov = c("simple", "complex"),
                near.pd = FALSE,
                use.rve = TRUE,
                
                # Output
                html = TRUE,
                
                # Additional arguments
                ...)
```

## Arguments

- data:

  `data.frame`. Effect size data in the wide format, as created by
  [`calculateEffectSizes`](calculateEffectSizes.md).

- which.run:

  `character`. Selection of models to be calculated. See 'Details'.

- es.measure:

  `character`. Should meta-analyses be calculated using the
  bias-corrected standardized mean difference (`"g"`; default), risk
  ratios (`"RR"`), or logit-transformed experimental/control group event
  rates (`"EER"` or `"CER"`)? Meta-analyses will only be conducted using
  comparisons that contain non-missing values in the `es.var` and
  `se.var` columns.

- es.type:

  `character`. Should pre-calculated or raw event data (i.e. the
  Mantel-Haenszel method) be used for meta-analyses of risk ratios? Can
  be set to `"precalculated"` (default) or `"raw"`.

- es.var:

  `character`. Specifies the name of the variable containing the (pre-
  calculated) effect size data in `data`. When `es.measure = "g"`,
  `"EER"`, or `"CER"`, this is set to `.g` by default; `".log_rr"` is
  used when `es.measure = "RR"`. The default settings correspond with
  the standard output of
  [`calculateEffectSizes`](calculateEffectSizes.md).

- se.var:

  `character`. Specifies the name of the variable containing the
  (pre-calculated) standard errors (square root of the variance) of the
  effect size metric defined in `es.var`. If `es.measure = "g"`, this is
  automatically set to `.g_se`; if `es.measure = "RR"`, `".log_rr_se"`
  is used. The default settings correspond with the standard output of
  [`calculateEffectSizes`](calculateEffectSizes.md).

- es.binary.raw.vars:

  `character`. A vector defining the column names in `data` in which
  the (1) raw event counts in the experimental group, (2) raw event
  counts in the control/reference group, (3) sample size in the
  experimental group, and (4) sample size of the control/reference group
  are stored. Defaults correspond with the standard output of
  [`calculateEffectSizes`](calculateEffectSizes.md).

- method.tau:

  `character`. A character string indicating which method is used to
  estimate the between-study variance (tau-squared) and its square root
  (tau). Either `"REML"` (default), `"DL"`, `"PM"`, `"ML"`, `"HS"`,
  `"SJ"`, `"HE"`, or `"EB"`, can be abbreviated (see
  [`metagen`](https://rdrr.io/pkg/meta/man/metagen.html)). Use `"FE"` to
  use a fixed-effect/"common effect" model.

- method.tau.ci:

  `character`. A character string indicating which method is used to
  estimate the confidence interval of the between-study heterogeneity
  variance \\\tau^2\\. Either `"Q-Profile"` (default and recommended;
  [Viechtbauer,
  2017](http://www.wvbauer.com/lib/exe/fetch.php/articles:viechtbauer2007b.pdf)),
  `"BJ"`, `"J"`, or `"PL"` can be abbreviated. See
  [`metagen`](https://rdrr.io/pkg/meta/man/metagen.html) and
  [`rma.uni`](https://wviechtb.github.io/metafor/reference/rma.uni.html)
  for details.

- i2.ci.boot:

  `logical`. Confidence intervals for \\\tau^2\\ as calculated by the
  Q-Profile method are not directly applicable for three-level models,
  in which two heterogeneity variance components are estimated. By
  default, this argument is therefore set to `FALSE`, and no confidence
  intervals around \\\tau^2\\ and \\I^2\\ are provided for the
  `"threelevel"` and `"threelevel.che"` model. If this argument is set
  to `TRUE`, parametric bootstrapping is used to calculate confidence
  intervals around the between- and within-study heterogeneity estimates
  (\\\tau\\ and \\I^2\\). Please note that this can take several
  minutes, depending on the number of effect sizes. If
  [`correctPublicationBias`](correctPublicationBias.md) is used and
  `i2.ci.boot` is `TRUE`, bootstrapping will also be used to calculate
  confidence intervals around the \\G^2\\ statistic ([Rücker et al.,
  2011](https://academic.oup.com/biostatistics/article/12/1/122/391113))
  used in the limit meta-analysis (note that \\G^2\\ is printed as
  \\I^2\\ in this package).

- nsim.boot:

  `numeric` Number of bootstrap samples to be drawn when `i2.ci.boot` is
  `TRUE`. Defaults to 5000.

- hakn:

  `logical`. Should the Knapp-Hartung adjustment for effect size
  significance tests be used? Default is `TRUE`.

- study.var:

  `character`. The name of the variable in `data` in which the study IDs
  are stored.

- arm.var.1:

  `character`. The name of the variable in `data` in which the condition
  (e.g. "guided iCBT") of the *first* arm within a comparison are
  stored.

- arm.var.2:

  `character`. The name of the variable in `data` in which the condition
  (e.g. "wlc") of the *second* arm within a comparison are stored.

- measure.var:

  `character`. The name of the variable in `data` in which the
  instrument used for the comparison is stored.

- low.rob.filter:

  `character`. A filtering statement by which to include studies for the
  "low RoB only" analysis. Please note that the name of the variable
  must be included as a column in `data`.

- round.digits:

  `numeric`. Number of digits to round the (presented) results by.
  Default is `2`.

- rob.data:

  `list`. Optional list detailing how risk of bias data should be
  appended to the model (see Details).

- which.combine:

  `character`. Should multiple effect sizes within one study be pooled
  on an `"arms"` (default) or `"studies"` level? When a study is a
  multi-arm trial, setting `which.combine = "arms"` will aggregate the
  effect sizes for each trial arm individually before pooling; the
  `which.combine.var` argument can be used to control which effects
  within a study should be aggregated. When `which.combine = "studies"`,
  one overall aggregated effect is created for each study. This setting
  is preferable from a statistical perspective, since it ensures that
  all pooled effects can be assumed to be independent.

- which.combine.var:

  `character`. Additional grouping variable within studies to be used
  for the `"combined"` analysis when `which.combine = "arms"`. If the
  specified variable differs within one study (as defined by
  `study.var`), effects will be aggregated separately for each unique
  value in `which.combine.var`. Defaults to `"multi_arm1"`, the variable
  which encodes multi-arm intervention conditions in the [Metapsy data
  standard](https://docs.metapsy.org/data-preparation/format/#standard-variables).

- which.outliers:

  `character`. Which model should be used to conduct outlier analyses?
  Must be `"overall"` or `"combined"`, with `"overall"` being the
  default.

- which.influence:

  `character`. Which model should be used to conduct influence analyses?
  Must be `"overall"` or `"combined"`, with `"overall"` being the
  default.

- which.rob:

  `character`. Which model should be used to conduct the "low risk of
  bias only" analyses? Must be `"overall"` or `"combined"`, with
  `"overall"` being the default.

- which.waap.wls:

  `character`. Which model should be used to run the (non-default)
  WAAP-WLS model? Must be `"overall"` or `"combined"`, with `"overall"`
  being the default.

- nnt.cer:

  `numeric`. Value between 0 and 1, indicating the assumed control group
  event rate to be used for calculating NNTs via the Furukawa-Leucht
  method.

- rho.within.study:

  `numeric`. Value between 0 and 1, indicating the assumed correlation
  of effect sizes within studies. This is relevant to combine effect
  sizes for the `"combined"` analysis type, and used to estimate the
  variance-covariance matrices needed for the conditional and
  hierarchical effects three-level model. Default is `0.6`.

- phi.within.study:

  `numeric`. Value between 0 and 1, indicating the assumed one-week
  autocorrelation of effect sizes. This is only used when
  `vcov="complex"` to approximate the variance-covariance matrices
  needed for the `"combined"` and `"threelevel.che"` model. Default is
  0.9. See "Details".

- power.within.study:

  `numeric`. Value between 0 and 1, indicating the minimum power
  (\\1-\beta\\) required for a study to be considered sufficiently
  powered. This value is used to identify adequately powered studies in
  the WAAP-WLS analysis. See "Details".

- w1.var:

  `character`. Name of the variable in `data` in which the sample sizes
  of the first arm are stored. See "Details".

- w2.var:

  `character`. Name of the variable in `data` in which the sample sizes
  of the second arm are stored. See "Details".

- time.var:

  `character`. Name of the variable in `data` in which the assessment
  time point is stored. Should be expressed as weeks since
  randomization; other units (e.g. days, months) are also possible if
  `phi.within.study` is specified accordingly. See "Details".

- vcov:

  `character`. For the `"combined"` and `"threelevel.che"` model, should
  variance-covariance matrices (representing the dependency structure of
  the data) be approximated using a heterogeneous compound symmetry
  (`"simple"`; default) or unstructured matrix structure (`"complex"`)?
  See "Details".

- near.pd:

  `logical`. If at least one of the study variance-covariance matrices
  constructed when `vcov="complex"` is not positive definite/invertible,
  should the [`nearPD`](https://rdrr.io/pkg/Matrix/man/nearPD.html)
  function be used to compute the nearest positive definite matrix?
  Default is `FALSE`.

- use.rve:

  `logical`. Should robust variance estimation be used to calculate
  confidence intervals and tests of three-level models? `TRUE` by
  default.

- html:

  `logical`. Should an HTML table be created for the results? Default is
  `TRUE`.

- ...:

  Additional arguments.

## Value

Returns an object of class `"runMetaAnalysis"`. This object includes,
among other things, a `data.frame` with the name `summary`, in which all
results are summarized - including the studies which were removed for
some analysis steps. Other objects are the "raw" model objects returned
by all selected analysis types. This allows to conduct further
operations on some models specifically (e.g. run a meta-regression by
plugging one of the model objects in `meta:::update.meta`.

## Details

The `runMetaAnalysis` function is a wrapper for several types of
meta-analytic models that are typically used. It allows to run all of
these models in one step in order to generate results that are somewhat
closer to being "publication-ready". See
[`plot.runMetaAnalysis`](plot.runMetaAnalysis.md) and
[`createRobSummary`](createRobSummary.md) for additional functionality.

By *default*, the following models are calculated:

- `"overall"`. Runs a generic inverse-variance (random-effects) model.
  All included effect sizes are treated as independent. When
  `es.measure = "RR"` and `es.type = "raw"`, the Mantel-Haenszel method
  is used for pooling instead.

- `"combined"`. Pools all effect sizes within one study (defined by
  `study.var`) before pooling. This ensures that all effect sizes are
  independent (i.e., unit-of-analysis error & double-counting is
  avoided). To combine the effects, one has to assume a correlation of
  effect sizes within studies, empirical estimates of which are
  typically not available.

- `"lowest.highest"`. Runs a meta-analysis, but with only (i) the lowest
  and (ii) highest effect size within each study included.

- `"outlier"`. Runs a meta-analysis without statistical outliers (i.e.
  effect sizes for which the confidence interval does not overlap with
  the confidence intervall of the overall effect).

- `"influence"`. Runs a meta-analysis without influential cases (see
  [`influence.rma.uni`](https://wviechtb.github.io/metafor/reference/influence.rma.uni.html)
  for details).

- `"rob"`. Runs a meta-analysis with only low-RoB studies included.

- `"threelevel"`. Runs a multilevel (three-level) meta-analysis model,
  with effect sizes nested in studies.

- `"threelevel.che"`. Runs a multilevel (three-level) meta-analysis
  model, with effect sizes nested in studies. Variance-covariance
  matrices of each study with two or more effect sizes are estimated
  using `rho.within.study` as the assumed overall within-study
  correlation. This imputation allows to run a "correlated and
  hierarchical effects" (CHE) model, which is typically a good
  approximation for data sets with unknown and/or complex dependence
  structures.

The following *non-default* models can also be specified:

- `"ccrem"`. Runs a multilevel (three-level) cross-classified random
  effects model (CCREM), with `measure.var` modeled as a crossed random
  effect. Such models can be more appropriate when, across studies,
  effect sizes are nested in various instruments, subscales, or
  outcomes; and if we want to estimate (and generalize over) the
  variation due different instruments (or subscales, outcomes;
  Fernandez-Castilla et al.,
  [2018](https://doi.org/10.3758/s13428-018-1063-2)).

- `"ccrem.che"`. Runs a multilevel (three-level) cross-classified random
  effects model (CCREM), with `measure.var` modeled as a crossed random
  effect. Additionally, variance-covariance matrices of each study with
  two or more effect sizes are approximated using `rho.within.study` as
  the assumed overall within-study correlation, similar to the
  "correlated and hierarchical effects" (CHE, `threelevel.che`) model.

- `"waap.wls"`. Runs a weighted average of adequately powered studies
  (WAAP) weighted least squares model (WLS). If at least 3 effect sizes
  are available, effects are only computed among adequately powered
  studies, which are determined using the pooled common-effect estimate
  as basis. If less than 3 studies show adequate power, a simple WLS
  model among all included effect sizes is returned. This type of
  analysis employing multiplicative error models has been proposed to
  obtain more robust effect estimates that guard against common biases
  introduced by selective publication and QRPs (Stanley, Doucouliagos &
  Ioannidis,
  [2017](https://onlinelibrary.wiley.com/doi/10.1002/sim.7228); Stanley
  & Doucouliagos,
  [2017](https://onlinelibrary.wiley.com/doi/10.1002/jrsm.1211); Carter
  et al.,
  [2019](https://journals.sagepub.com/doi/full/10.1177/2515245919847196))

Internally, the `overall`, `combined`, `lowest.highest`, `outlier`,
`influence` and `rob` models are fitted by calling the
[`meta::metagen()`](https://rdrr.io/pkg/meta/man/metagen.html) or
[`meta::metabin()`](https://rdrr.io/pkg/meta/man/metabin.html) function,
respectively, in **{meta}** (Balduzzi, Rücker & Schwarzer,
[2019](https://pubmed.ncbi.nlm.nih.gov/31563865/)). The `threelevel`,
`threelevel.che`, `ccrem`, and `ccrem.che` models are implemented using
[`metafor::rma.mv()`](https://wviechtb.github.io/metafor/reference/rma.mv.html)
in **{metafor}** (Viechtbauer,
[2005](https://www.jstatsoft.org/article/view/v036i03)).

Outlier selection is implemented using the
[`dmetar::find.outliers()`](http://dmetar.protectlab.org/reference/find.outliers.md)
function, and influence analyses using the
[`dmetar::InfluenceAnalysis()`](http://dmetar.protectlab.org/reference/InfluenceAnalysis.md)
function. The latter function is a wrapper for
[`metafor::influence.rma.uni()`](https://wviechtb.github.io/metafor/reference/influence.rma.uni.html).

\\~\\

### Simple or complex variance-covariance approximation

The `vcov` argument controls if the effect size dependencies within the
data should be approximated using a `"simple"` (default) or more
`"complex"` (but potentially more accurate) method. This argument is
only relevant for the `"combined"`, `"threelevel.che"`, and
`"ccrem.che"` models. The default "simple" method constructs
variance-covariance matrices \\\Sigma_k\\ for each study using a
constant sampling correlation \\\rho\\ (defined by `rho.within.study`),
which is identical across all studies, outcomes, and time points. This
simplifying assumption is part of the formulation of the CHE model
originally provided by Pustejovsky and Tipton
([2022](https://link.springer.com/article/10.1007/s11121-021-01246-3)).

Naturally, employing a common value of \\\rho\\ across all studies may
not be reasonable in some analyses, and other information may be
available to better approximate the effect size dependencies in the
collected data. Setting `vcov` to `"complex"` allows to assume that
correlations between effect sizes may differ conditional on the type of
dependency. This means that the variance-covariance matrix \\\Sigma_k\\
of some study \\k\\ is approximated by an unstructured matrix with
varying \\\rho\_{ij}\\ (instead of a heterogeneous compound symmetry
matrix with fixed \\\rho\\, as is used when `vcov="simple"`).

\\\begin{array}{ccc}\texttt{vcov="simple"} & \texttt{vcov="complex"} &
\\ \Sigma_k = \begin{bmatrix} \sigma^2_1 \\ \rho \sigma_2 \sigma_1 &
\sigma^2_2 & & \\ \rho \sigma_3 \sigma_1 & \rho \sigma_3 \sigma_2 &
\sigma^2_3 & \\ \rho \sigma_4 \sigma_1 & \rho \sigma_4 \sigma_2 & \rho
\sigma_4 \sigma_3 & \sigma^2_4 \end{bmatrix} & \Sigma_k =
\begin{bmatrix} \sigma^2_1 & & & \\ \rho\_{21} \sigma_2 \sigma_1 &
\sigma^2_2 & & \\ \rho\_{31} \sigma_3 \sigma_1 & \rho\_{32} \sigma_3
\sigma_2 & \sigma^2_3 & \\ \rho\_{41} \sigma_4 \sigma_1 & \rho\_{42}
\sigma_4 \sigma_2 & \rho\_{43} \sigma_4 \sigma_3 & \sigma^2_4
\end{bmatrix} & \end{array}\\

For example, setting `vcov = "complex"` allows to additionally
incorporate assumed correlations specific to multiple testing over time
(e.g. correlations between effects at post-test and long-term
follow-up). The value provided in `phi.within.study` represents the
(auto-)correlation coefficient \\\phi\\, which serves as a rough
estimate of the re-test correlation after 1 week. When a vector of
follow-up lengths is provided in `time.var`, this allows to model the
gradual decrease in correlation between measurements over time.
Furthermore, it is possible to calculate a correlation coefficient
\\\rho_w\\ for multi-arm trials, which is directly proportional to the
size of each individual trial arm. When all trial arms have the same
size, meaning that each arm's weight \\w\\ is identical, \\\rho_w\\ is
known to be 0.5. Multiarm weights \\w\\ (and thus \\\rho_w\\) can be
derived if the `w1.var` and `w2.var` variables, containing the sample
size of each study arm, are provided.

Using the complex approximation method increases the risk that at least
one studies' \\\Sigma_k\\ matrix is not positive definite. In this case,
the function automatically switches back to the constant sampling
correlation approximation.

\\~\\

### Replacement functions

Once a model has been fitted using `runMetaAnalysis`, **replacement
functions** are defined for each function argument. This allows to
quickly tweak one or more analysis settings, which are implemented once
the `rerun` function is called. Say that we saved the results of
`runMetaAnalysis` in an object `m`. If, for example, we want to check
the results using a different estimator of \\\tau^2\\, leaving all other
settings the same, we could run e.g. `method.tau(m) <- "PM"`, followed
by `rerun(m)`. This would provide results using the Paule-Mandel
estimator. A list of all available setting replacement functions is
provided
[here](https://tools.metapsy.org/reference/replacement-functions).

\\~\\

### Risk of Bias data

Using the `rob.data` argument, it is possible to specify which variables
in the data set provided in `data` contain Risk of Bias (ROB) assessment
information. If specified, this allows to generate forest plots with
added ROB information when using
[`plot.runMetaAnalysis`](plot.runMetaAnalysis.md). ROB information added
this way can also be used to create ROB summary plots using
[`createRobSummary`](createRobSummary.md).

The object provided to `rob.data` must be a `list` element with these
elements:

- `domains`: A vector of characters, specifying variables that contain
  ratings for different ROB domains (e.g. allocation concealment,
  missing data handling, ...).

- `domain.names` (*Optional*): A vector of characters with the same
  length as `domains`, which provides long-format labels for each
  included domain.

- `overall.rob` (*Optional*): A single character specifying the variable
  in `data` that contains the overall ROB rating.

- `categories`: A vector of characters, specifying how ROB ratings were
  coded in the selected variables (e.g., `"low`, `"high"`, `"unclear"`).

- `symbols` (*Optional*): A vector of single-letter characters (or
  symbols) that should be used when plotting ROB rating in the forest
  plot.

- `colors` (*Optional*): A vector of characters of the same length as
  `categories`, specifying colors for each rating code.

For concrete usage examples with this functionality, see "Examples".

For more details see the [Get
Started](https://tools.metapsy.org/articles/metapsytools) vignette.

## See also

[`plot.runMetaAnalysis`](plot.runMetaAnalysis.md),
[`summary.runMetaAnalysis`](summary.runMetaAnalysis.md),
[`eb.runMetaAnalysis`](eb.runMetaAnalysis.md),
[`profile.runMetaAnalysis`](profile.runMetaAnalysis.md),
[`subgroupAnalysis`](subgroupAnalysis.md).
[`correctPublicationBias`](correctPublicationBias.md),
[`metaRegression`](metaRegression.md),
[`metaRegression`](metaRegression.md),
[`calculateEffectSizes`](calculateEffectSizes.md),

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Paula Kuper
<paula.r.kuper@gmail.com>, Pim Cuijpers <p.cuijpers@vu.nl>

## Examples

``` r
if (FALSE) { # \dontrun{
data("depressionPsyCtr")
library(meta)

depressionPsyCtr %>%
  checkDataFormat() %>%
  checkConflicts() %>%
  calculateEffectSizes() %>% 
  filterPoolingData(condition_arm2 %in% 
                      c("wl", "other ctr")) -> data

# Run the meta-analyses
runMetaAnalysis(data) -> res

# Check if variance components of "threelevel.che" model
# are identifiable
profile(res, "threelevel.che")

# Run non-default models
# - Cross-classified random effects model
runMetaAnalysis(data, 
                which.run = c("ccrem", "ccrem.che"), 
                vcov = "complex")

# - Weighted average of adequately powered studies (WAAP-WLS)
runMetaAnalysis(data, which.run = "waap.wls")
runMetaAnalysis(data, which.run = "waap.wls", power.within.study = .9)

# - Meta-analysis using raw response rate data
runMetaAnalysis(data, 
                es.measure = "RR",
                es.type = "raw")

# - Estimate intervention response rates, then pool 
data %>% 
  calculateEffectSizes(impute.response = TRUE) %>% 
  runMetaAnalysis(which.run = "combined", 
                  es.measure = "EER")

# Use replacement function to show results for
# differing settings
method.tau(res) <- "PM"
hakn(res) <- FALSE
rerun(res)

# Show forest plot (by default, "overall" is used)
plot(res)

# Comapare effects across models
plot(res, "summary")

# Show forest plot of specific analysis
plot(res, "outliers")
plot(res, "threelevel")
plot(res, "baujat")
plot(res, "influence")

# Show forest plot with empirical Bayes estimates
plot(res, "combined", eb = TRUE)

# Extract specific model and do further calculations
# (e.g. meta-regression on 'year')
metaRegression(res$model.overall, ~ scale(year))

# Conduct a subgroup analysis
subgroupAnalysis(res, country)

# Correct for publication bias/small-study effects
(correctPublicationBias(res) -> res.pb)
plot(res.pb, "trimfill")
plot(res.pb, "limitmeta")
plot(res.pb, "selection")

# For the combined analysis, set which.combine to
# "studies" here, so that all effects in a study are aggregated
# first before pooling
data %>% 
  runMetaAnalysis(which.combine = "studies") %>% 
  plot("combined")

# Define ROB data to be added to the models
robData <- list(
  # Names of ROB variables included in 'data'
  domains = c("sg", "ac", "ba", "itt"),
  # Long-format labels for each ROB domain
  domain.names = c("Sequence Generation", 
                   "Allocation Concealment", 
                   "Blinding of Assessors", 
                   "ITT Analyses"),
  # Codes used to rate the risk of bias (sr=self-report)
  categories = c("0", "1", "sr"),
  # Symbols that should be used for these codes in forest plots
  symbols = c("-", "+", "s"),
  # Colors to be used in forest plots for each of these codes
  colors = c("red", "green", "yellow"))

# Re-run model with appended ROB data
res <- runMetaAnalysis(data, rob.data = robData) 

# Generate forest plot with ROB data
plot(res, "combined")

# Create a summary plot
createRobSummary(res, 
                 name.low = "1", 
                 name.high = "0", 
                 name.unclear = "sr",
                 which.run = "combined")
} # }
```
