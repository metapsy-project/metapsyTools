# Changelog

## metapsyTools 1.0.15

Results of earlier analyses can change in the situations listed under
“Changed results”.

### Changed results

- [`subgroupAnalysis()`](../reference/subgroupAnalysis.md), multilevel
  models: the *p*-value tested if all subgroup means are zero, not if
  subgroups differ. It now tests subgroup differences.
- [`subgroupAnalysis()`](../reference/subgroupAnalysis.md), CHE and
  CCREM models: the analysis ignored the correlated sampling errors
  and/or the crossed random effect of the selected model. It now keeps
  them.
- [`subgroupAnalysis()`](../reference/subgroupAnalysis.md): effects
  excluded in the selected model (e.g. `"outliers"`, `"rob"`,
  `"lowest"`) were included again.
- [`imputeResponse()`](../reference/imputeResponse.md),
  `calculateEffectSizes(impute.response = TRUE)`: with
  `lower.is.better = FALSE` and no `cutoff`, non-responders were counted
  as responders.
- [`runMetaAnalysis()`](../reference/runMetaAnalysis.md),
  `es.measure = "EER"`/`"CER"`: an event rate of 0 was treated as 50%.
  Rates of 0 and 1 now use a continuity correction of 0.5.
- [`runMetaAnalysis()`](../reference/runMetaAnalysis.md),
  `es.measure = "ROM"`: estimates were back-transformed from values
  already rounded on the log scale.
- [`correctPublicationBias()`](../reference/correctPublicationBias.md):
  `lower.is.better` now defaults to `NULL`, which uses the sign of the
  pooled effect. The former default `TRUE` estimated the selection model
  in the wrong direction when effects were positive, e.g. log-risk
  ratios of response or positively coded *g* values, and could lead to
  extreme estimates. Explicit values that contradict the pooled effect
  now cause a warning.
- [`correctPublicationBias()`](../reference/correctPublicationBias.md),
  `es.measure = "EER"`/`"CER"`: the limit meta-analysis row showed the
  unadjusted estimate. Selection models are no longer estimated for
  event rates.
- [`correctPublicationBias()`](../reference/correctPublicationBias.md),
  `method.tau = "FE"`: trim-and-fill and the selection model used a
  random-effects model. Both now use a fixed-effect model.
- [`proportionMID()`](../reference/proportionMID.md): multilevel models
  now use the total heterogeneity variance instead of only the
  between-study variance; `method.tau = "FE"` models use the
  fixed-effect estimate; `mid` is log-transformed for ROM models.
- [`simulateTreatmentCycles()`](../reference/simulateTreatmentCycles.md):
  `avg.no.treatments` did not account for response rate decay.
- `calculateEffectSizes(change.sign = )`: missing values in the column
  turned `.g` into `NA`; they now mean no sign change.
- [`includeSwitchedArms()`](../reference/includeSwitchedArms.md): the
  sign of `.log_rom` is now reversed as well.
- [`rom.m.sd()`](../reference/rom.m.sd.md),
  [`rom.change.m.sd()`](../reference/rom.change.m.sd.md): means with
  opposite signs now return `NA` for both the effect and its standard
  error.
- NNT of `es.type = "raw"` in the `"waap.wls"` model is now based on the
  effects used in this model.

### Fixed errors

- [`plot.runMetaAnalysis()`](../reference/plot.runMetaAnalysis.md):
  fixed errors for `"ccrem"`, `"ccrem.che"` and, if it was the first
  model, `"threelevel.che"` and `"waap.wls"`; `plot(x)` now shows both
  models for `"lowest.highest"`.
- [`plot.runMetaAnalysis()`](../reference/plot.runMetaAnalysis.md):
  forest plots of multilevel models now show the same CI and prediction
  interval as the summary table.
- [`print()`](https://rdrr.io/r/base/print.html) for `runMetaAnalysis`,
  `subgroupAnalysis` and `proportionMID` objects cut off model names and
  truncated columns on consoles without color support or with limited
  width.
- [`runMetaAnalysis()`](../reference/runMetaAnalysis.md): fixed errors
  when `study.var` is not `"study"` (`"combined"` model), for
  `es.measure = "EER"`/`"CER"` with `use.rve = FALSE`, and for
  `vcov = "complex"` without `which.combine.var` in the data.
- [`runMetaAnalysis()`](../reference/runMetaAnalysis.md): prediction
  intervals of multilevel models were missing with `hakn = FALSE`.
- [`runMetaAnalysis()`](../reference/runMetaAnalysis.md): `$data` now
  contains the input columns; `.g` was overwritten for ROM, EER and CER
  models, which affected later
  [`rerun()`](../reference/Replacement-functions.md) calls.
- [`runMetaAnalysis()`](../reference/runMetaAnalysis.md): risk of bias
  data was not added to the `"lowest"` and `"highest"` models.
- `low.rob.filter` and `high.rob.filter` now accept all valid
  statements, e.g. `"rob>2"` or statements using several variables.
- [`summary.runMetaAnalysis()`](../reference/summary.runMetaAnalysis.md):
  fixed the missing correlation and aggregation level for the
  `"combined"` model, and errors for `"ccrem.che"`, `"waap.wls"`, event
  rates and some `method.tau = "FE"` models.
- [`correctPublicationBias()`](../reference/correctPublicationBias.md),
  [`createRobSummary()`](../reference/createRobSummary.md): `"lowest"`
  and `"highest"` can now be selected.
- [`correctPublicationBias()`](../reference/correctPublicationBias.md):
  the bootstrapped G-squared confidence interval was always empty; NNTs
  are now provided for pre-calculated risk ratios.
- [`eb()`](../reference/eb.md), [`blup()`](../reference/blup.md): fixed
  error for `es.type = "raw"`.
- [`createRobSummary()`](../reference/createRobSummary.md): now works
  with a single domain; errors are no longer hidden.
- [`metapsyNNT()`](../reference/metapsyNNT.md): fixed error for integer
  event counts.
- [`imputeResponse()`](../reference/imputeResponse.md): fixed the error
  handling of `rho`, and the return value of `method = "rci"` without
  control group data.
- `filterPoolingData(.filter.missing.rows = TRUE)`: default column is
  now `.g`.
- [`subgroupAnalysis()`](../reference/subgroupAnalysis.md): clear errors
  for `"lowest.highest"` and `"waap.wls"`; labels are correct if
  variable names are contained in other variable names.
- `i2.ci.threelevel<-` set an argument that no longer exists; it now
  sets `i2.ci.boot`. `i2.ci.boot<-` was added.
- [`print()`](https://rdrr.io/r/base/print.html) is now registered for
  `flagStudy` objects, and `flagEffectSizes` results with missing effect
  sizes can be used in
  [`runMetaAnalysis()`](../reference/runMetaAnalysis.md).
- [`calculateEffectSizes()`](../reference/calculateEffectSizes.md) now
  reports rows without any effect size, changed values in an existing
  `.g` column, and old change score column names (`change_m_arm1`, …).

### Internal changes and documentation

- Replaced
  [`clubSandwich::impute_covariance_matrix()`](http://jepusto.github.io/clubSandwich/reference/impute_covariance_matrix.md),
  which is deprecated, with an internal function that gives the same
  results.
- Corrected and extended the documentation, including effect size
  direction, NNTs, I-squared definitions, robust variance estimators and
  the `"combined"` model.
- Added tests that compare model results with independent calculations
  in [meta](https://github.com/guido-s/meta/) and
  [metafor](https://www.metafor-project.org).
