# metapsyTools 1.0.15

Results of earlier analyses can change in the situations listed under "Changed results".

## Changed results

* `subgroupAnalysis()`, multilevel models: the _p_-value tested if all subgroup means are zero, not if subgroups differ. It now tests subgroup differences.
* `subgroupAnalysis()`, CHE and CCREM models: the analysis ignored the correlated sampling errors and/or the crossed random effect of the selected model. It now keeps them.
* `subgroupAnalysis()`: effects excluded in the selected model (e.g. `"outliers"`, `"rob"`, `"lowest"`) were included again.
* `imputeResponse()`, `calculateEffectSizes(impute.response = TRUE)`: with `lower.is.better = FALSE` and no `cutoff`, non-responders were counted as responders.
* `runMetaAnalysis()`, `es.measure = "EER"`/`"CER"`: an event rate of 0 was treated as 50%. Rates of 0 and 1 now use a continuity correction of 0.5.
* `runMetaAnalysis()`, `es.measure = "ROM"`: estimates were back-transformed from values already rounded on the log scale.
* `correctPublicationBias()`: `lower.is.better` now defaults to `NULL`, which uses the sign of the pooled effect. The former default `TRUE` estimated the selection model in the wrong direction when effects were positive, e.g. log-risk ratios of response or positively coded _g_ values, and could lead to extreme estimates. Explicit values that contradict the pooled effect now cause a warning.
* `correctPublicationBias()`, `es.measure = "EER"`/`"CER"`: the limit meta-analysis row showed the unadjusted estimate. Selection models are no longer estimated for event rates.
* `correctPublicationBias()`, `method.tau = "FE"`: trim-and-fill and the selection model used a random-effects model. Both now use a fixed-effect model.
* `proportionMID()`: multilevel models now use the total heterogeneity variance instead of only the between-study variance; `method.tau = "FE"` models use the fixed-effect estimate; `mid` is log-transformed for ROM models.
* `simulateTreatmentCycles()`: `avg.no.treatments` did not account for response rate decay.
* `calculateEffectSizes(change.sign = )`: missing values in the column turned `.g` into `NA`; they now mean no sign change.
* `includeSwitchedArms()`: the sign of `.log_rom` is now reversed as well.
* `rom.m.sd()`, `rom.change.m.sd()`: means with opposite signs now return `NA` for both the effect and its standard error.
* NNT of `es.type = "raw"` in the `"waap.wls"` model is now based on the effects used in this model.

## Fixed errors

* `plot.runMetaAnalysis()`: fixed errors for `"ccrem"`, `"ccrem.che"` and, if it was the first model, `"threelevel.che"` and `"waap.wls"`; `plot(x)` now shows both models for `"lowest.highest"`.
* `plot.runMetaAnalysis()`: forest plots of multilevel models now show the same CI and prediction interval as the summary table.
* `print()` for `runMetaAnalysis`, `subgroupAnalysis` and `proportionMID` objects cut off model names and truncated columns on consoles without color support or with limited width.
* `runMetaAnalysis()`: fixed errors when `study.var` is not `"study"` (`"combined"` model), for `es.measure = "EER"`/`"CER"` with `use.rve = FALSE`, and for `vcov = "complex"` without `which.combine.var` in the data.
* `runMetaAnalysis()`: prediction intervals of multilevel models were missing with `hakn = FALSE`.
* `runMetaAnalysis()`: `$data` now contains the input columns; `.g` was overwritten for ROM, EER and CER models, which affected later `rerun()` calls.
* `runMetaAnalysis()`: risk of bias data was not added to the `"lowest"` and `"highest"` models.
* `low.rob.filter` and `high.rob.filter` now accept all valid statements, e.g. `"rob>2"` or statements using several variables.
* `summary.runMetaAnalysis()`: fixed the missing correlation and aggregation level for the `"combined"` model, and errors for `"ccrem.che"`, `"waap.wls"`, event rates and some `method.tau = "FE"` models.
* `correctPublicationBias()`, `createRobSummary()`: `"lowest"` and `"highest"` can now be selected.
* `correctPublicationBias()`: the bootstrapped G-squared confidence interval was always empty; NNTs are now provided for pre-calculated risk ratios.
* `eb()`, `blup()`: fixed error for `es.type = "raw"`.
* `createRobSummary()`: now works with a single domain; errors are no longer hidden.
* `metapsyNNT()`: fixed error for integer event counts.
* `imputeResponse()`: fixed the error handling of `rho`, and the return value of `method = "rci"` without control group data.
* `filterPoolingData(.filter.missing.rows = TRUE)`: default column is now `.g`.
* `subgroupAnalysis()`: clear errors for `"lowest.highest"` and `"waap.wls"`; labels are correct if variable names are contained in other variable names.
* `i2.ci.threelevel<-` set an argument that no longer exists; it now sets `i2.ci.boot`. `i2.ci.boot<-` was added.
* `print()` is now registered for `flagStudy` objects, and `flagEffectSizes` results with missing effect sizes can be used in `runMetaAnalysis()`.
* `calculateEffectSizes()` now reports rows without any effect size, changed values in an existing `.g` column, and old change score column names (`change_m_arm1`, ...).

## Internal changes and documentation

* Replaced `clubSandwich::impute_covariance_matrix()`, which is deprecated, with an internal function that gives the same results.
* Corrected and extended the documentation, including effect size direction, NNTs, I-squared definitions, robust variance estimators and the `"combined"` model.
* Added tests that compare model results with independent calculations in `{meta}` and `{metafor}`.
