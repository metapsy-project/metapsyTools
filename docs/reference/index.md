# Package index

## Preparation Module

Functions to prepare the data and calculate effect sizes.

- [`checkDataFormat()`](checkDataFormat.md) : Check data format
- [`checkConflicts()`](checkConflicts.md) : Check for (potential) data
  format conflicts
- [`filterPoolingData()`](filterPoolingData.md) : Filter data to be
  pooled for meta-analysis
- [`filterPriorityRule()`](filterPriorityRule.md) : Filter data based on
  a priority rule
- [`calculateEffectSizes()`](calculateEffectSizes.md) : Calculate effect
  sizes
- [`flagEffectSizes()`](flagEffectSizes.md) : Flag extreme and/or
  implausible effect sizes

## Analysis Module

Functions to run meta-analyses.

- [`runMetaAnalysis()`](runMetaAnalysis.md) : Run different types of
  meta-analyses
- [`correctPublicationBias()`](correctPublicationBias.md) : Correct the
  effect size for publication bias/small-study effects
- [`subgroupAnalysis()`](subgroupAnalysis.md) : Run subgroup analyses
- [`metaRegression()`](metaRegression.md) : Meta-Regression method for
  objects of class 'runMetaAnalysis'
- [`metaRegression(`*`<meta>`*`)`](metaRegression.meta.md) :
  Meta-Regression method for objects of class 'runMetaAnalysis'
- [`metaRegression(`*`<rma>`*`)`](metaRegression.rma.md) :
  Meta-Regression method for objects of class 'runMetaAnalysis'
- [`metaRegression(`*`<runMetaAnalysis>`*`)`](metaRegression.runMetaAnalysis.md)
  : Meta-Regression method for objects of class 'runMetaAnalysis'
- [`exploreStudies()`](exploreStudies.md) : Explore included treatments
  and comparisons
- [`createStudyTable()`](createStudyTable.md) : Create study table
- [`createRobSummary()`](createRobSummary.md) : Create a summary risk of
  bias plot

## Effect Size Calculators

Default plug-in functions used by `calculateEffectSizes`.

- [`g.m.sd()`](g.m.sd.md) : Calculate Hedges' g using means and standard
  deviations
- [`g.change.m.sd()`](g.change.m.sd.md) : Calculate Hedges' g using
  within-group change data
- [`g.binary()`](g.binary.md) : Calculate Hedges' g using binary outcome
  data
- [`g.precalc()`](g.precalc.md) : Forward pre-calculated values of
  Hedges' g
- [`rr.precalc()`](rr.precalc.md) : Forward pre-calculated log-risk
  ratios
- [`rr.binary()`](rr.binary.md) : Calculate the log-risk ratio using
  binary outcome data

## S3 Methods & Helper Functions

Additional functionality for core functions.

- [`blup(`*`<runMetaAnalysis>`*`)`](blup.runMetaAnalysis.md) : Best
  Linear Unbiased Predictions (BLUPs) for 'runMetaAnalysis' models.
- [`blup()`](blup.md) : blup: Empirical Bayes estimates
- [`eb(`*`<runMetaAnalysis>`*`)`](eb.runMetaAnalysis.md) : Best Linear
  Unbiased Predictions (BLUPs) for 'runMetaAnalysis' models.
- [`eb()`](eb.md) : eb: Empirical Bayes estimates
- [`imputeResponse()`](imputeResponse.md) : Impute response rates based
  on continuous outcome data
- [`metapsyFindOutliers()`](metapsyFindOutliers.md) : Find Statistical
  Outliers in a Meta-Analysis
- [`metapsyInfluenceAnalysis()`](metapsyInfluenceAnalysis.md) :
  Influence Diagnostics
- [`profile(`*`<runMetaAnalysis>`*`)`](profile.runMetaAnalysis.md) :
  Profile Likelihood Plots for 'runMetaAnalysis' models.
- [`proportionMID()`](proportionMID.md) : Calculate the proportion of
  true effect sizes above a meaningful threshold
- [`simulateTreatmentCycles()`](simulateTreatmentCycles.md) : Simulate
  the number of treatment cycles and "excess treatments"
- [`addTrialArmInfo()`](addTrialArmInfo.md) : Add information that
  varies between trial arms as extra columns to your meta-analysis
  dataset
- [`plot(`*`<flagEffectSizes>`*`)`](plot.flagEffectSizes.md) : Plot
  method for objects of class 'flagEffectSizes'
- [`plot(`*`<proportionMID>`*`)`](plot.proportionMID.md) : Plot method
  for objects of class 'proportionMID'
- [`plot(`*`<runMetaAnalysis>`*`)`](plot.runMetaAnalysis.md) : Plot
  method for objects of class 'runMetaAnalysis'
- [`plot(`*`<simulateTreatmentCycles>`*`)`](plot.simulateTreatmentCycles.md)
  : Plot method for objects of class 'simulateTreatmentCycles'
- [`plot(`*`<subgroupAnalysis>`*`)`](plot.subgroupAnalysis.md) : Plot
  method for objects of class 'runMetaAnalysis'
- [`print(`*`<checkConflicts>`*`)`](print.checkConflicts.md) : Print
  method for the 'checkConflicts' function
- [`print(`*`<exploreStudies>`*`)`](print.exploreStudies.md) : Print
  method for objects of class 'exploreStudies'
- [`print(`*`<flagEffectSizes>`*`)`](print.flagEffectSizes.md) : Print
  method for objects of class 'flagEffectSizes'
- [`print(`*`<proportionMID>`*`)`](print.proportionMID.md) : Print
  method for objects of class 'proportionMID'
- [`print(`*`<runMetaAnalysis>`*`)`](print.runMetaAnalysis.md) : Print
  method for objects of class 'runMetaAnalysis'
- [`print(`*`<simulateTreatmentCycles>`*`)`](print.simulateTreatmentCycles.md)
  : Print method for objects of class 'simulateTreatmentCycles'
- [`print(`*`<subgroupAnalysis>`*`)`](print.subgroupAnalysis.md) : Print
  method for objects of class 'subgroupAnalysis'
- [`summary(`*`<runMetaAnalysis>`*`)`](summary.runMetaAnalysis.md) :
  Show details of 'runMetaAnalysis' class objects
- [`` `data<-`() ``](Replacement-functions.md)
  [`` `which.run<-`() ``](Replacement-functions.md)
  [`` `es.measure<-`() ``](Replacement-functions.md)
  [`` `es.type<-`() ``](Replacement-functions.md)
  [`` `es.var<-`() ``](Replacement-functions.md)
  [`` `se.var<-`() ``](Replacement-functions.md)
  [`` `es.binary.raw.vars<-`() ``](Replacement-functions.md)
  [`` `method.tau<-`() ``](Replacement-functions.md)
  [`` `i2.ci.threelevel<-`() ``](Replacement-functions.md)
  [`` `nsim.boot<-`() ``](Replacement-functions.md)
  [`` `hakn<-`() ``](Replacement-functions.md)
  [`` `study.var<-`() ``](Replacement-functions.md)
  [`` `arm.var.1<-`() ``](Replacement-functions.md)
  [`` `arm.var.2<-`() ``](Replacement-functions.md)
  [`` `measure.var<-`() ``](Replacement-functions.md)
  [`` `low.rob.filter<-`() ``](Replacement-functions.md)
  [`` `method.tau.ci<-`() ``](Replacement-functions.md)
  [`` `which.combine<-`() ``](Replacement-functions.md)
  [`` `which.combine.var<-`() ``](Replacement-functions.md)
  [`` `which.outliers<-`() ``](Replacement-functions.md)
  [`` `which.influence<-`() ``](Replacement-functions.md)
  [`` `which.rob<-`() ``](Replacement-functions.md)
  [`` `nntCer<-`() ``](Replacement-functions.md)
  [`` `nnt<-`( ``*`<runMetaAnalysis>`*`)`](Replacement-functions.md)
  [`` `rho.within.study<-`() ``](Replacement-functions.md)
  [`` `phi.within.study<-`() ``](Replacement-functions.md)
  [`` `power.within.study<-`() ``](Replacement-functions.md)
  [`` `w1.var<-`() ``](Replacement-functions.md)
  [`` `w2.var<-`() ``](Replacement-functions.md)
  [`` `time.var<-`() ``](Replacement-functions.md)
  [`` `vcov<-`() ``](Replacement-functions.md)
  [`` `near.pd<-`() ``](Replacement-functions.md)
  [`` `use.rve<-`() ``](Replacement-functions.md)
  [`` `html<-`() ``](Replacement-functions.md)
  [`` `lower.is.better<-`() ``](Replacement-functions.md)
  [`` `selmodelSteps<-`() ``](Replacement-functions.md)
  [`` `selmodel<-`( ``*`<runMetaAnalysis>`*`)`](Replacement-functions.md)
  [`rerun()`](Replacement-functions.md) : Replacement functions for
  "runMetaAnalysis" results objects

## Datasets

- [`depressionPsyCtr`](depressionPsyCtr.md) : The 'depressionPsyCtr'
  dataset
