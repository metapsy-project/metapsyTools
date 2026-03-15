# Replacement functions for "runMetaAnalysis" results objects

Once a model has been fitted using `runMetaAnalysis`, replacement
functions are defined for each function argument. This allows to quickly
tweak one or more analysis settings, which are implemented once the
`rerun` function is called.

## Usage

``` r
data(x) <- value
which.run(x) <- value
es.measure(x) <- value
es.type(x) <- value
es.var(x) <- value
se.var(x) <- value
es.binary.raw.vars(x) <- value
method.tau(x) <- value
i2.ci.threelevel(x) <- value
nsim.boot(x) <- value
hakn(x) <- value
study.var(x) <- value
arm.var.1(x) <- value
arm.var.2(x) <- value
measure.var(x) <- value
low.rob.filter(x) <- value
method.tau.ci(x) <- value
which.combine(x) <- value
which.combine.var(x) <- value
which.outliers(x) <- value
which.influence(x) <- value
which.rob(x) <- value
nntCer(x) <- value
# S3 method for class 'runMetaAnalysis'
nnt(x) <- value
rho.within.study(x) <- value
phi.within.study(x) <- value
power.within.study(x) <- value
w1.var(x) <- value
w2.var(x) <- value
time.var(x) <- value
vcov(x) <- value
near.pd(x) <- value
use.rve(x) <- value
html(x) <- value
lower.is.better(x) <- value
selmodelSteps(x) <- value
# S3 method for class 'runMetaAnalysis'
selmodel(x) <- value
rerun(m)
```

## Arguments

- x:

  Object of class `runMetaAnalysis`.

- value:

  Value of one of the arguments in `runMetaAnalysis` or
  `correctPublicationBias`

- m:

  (Adapted) object of class `runMetaAnalysis`.

## Details

`nntCer` / `nntCer<-` and `selmodelSteps` / `selmodelSteps<-` are also
available as S3 methods `nnt`, `nnt<-`, `selmodel`, `selmodel<-` for
`runMetaAnalysis` objects. The function `time.var` (and `time.var<-`) is
a regular replacement function, not an S3 method for `time`.

## See also

[`runMetaAnalysis`](runMetaAnalysis.md),
[`correctPublicationBias`](correctPublicationBias.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data("depressionPsyCtr")
depressionPsyCtr %>%
  checkDataFormat() %>%
  checkConflicts() %>%
  calculateEffectSizes() %>% 
  filterPoolingData(condition_arm2 %in% 
                      c("wl", "other ctr")) -> data

m <- runMetaAnalysis(data, "combined")

# Compare results when other tau^2 estimator is used
method.tau(m) <- "DL"
rerun(m)
} # }
```
