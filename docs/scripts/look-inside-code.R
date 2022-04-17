# Load package
library(metapsyTools)
library(meta)
data("inpatients")

# checkDataFormat -------------------------------------------------------------
res.cdf <- inpatients %>% checkDataFormat()

# checkConflicts --------------------------------------------------------------
res.cc <- res.cdf %>% checkConflicts()

# expandMultiarmTrials --------------------------------------------------------
res.emt <- res.cc %>% expandMultiarmTrials()
nrow(res.cc)
nrow(res.emt)

# calculateEffectSizes --------------------------------------------------------
res.ces <- res.emt %>% calculateEffectSizes()
res.ces[,c("study", "es", "se")]

res.ces <- res.emt %>% calculateEffectSizes(include.switched.arms = TRUE)
res.ces[,c("study", "es", "se")]

# runMetaAnalysis -------------------------------------------------------------xw
res.fpd <- res.ces %>% filterPoolingData(Cond_spec_trt1 == "cbt")

res.rma <- res.fpd %>% runMetaAnalysis()
funnel.meta(res.rma$model.combined)
summary(res.rma)