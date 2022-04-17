# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                             #
#   A Look Inside {metapsyTools}                                              #
#   Coding Example                                                            #
#   Last Update: May 2022                                                     #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# 0. Dependencies -------------------------------------------------------------

# Typcially, it is convenient to also load {meta} and {metafor} from the start,
# especially if we want to apply other meta-analytic techniques later on.
library(metapsyTools)
library(meta)
library(metafor)

# We also use functions included in {dplyr}, {dmetar} and {metasens}.
library(dplyr)
library(dmetar)
library(metasens)


# 1. Dataset ------------------------------------------------------------------

# In this example, we use the "psyCtrSubset" dataset
# This dataset is pre-installed in {metapsyTools}

# Eyeball the data
# Our data is in the "long" format:
data("psyCtrSubset")
glimpse(psyCtrSubset)

# Show the number of studies
psyCtrSubset %>% 
    pull(study) %>% 
    unique() %>% length()



# 2. Preparation --------------------------------------------------------------

# Data format check
psyCtrSubset %>% 
    checkDataFormat() -> psyCtrSubset

# Let's say we want the rob columns to have the class "numeric"
psyCtrSubset %>% 
    checkDataFormat(variable.class = list("rob" = "numeric")) -> psyCtrSubset


# Check data for conflicts
# Since checkDataFormat() has already been applied, the function knows that
# the data format is long
psyCtrSubset %>% 
    checkConflicts() -> psyCtrSubset

# Change vars.for.id so that non-unique IDs are generated
# This throws an error:
psyCtrSubset %>% 
    checkConflicts(vars.for.id = c("study", "Time")) -> psyCtrSubsetError


# Now, we expand the multiarm trials
# First we check the number of rows *before* the expansion
nrow(psyCtrSubset)
psyCtrSubset %>% 
    expandMultiarmTrials() -> psyCtrSubset
nrow(psyCtrSubset)

# For the object containing ID conflicts, this will not work:
psyCtrSubsetError %>% 
    expandMultiarmTrials()

# Now, we calculate the effect sizes ...
# ... using default settings
psyCtrSubset %>% 
    calculateEffectSizes() %>% 
    select(study, es, se)

# ... using changed signs (positive effects)
psyCtrSubset$changeSign <- rep(TRUE, nrow(psyCtrSubset))
psyCtrSubset %>% 
    calculateEffectSizes(change.sign = "changeSign") %>% 
    select(study, es, se)

# ... including switched arms
psyCtrSubset %>% 
    calculateEffectSizes(include.switched.arms = TRUE) %>% 
    select(study, es, se)


# Typically, one would not apply these function step by step, 
# but include them in a single pipe:
data("psyCtrSubset")
psyCtrSubset %>% 
    checkDataFormat(variable.class = list("rob" = "numeric"),
                    data.format = "long") %>% 
    checkConflicts() %>% 
    expandMultiarmTrials() %>% 
    calculateEffectSizes(include.switched = TRUE) -> data.prep

write.csv(data.prep, "data_prep.csv")



# 3. Analysis -----------------------------------------------------------------


## 3.1 Filtering Data ---------------------------------------------------------

# Filter data for our planned analysis
instrumentPriority = c("hdrs", "bdi-2", "ces-d", "phq-9")

data.prep %>% 
    filterPoolingData(Detect(Cond_spec_trt1, "cbt|pst"),
                      Detect(Cond_spec_trt2, "cau|wl")) -> data

# For illustration: filter using a priority rule;
# In this case, only studies that employed one of the instruments in 
# 'instrumentPriority' will be included, and the rest is discarded!
instrumentPriority = c("hdrs", "bdi-2", "ces-d", "phq-9") 

data.prep %>% 
    filterPoolingData(Detect(Cond_spec_trt1, "cbt|pst"),
                      Detect(Cond_spec_trt2, "cau|wl")) %>% 
    filterPriorityRule(Outc_measure = instrumentPriority) %>% 
    pull(Outc_measure)



## 3.2 Effect Size Pooling ----------------------------------------------------

# Now we pool the effects using runMetaAnalysis.
# We change a few of the default settings to showcase the functiontality.
res <- runMetaAnalysis(data,
                       method.tau = "PM",
                       extra.grouping.var = "Cond_spec_trt1",
                       low.rob.filter = "rob >= 3",
                       which.outliers = "combined",
                       which.influence = "combined",
                       which.rob = "combined",
                       nnt.cer = 0.17,
                       rho.within.study = 0.65)
res

# Print technical details and summary forest plot
summary(res)



## 3.3 Forest Plots -----------------------------------------------------------

# Forest plot of the 'combined' analysis - check if additional grouping worked
# Width and height are specified in *inches* (sigh!)
pdf("forest_plot.pdf", width=15, height=20)
plot(res, "combined")
dev.off()

# More styling features are available if we extract the fitted {meta} object
# directly and plug it into 'forest'.
# Run names(pdfFonts()) for available fonts.
pdf("forest_plot_styled.pdf", width=15, height=20)
forest.meta(res$model.combined,
            col.square = "dodgerblue3",
            fontfamily = "Palatino",
            sortvar = TE,
            print.tau2 = FALSE,
            col.predict = "black",
            smlab = "Effect on Depressive \nSymptom Severity",
            leftlabs = c("Study", "g", "S.E.")) 
dev.off()



## 3.4 Model Diagnostics ------------------------------------------------------

# Plot a baujat plot
plot(res, "baujat")

# Plot results of the 'Leave-One-Out' influence analysis, sorted by I2
plot(res, "loo-i2")

# For the multilevel-level (CHE) models, a the profile log-likelihoods
# of the variance components (tau2 on both levels) should be checked
profile(res$model.threelevel.che)



## 3.5 Small Study-Effects ----------------------------------------------------

# Funnel plot
col.contour <- c("gray75", "gray85", "gray95")
funnel.meta(res$model.overall, 
            contour = c(0.9, 0.95, 0.99),
            col.contour = col.contour)
legend(x = -8, y = 0.01, 
       legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
       fill = col.contour)


# Eggers' test using the Pustejovsky-Rodgers adjustment:
res$model.combined$n.e <- res$model.combined$data$Post_N_trt1
res$model.combined$n.c <- res$model.combined$data$Post_N_trt2
metabias(res$model.combined,
         method = "Pustejovsky")


# Duval & tweedie trim-and-fill method
tf <- trimfill(res$model.combined) 
funnel.meta(tf, contour = c(0.9, 0.95, 0.99),
            col.contour = col.contour)


# Pcurve method: this method should only be applied for analyses
# with low heterogeneity; therefore we use the outliers removed analysis
sampleSize <- with(res$model.outliers$data,{Post_N_trt1 + Post_N_trt2})

pc <- pcurve(res$model.outliers,
             effect = TRUE,
             dmin = -1, dmax = 1,
             N = sampleSize)

# Rücker’s Limit Meta-Analysis Method
lim <- limitmeta(res$model.combined)
funnel.limitmeta(lim, shrunken = TRUE)



## 3.6 Meta-Regression --------------------------------------------------------

# Meta-Regression using the 'metaRegression' function
# Here we conduct a multiple meta-regression using year and risk of bias
# as predictor. We center and scale both variables in the model.

# ... using the 'combined' model:
metaRegression(res$model.combined, ~ scale(year) + scale(rob))

# ... using the multilevel CHE model:
metaRegression(res$model.threelevel.che, ~ scale(year) + scale(rob))

# PET model (small study effects-corrected effect size)
metaRegression(res$model.combined, ~ se)



## 3.7 Subgroup Analysis ------------------------------------------------------

# We want to create a table with subgroup analysis results for 
# 'type_format_trt2', 'country', 'age_group'.
# We use the 'combined' analysis for subgroup analyses
sg <- subgroupAnalysis(res, format_trt2, country, age_group,
                       .which.run = "combined")

# Generate forest plot for "country"
plot(sg, "country")



## 3.8 Study Table ------------------------------------------------------------

# Create a study table for the meta-analysis
createStudyTable(res, study, 
                 Cond_spec_trt1, Cond_spec_trt2,
                 Post_N_trt1, Post_N_trt2,
                 n_sessions_trt1, 
                 mean_age_trt1, percent_women_trt1,
                 country = c("Canada" = "4", "Europe" = "3", "USA" = "1",
                             "Asia" = "6", "Middle East" = "7", 
                             "Australia" = "5"), 
                 ac, ba, itt, sg,
                 .round.by.digits = list(mean_age_trt1 = 1),
                 .na.replace = "n.r.",
                 .column.names = list(Cond_spec_trt1 = "IG",
                                      Cond_spec_trt2 = "CG",
                                      Post_N_trt1 = "N.ig",
                                      Post_N_trt2 = "N.cg",
                                      n_sessions_trt1 = "n.sessions",
                                      mean_age_trt1 = "age",
                                      percent_women_trt1 = "perc. women")
)

