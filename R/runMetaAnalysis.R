#' Run different types of meta-analyses
#'
#' @description  This wrapper function allows to simultaneously pool effect sizes using
#' different meta-analytic approaches.
#' \loadmathjax
#'
#' @usage runMetaAnalysis(data,
#' 
#'                 # Models to run
#'                 which.run = c("overall", "combined",
#'                               "lowest.highest", "outliers",
#'                               "influence", "rob", "threelevel",
#'                               "threelevel.che"),
#'                               
#'                 # Effect size measure
#'                 es.measure = c("g", "RR"),
#'                 es.type = c("precalculated", "raw"),
#'                 es.var = ifelse(identical(es.measure[1], "g"), 
#'                                 ".g", ".log_rr"),
#'                 se.var = ifelse(identical(es.measure[1], "g"), 
#'                                 ".g_se", ".log_rr_se"),
#'                 es.binary.raw.vars = 
#'                   c(".event_arm1", ".event_arm2",
#'                     ".totaln_arm1", ".totaln_arm2"),
#'                     
#'                 # Estimator of the heterogeneity variance
#'                 method.tau = "REML",
#'                 method.tau.ci = "Q-Profile",
#'                 i2.ci.boot = FALSE,
#'                 nsim.boot = 5e3,
#'                 hakn = TRUE,
#'                 
#'                 # Data specifications
#'                 study.var = "study",
#'                 arm.var.1 = "condition_arm1",
#'                 arm.var.2 = "condition_arm2",
#'                 measure.var = "instrument",
#'                 low.rob.filter = "rob > 2",
#'                 round.digits = 2,
#'                 
#'                 # Model specifications
#'                 which.combine = c("arms", "studies"),
#'                 which.combine.var = "multi_arm1",
#'                 which.outliers = c("overall", "combined"),
#'                 which.influence = c("overall", "combined"),
#'                 which.rob = c("overall", "combined"),
#'                 nnt.cer = 0.2,
#'                 rho.within.study = 0.6,
#'                 phi.within.study = 0.9,
#'                 w1.var = ifelse(identical(es.measure[1], "g"), 
#'                                 "n_arm1", "totaln_arm1"),
#'                 w2.var = ifelse(identical(es.measure[1], "g"), 
#'                                 "n_arm2", "totaln_arm2"),
#'                 time.var = "time_weeks",
#'                 vcov = c("simple", "complex"),
#'                 near.pd = FALSE,
#'                 use.rve = TRUE,
#'                 
#'                 # Output
#'                 html = TRUE,
#'                 
#'                 # Additional arguments
#'                 ...)
#'
#' @param data \code{data.frame}. Effect size data in the wide format, as 
#' created by \code{\link{calculateEffectSizes}}.
#' @param which.run \code{character}. Selection of models to be calculated. 
#' See 'Details'.
#' @param es.measure `character`. Should meta-analyses be calculated using the
#' bias-corrected standardized mean difference (`"g"`; default), or using 
#' risk ratios (`"RR"`)? Meta-analyses will only be conducted using comparisons
#' that contain non `NA` values in the `es.var` and `se.var` columns.
#' @param es.type `character`. Should pre-calculated or raw event data (i.e.
#' the Mantel-Haenszel method) be used for meta-analyses of risk ratios?
#' Can be set to `"precalculated"` (default) or `"raw"`.
#' @param es.var `character`. Specifies the name of the variable containing the (pre-
#' calculated) effect size data in `data`. When `es.measure = "g"`, this is set to
#' `.g` by default; `".log_rr"` is used when `es.measure = "RR"`. The default settings
#' correspond with the standard output of [calculateEffectSizes()].
#' @param se.var `character`. Specifies the name of the variable containing the (pre-calculated)
#' standard errors (square root of the variance) of the effect size metric defined
#' in `es.var`. If `es.measure = "g"`, this is automatically set to `.g_se`; if 
#' `es.measure = "RR"`, `".log_rr_se"` is used. The default settings
#' correspond with the standard output of [calculateEffectSizes()].
#' @param es.binary.raw.vars `character` vector, defining the column names in
#' which the (1) raw event counts in the experimental group, (2) raw event counts in the
#' control/reference group, (3) sample size in the experimental group, and (4) sample size
#' of the control/reference group are stored. 
#' Defaults correspond with the standard output of [calculateEffectSizes()].
#' @param method.tau \code{character}. A character string indicating 
#' which method is used to estimate the between-study variance (tau-squared) 
#' and its square root (tau). Either \code{"REML"} (default), \code{"DL"},
#' \code{"PM"}, \code{"ML"}, \code{"HS"}, \code{"SJ"}, \code{"HE"}, or \code{"EB"}, can be abbreviated (see
#' \code{\link[meta]{metagen}}). Use \code{"FE"} to use a fixed-effect/"common effect" model.
#' @param method.tau.ci \code{character}. A character string indicating which method is used to estimate the
#' confidence interval of the between-study heterogeneity variance \mjeqn{\tau^2}{\tau^2}. Either \code{"Q-Profile"} (default and recommended; [Viechtbauer, 2017](http://www.wvbauer.com/lib/exe/fetch.php/articles:viechtbauer2007b.pdf)),
#' \code{"BJ"}, \code{"J"}, or \code{"PL"} can be abbreviated. See \code{\link[meta]{metagen}} and \code{\link[metafor]{rma.uni}} for details.
#' @param i2.ci.boot `logical`. Confidence intervals for \mjeqn{\tau^2}{\tau^2} as calculated by the Q-Profile method are not 
#' directly applicable for three-level models, in which two heterogeneity variance components are estimated. By default,
#' this argument is therefore set to `FALSE`, and no confidence intervals around \mjeqn{\tau^2}{\tau^2} and \mjeqn{I^2}{I^2} are provided for the
#' `"threelevel"` and `"threelevel.che"` model.
#' If this argument is set to `TRUE`, parametric bootstrapping is used to calculate confidence intervals around the between- and
#' within-study heterogeneity estimates (\mjeqn{\tau}{\tau} and \mjeqn{I^2}{I^2}). Please note that this can take several minutes,
#' depending on the number of effect sizes. If [correctPublicationBias()] is used and `i2.ci.boot` is `TRUE`,
#' bootstrapping will also be used to calculate confidence intervals around the \mjeqn{G^2}{G^2} statistic 
#' ([Rücker et al., 2011](https://academic.oup.com/biostatistics/article/12/1/122/391113))
#' used in the limit meta-analysis (note that \mjeqn{G^2}{G^2} is printed as \mjeqn{I^2}{I^2} in this package).
#' @param nsim.boot `numeric` Number of bootstrap samples to be drawn when `i2.ci.boot` is `TRUE`. Defaults to
#' 5000.
#' @param hakn \code{logical}. Should the Knapp-Hartung adjustment for effect size significance tests be used? Default is \code{TRUE}.
#' @param study.var \code{character}. The name of the variable in \code{data} in which the study IDs are stored.
#' @param arm.var.1 \code{character}. The name of the variable in \code{data} in which the condition (e.g. "guided iCBT")
#' of the \emph{first} arm within a comparison are stored.
#' @param arm.var.2 \code{character}. The name of the variable in \code{data} in which the condition (e.g. "wlc")
#' of the \emph{second} arm within a comparison are stored.
#' @param measure.var \code{character}. The name of the variable in \code{data} 
#' in which the instrument used for the comparison is stored.
#' @param low.rob.filter \code{character}. A filtering statement by which to 
#' include studies for the "low RoB only" analysis. Please note that the name of 
#' the variable must be included as a column in \code{data}.
#' @param round.digits \code{numeric}. Number of digits to round the (presented) results by. Default is \code{2}.
#' @param which.combine `character`. Should multiple effect sizes within one study be pooled
#' on an `"arms"` (default) or `"studies"` level? When a study is a multi-arm trial, setting
#' `which.combine = "arms"` will aggregate the effect sizes for each trial arm individually before pooling;
#' the `which.combine.var` argument can be used to control which effects within a study should be aggregated.
#' When `which.combine = "studies"`, one overall aggregated effect is created for each study. This setting is preferable
#' from a statistical perspective, since it ensures that all pooled effects can be assumed to be independent.
#' @param which.combine.var \code{character}. Additional grouping variable within studies to be used for the \code{"combined"}
#' analysis when `which.combine = "arms"`. If the specified variable differs within one study (as defined
#' by `study.var`), effects will be aggregated separately for each unique value in `which.combine.var`.
#' Defaults to `"multi_arm1"`, the variable which encodes multi-arm intervention conditions in the
#' [Metapsy data standard](https://docs.metapsy.org/data-preparation/format/#standard-variables).
#' @param which.outliers \code{character}. Which model should be used to conduct outlier analyses? Must be
#' \code{"overall"} or \code{"combined"}, with \code{"overall"} being the default.
#' @param which.influence \code{character}. Which model should be used to conduct influence analyses? Must be
#' \code{"overall"} or \code{"combined"}, with \code{"overall"} being the default.
#' @param which.rob \code{character}. Which model should be used to conduct the "low risk of bias only" analyses? Must be
#' \code{"overall"} or \code{"combined"}, with \code{"overall"} being the default.
#' @param nnt.cer \code{numeric}. Value between 0 and 1, indicating the assumed control group event rate to be used
#' for calculating NNTs via the Furukawa-Leucht method.
#' @param rho.within.study \code{numeric}. Value between 0 and 1, indicating the assumed correlation of effect sizes
#' within studies. This is relevant to combine effect sizes for the \code{"combined"} analysis type, and used to estimate
#' the variance-covariance matrices needed for the conditional and hierarchical effects three-level model. Default is \code{0.6}.
#' @param phi.within.study \code{numeric}. Value between 0 and 1, indicating the assumed one-week autocorrelation of effect sizes. 
#' This is only used when `vcov="complex"` to approximate the variance-covariance matrices needed for the `"combined"` and 
#' `"threelevel.che"` model. Default is 0.9. See "Details".
#' @param w1.var `character`. Name of the variable in `data` in which the sample 
#' sizes of the first arm are stored. See "Details".
#' @param w2.var `character`. Name of the variable in `data` in which 
#' the sample sizes of the second arm are stored. See "Details".
#' @param time.var `character`. Name of the variable in `data` in which the assessment time point is stored.
#' Should be expressed as weeks since randomization; other units (e.g. days, months) are also possible if `phi.within.study`
#' is specified accordingly. See "Details".
#' @param vcov `character`. For the `"combined"` and `"threelevel.che"` model, should variance-covariance matrices
#' (representing the dependency structure of the data) be approximated using a heterogeneous compound symmetry (`"simple"`; default)
#' or unstructured matrix structure (`"complex"`)? See "Details".
#' @param near.pd `logical`. If at least one of the study variance-covariance matrices constructed
#' when `vcov="complex"` is not positive definite/invertible, should the \code{\link[Matrix]{nearPD}} function
#' be used to compute the nearest positive definite matrix? Default is \code{FALSE}.
#' @param use.rve \code{logical}. Should robust variance estimation be used to calculate confidence intervals and tests of three-level models?
#' \code{TRUE} by default.
#' @param html \code{logical}. Should an HTML table be created for the results? Default is \code{TRUE}.
#' @param ... Additional arguments.
#'
#'
#' @return Returns an object of class \code{"runMetaAnalysis"}. This object includes, among other things,
#' a \code{data.frame} with the name \code{summary}, in which all results are summarized - including the
#' studies which were removed for some analysis steps. Other objects are the "raw" model objects returned
#' by all selected analysis types. This allows to conduct further operations on some models specifically (e.g.
#' run a meta-regression by plugging one of the model objects in \code{update.meta}.
#'
#'
#' @examples
#' \dontrun{
#' data("depressionPsyCtr")
#' library(meta)
#' 
#' depressionPsyCtr %>%
#'   checkDataFormat() %>%
#'   checkConflicts() %>%
#'   calculateEffectSizes() %>% 
#'   filterPoolingData(condition_arm2 %in% 
#'                       c("wl", "other ctr")) -> data
#' 
#' # Run the meta-analyses
#' runMetaAnalysis(data) -> res
#' 
#' # Use replacement function to show results for
#' # differing settings
#' method.tau(res) <- "PM"
#' hakn(res) <- FALSE
#' rerun(res)
#' 
#' # Show summary
#' res
#' 
#' # Show forest plot (by default, "overall" is used)
#' plot(res)
#' 
#' # Show forest plot of specific analysis
#' plot(res, "outliers")
#' plot(res, "threelevel")
#' plot(res, "baujat")
#' plot(res, "influence")
#' plot(res, "lowest.highest")
#' 
#' # Extract specific model and do further calculations
#' # (e.g. meta-regression on 'year')
#' metaRegression(res$model.overall, ~ scale(year))
#' 
#' # Conduct a subgroup analysis
#' subgroupAnalysis(res, country)
#' 
#' # Correct for publication bias/small-study effects
#' correctPublicationBias(res)
#' 
#' # For the combined analysis, set which.combine to
#' # "studies" here, so that all effects in a study are aggregated
#' # first before pooling
#' data %>% 
#'   runMetaAnalysis(which.combine = "studies") %>% 
#'   plot("combined")
#' 
#' # Run meta-analysis using raw response rate data
#' data %>% 
#'   runMetaAnalysis(es.measure = "RR",
#'                   es.type = "raw")
#' }
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{calculateEffectSizes}}, 
#' \code{\link{subgroupAnalysis}}. \code{\link{correctPublicationBias}},
#' \code{\link{metaRegression}}
#'
#' @details The \code{runMetaAnalysis} function is a wrapper for several types of meta-analytic models
#' that are typically used. It allows to run all of these models in one step in order to generate results
#' that are somewhat closer to being "publication-ready".
#'
#' By default, the following models are calculated:
#'
#' \itemize{
#'   \item \code{"overall"}. Runs a generic inverse-variance (random-effects) model. All included
#'   effect sizes are treated as independent. When `es.measure = "RR"` and `es.type = "raw"`,
#'   the Mantel-Haenszel method is used for pooling instead.
#'   \item \code{"combined"}. Pools all effect sizes within one study (defined by \code{study.var}) before
#'   pooling. This ensures that all effect sizes are independent (i.e., unit-of-analysis error &
#'   double-counting is avoided). To combine the effects, one has to assume a correlation of effect sizes
#'   within studies, empirical estimates of which are typically not available.
#'   \item \code{"lowest.highest"}. Runs a meta-analysis, but with only (i) the lowest and (ii) highest
#'   effect size within each study included.
#'   \item \code{"outlier"}. Runs a meta-analysis without statistical outliers (i.e. effect sizes for which
#'   the confidence interval does not overlap with the confidence intervall of the overall effect).
#'   \item \code{"influence"}. Runs a meta-analysis without influential cases (see \code{\link[metafor]{influence.rma.uni}} for
#'   details).
#'   \item \code{"rob"}. Runs a meta-analysis with only low-RoB studies included.
#'   \item \code{"threelevel"}. Runs a multilevel (three-level) meta-analysis model, with effect sizes nested
#'   in studies.
#'   \item \code{"threelevel.che"}. Runs a multilevel (three-level) meta-analysis model, with effect sizes nested
#'   in studies. Variance-covariance matrices of each study with two or more effect sizes are estimated using
#'   \code{rho.within.study} as the assumed overall within-study correlation. This imputation allows to run a "correlated and
#'   hierarchical effects" (CHE) model, which is typically a good approximation for data sets with unknown and/or
#'   complex dependence structures.
#' }
#' 
#' Internally, the `overall`, `combined`, `lowest.highest`, `outlier`, `influence` and `rob`
#' models are fitted by calling the [meta::metagen()] or [meta::metabin()] function,
#' respectively, in **\{meta\}** (Balduzzi, Rücker & Schwarzer, [2019](https://pubmed.ncbi.nlm.nih.gov/31563865/)). 
#' The `threelevel` and `threelevel.che` models are implemented using [metafor::rma.mv()]
#' in **\{metafor\}** (Viechtbauer, [2005](https://www.jstatsoft.org/article/view/v036i03)). 
#' 
#' Outlier selection is implemented using the [dmetar::find.outliers()] function,
#' and influence analyses using the [dmetar::InfluenceAnalysis()] function. The latter function is a wrapper for 
#' [metafor::influence.rma.uni()].
#' 
#' \mjeqn{~}{~}
#'  
#' ## Simple or complex variance-covariance approximation
#' 
#' The `vcov` argument controls if the effect size dependencies within the data
#' should be approximated using a `"simple"` (default) or more `"complex"` (but potentially more accurate)
#' method. This argument is only relevant for the `"combined"` and `"threelevel.che"` models. The default 
#' "simple" method constructs variance-covariance matrices \mjeqn{\Sigma_k}{\Sigma_k} for each study using a 
#' constant sampling correlation \mjeqn{\rho}{\rho} (defined by `rho.within.study`), which is identical across all studies, outcomes, and time points.
#' This simplifying assumption is part of the formulation of the CHE model originally provided by 
#' Pustejovsky and Tipton ([2022](https://link.springer.com/article/10.1007/s11121-021-01246-3)). 
#' 
#' Naturally, employing a common value of \mjeqn{\rho}{\rho} across all studies may not be reasonable
#' in some analyses, and other information may be available to better approximate the effect size dependencies
#' in the collected data. Setting `vcov` to `"complex"` allows to assume that correlations between effect 
#' sizes may differ conditional on the type of dependency. This means that the variance-covariance matrix \mjeqn{\Sigma_k}{\Sigma_k} of some study \mjeqn{k}{k} is
#' approximated by an unstructured matrix with varying \mjeqn{\rho_{ij}}{\rho_{ij}} (instead of a 
#' heterogeneous compound symmetry matrix with fixed \mjeqn{\rho}{\rho}, as is used when `vcov="simple"`).
#' 
#' 
#' \mjtdeqn{\small \begin{array}{ccc} \texttt{vcov="simple"} & \texttt{vcov="complex"} & \\\ \Sigma_k =\left[ \begin{array}{cccc} \sigma^2_1 & & & \\\ \rho \sigma_2 \sigma_1 & \sigma^2_2 & & \\\ \rho \sigma_3 \sigma_1 & \rho \sigma_3 \sigma_2 & \sigma^2_3 & \\\ \rho \sigma_4 \sigma_1 & \rho \sigma_4 \sigma_2 & \rho \sigma_4 \sigma_3 & \sigma^2_4 \end{array} \right] & \left[ \begin{array}{cccc} \sigma^2_1 & & & \\\ \rho_{21} \sigma_2 \sigma_1 & \sigma^2_2 & & \\\ \rho_{31} \sigma_3 \sigma_1 & \rho_{32} \sigma_3 \sigma_2 & \sigma_3 & \\\ \rho_{41} \sigma_4 \sigma_1 & \rho_{42} \sigma_4 \sigma_2 & \rho_{43} \sigma_4 \sigma_3 & \sigma^2_4 \end{array} \right] &  \end{array}}{\begin{array}{ccc}\texttt{vcov="simple"} & \texttt{vcov="complex"} & \\\\\ \Sigma_k = \begin{bmatrix} \sigma^2_1 \\\\\ \rho \sigma_2 \sigma_1 & \sigma^2_2 & & \\\\\ \rho \sigma_3 \sigma_1 & \rho \sigma_3 \sigma_2 & \sigma^2_3 & \\\\\ \rho \sigma_4 \sigma_1 & \rho \sigma_4 \sigma_2 & \rho \sigma_4 \sigma_3 & \sigma^2_4 \end{bmatrix} & \Sigma_k = \begin{bmatrix} \sigma^2_1 & & & \\\\\ \rho_{21} \sigma_2 \sigma_1 & \sigma^2_2 & & \\\\\ \rho_{31} \sigma_3 \sigma_1 & \rho_{32} \sigma_3 \sigma_2 & \sigma^2_3 & \\\\\ \rho_{41} \sigma_4 \sigma_1 & \rho_{42} \sigma_4 \sigma_2 & \rho_{43} \sigma_4 \sigma_3 & \sigma^2_4 \end{bmatrix} &  \end{array}}{}
#'  
#' For example, setting `vcov = "complex"` allows to additionally incorporate assumed correlations specific to multiple testing over time 
#' (e.g. correlations between effects at post-test and long-term follow-up). The value provided in
#' `phi.within.study` represents the (auto-)correlation coefficient \mjeqn{\phi}{\phi}, which serves
#' as a rough estimate of the re-test correlation after 1 week. When a vector of follow-up lengths 
#' is provided in `time.var`, this allows to model the gradual decrease in correlation between measurements 
#' over time. Furthermore, it is possible to calculate a correlation coefficient \mjeqn{\rho_w}{\rho_w} for 
#' multi-arm trials, which is directly proportional to the size of each individual trial arm. When all trial 
#' arms have the same size, meaning that each arm's weight \mjeqn{w}{w} is identical, 
#' \mjeqn{\rho_w}{\rho_w} is known to be 0.5. Multiarm weights \mjeqn{w}{w} (and thus 
#' \mjeqn{\rho_w}{w}) can be derived if the `w1.var` and `w2.var` variables,
#' containing the sample size of each study arm, are provided. 
#' 
#' Using the complex approximation method increases the risk that at least one studies' \mjeqn{\Sigma_k}{\Sigma_k} matrix
#' is not positive definite. In this case, the function automatically switches back to the 
#' constant sampling correlation approximation.
#' 
#' \mjeqn{~}{~}
#' 
#' ## Replacement functions
#' Once a model has been fitted using `runMetaAnalysis`, **replacement functions** 
#' are defined for each function argument. This allows to quickly tweak one or 
#' more analysis settings, which are implemented once the `rerun` function is called.
#' Say that we saved the results of `runMetaAnalysis` in an object `m`. If, for example, we want 
#' to check the results using a different estimator of \mjeqn{\tau^2}{\tau^2}, leaving
#' all other settings the same, we could run e.g. `method.tau(m) <- "PM"`, followed by
#' `rerun(m)`. This would provide results using the Paule-Mandel estimator. A list of all available
#' setting replacement functions is provided [here](https://tools.metapsy.org/reference/replacement-functions).
#' 
#' 
#' For more details see the [Get Started](https://tools.metapsy.org/articles/metapsytools) vignette.
#'
#' @import dplyr
#' @import mathjaxr
#' @importFrom crayon green yellow cyan bold
#' @importFrom scales pvalue
#' @importFrom purrr map
#' @importFrom meta update.meta metagen
#' @importFrom metafor escalc aggregate.escalc rma.mv vcalc blsplit simulate.rma
#' @importFrom clubSandwich coef_test conf_int
#' @importFrom stats dffits model.matrix rnorm rstudent complete.cases median quantile
#' @importFrom utils combn
#' @export runMetaAnalysis

runMetaAnalysis = function(data,
                           which.run = c("overall", "combined",
                                         "lowest.highest", "outliers",
                                         "influence", "rob", "threelevel",
                                         "threelevel.che"),
                           es.measure = c("g", "RR"),
                           es.type = c("precalculated", "raw"),
                           es.var = ifelse(identical(es.measure[1], "g"), 
                                           ".g", ".log_rr"),
                           se.var = ifelse(identical(es.measure[1], "g"), 
                                           ".g_se", ".log_rr_se"),
                           es.binary.raw.vars = 
                             c(".event_arm1", ".event_arm2",
                               ".totaln_arm1", ".totaln_arm2"),
                           method.tau = "REML",
                           method.tau.ci = "Q-Profile",
                           i2.ci.boot = FALSE,
                           nsim.boot = 5e3,
                           hakn = TRUE,
                           study.var = "study",
                           arm.var.1 = "condition_arm1",
                           arm.var.2 = "condition_arm2",
                           measure.var = "instrument",
                           low.rob.filter = "rob > 2",
                           round.digits = 2,
                           which.combine = c("arms", "studies"),
                           which.combine.var = "multi_arm1",
                           which.outliers = c("overall", "combined"),
                           which.influence = c("overall", "combined"),
                           which.rob = c("overall", "combined"),
                           nnt.cer = 0.2,
                           rho.within.study = 0.6,
                           phi.within.study = 0.9,
                           w1.var = ifelse(identical(es.measure[1], "g"), 
                                           "n_arm1", "totaln_arm1"),
                           w2.var = ifelse(identical(es.measure[1], "g"), 
                                           "n_arm2", "totaln_arm2"),
                           time.var = "time_weeks",
                           vcov = c("simple", "complex"),
                           near.pd = FALSE,
                           use.rve = TRUE,
                           html = TRUE,
                           ...){
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  0. Input Checks                                                          #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  # Check class
  if (!inherits(data, "data.frame") & 
      !inherits(data, "metapsyDatabase")){
    stop("'data' must be a data.frame or 'metapsyDatabase' R6 object.")
  }
  
  if (inherits(data, "metapsyDatabase")){
    data = data[["data"]]
  }
  
  # If es.measure is "g", es.type must be "precalculated"
  if (!(identical(es.measure[1], "g") || 
        identical(es.measure[1], "RR"))){
    stop("'es.measure' must be 'g' or 'RR'.")
  }
  if (identical(es.measure[1], "g")){
    es.type = "precalculated"
  }
  
  # Throw out all NAs/Missings
  if (identical(es.type[1], "raw")){
    data = data[rowSums(is.na(data[es.binary.raw.vars])) == 0 &
                  rowSums(data[es.binary.raw.vars] == Inf) == 0 &
                  rowSums(data[es.binary.raw.vars] == -Inf) == 0,]
    .type.es = ifelse(identical(es.measure[1], "g"), "g", "RR")
    .raw.bin.es = TRUE
  } else {
    data = data[!is.na(data[[es.var[1]]]) & 
                  data[[es.var[1]]] != Inf & 
                  data[[es.var[1]]] != -Inf,]
    .type.es = ifelse(identical(es.measure[1], "g"), "g", "RR")
    .raw.bin.es = FALSE
  }
  
  if (nrow(data) < 3){
    stop("Meta-analyses can only be run with at least 3 effect sizes.")
  }
  
  if (!identical(which.combine[1], "studies") &&
      !identical(which.combine[1], "arms")){
    stop("'which.combine' must be either 'arms' or 'studies'.")
  }
  
  if (!is.numeric(nsim.boot)){
    stop("'nsim.boot' must numeric.")
  }
  
  if (nsim.boot < 1000 && i2.ci.boot){
    warning("Number of boostrap samples (nsim.boot) is low. ",
            "Use a higher number for reliable inferences.")
  }
  
  # 'subset' is not allowed in data
  if ("subset" %in% colnames(data)){
    stop("Columns with the name 'subset' are not allowed in the data set. Did you run checkDataFormat()?")
  }
  
  # 'subset' is not allowed in data
  if ("exclude" %in% colnames(data)){
    stop("Columns with the name 'exclude' are not allowed in the data set. Did you run checkDataFormat()?")
  }
  
  # Get three-dots arguments;
  # Initialize error model list
  dots = list(...)
  error.model.list = list()
  warn.end = FALSE
  data.original = data
  method.tau = method.tau
  
  # Send message (beginning of analyses)
  sendMessage("start", .type.es = .type.es, 
              es.type = es.type)
  
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  1. Overall Model                                                         #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  # Fit model
  mGeneral =
    fitOverallModel(data, es.var, se.var, arm.var.1, arm.var.2,
                    measure.var, study.var, .raw.bin.es, .type.es, hakn,
                    method.tau.meta, method.tau.ci, method.tau,
                    dots, es.binary.raw.vars, round.digits,
                    nnt.cer, which.run)
  rownames(mGeneral$res) = "Overall"
  sendMessage(mGeneral, "overall", which.run)
  
  # If model failed, add to error model list
  if (mGeneral$has.error){
    error.model.list = append(error.model.list, "overall")
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  2. Lowest Only                                                           #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  if ("lowest.highest" %in% which.run){
    
    mLowest = fitLowestModel(data, study.var, multi.study,
                             mGeneral, .type.es, round.digits,
                             .raw.bin.es, nnt.cer)
    rownames(mLowest$res) = "One ES/study (lowest)"
    
    # If model failed, add to error model list
    if (mLowest$has.error){
      error.model.list = append(error.model.list, "lowest")
      if ("overall" %in% error.model.list){
        mLowest$message = mGeneral$message
      }
    }
    sendMessage(mLowest, "lowest.highest", which.run)
  } else {
    mLowest = NULL
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  3. Highest Only                                                          #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  if ("lowest.highest" %in% which.run){
    mHighest = 
      fitHighestModel(
        data, study.var, multi.study,
        mGeneral, .type.es, round.digits,
        .raw.bin.es, nnt.cer)
    rownames(mHighest$res) = "One ES/study (highest)"
    # If model failed, add to error model list
    if (mHighest$has.error){
      error.model.list = append(error.model.list, "highest")
      if ("overall" %in% error.model.list){
        mHighest$message = mGeneral$message
      }
    }
    sendMessage(mHighest, "lowest.highest", which.run)
  } else {
    mHighest = NULL
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  4. Combined Effect Sizes                                                 #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  if ("combined" %in% which.run){
    mCombHasError = FALSE
    if (identical(vcov[1], "complex")){
      mComb = fitCombinedHACEModel(which.combine, which.combine.var, measure.var,
                                   data, study.var, multi.study, es.var, se.var,
                                   mGeneral, .type.es, round.digits, hakn,
                                   .raw.bin.es, nnt.cer, rho.within.study,
                                   method.tau, method.tau.ci, dots,
                                   es.binary.raw.vars, arm.var.1, arm.var.2,
                                   phi.within.study, n.var.arm1, 
                                   n.var.arm2, w1.var, w2.var, time.var,
                                   near.pd)
      sendMessage(mComb, "combined", which.run)
      if (mComb$has.error){
        message("- ", crayon::yellow("[!] "), 
                "model could not be fitted using",
                " vcov='complex'. Switching to 'simple'...")
        mCombHasError = TRUE
      }
    }
    if (identical(vcov[1], "simple") ||
        mCombHasError){
      mComb = fitCombinedModel(which.combine, which.combine.var,
                               data, study.var, multi.study, es.var, se.var,
                               mGeneral, .type.es, round.digits, hakn,
                               .raw.bin.es, nnt.cer, rho.within.study,
                               method.tau, method.tau.ci, dots,
                               es.binary.raw.vars, phi.within.study)
      # If model failed, add to error model list
      if (mComb$has.error){
        error.model.list = append(error.model.list, "combined")
        if ("combined" %in% error.model.list){
          mComb$message = mComb$message
        }
      }
      sendMessage(mComb, "combined", which.run)
    }
  } else {
    mComb = NULL
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  5. Outliers Removed                                                      #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  if ("outliers" %in% which.run){
    
    outlier.settings = 
      selectOutlierModel(
        which.run, mGeneral, mComb, which.outliers)
    m.for.outliers = outlier.settings$m.for.outliers
    which.outliers = outlier.settings$which.outliers
    
    mOutliers = fitOutliersModel(data, study.var, multi.study,
                                 mGeneral, .type.es, round.digits,
                                 .raw.bin.es, nnt.cer, which.run,
                                 which.outliers, method.tau,
                                 m.for.outliers)
    # If model failed, add to error model list
    if (mOutliers$has.error){
      error.model.list = append(error.model.list, "outliers")
    }
    sendMessage(mOutliers, "outliers", which.run)
  } else {
    mOutliers = NULL
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  6. Influential Cases                                                     #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  if ("influence" %in% which.run){
    mInfluence = fitInfluenceModel(which.influence, mComb, mGeneral,
                                   which.run, method.tau,
                                   .raw.bin.es, .type.es, round.digits,
                                   nnt.cer)
    influenceRes = mInfluence$influenceRes
    sendMessage(mInfluence, "influence", which.run)
    # If model failed, add to error model list
    if (mInfluence$has.error){
      error.model.list = append(error.model.list, "influence")
    }
  } else {
    mInfluence = NULL
    influenceRes = NULL
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  7. (Low) risk of bias                                                    #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  if ("rob" %in% which.run){
    mRob = fitRobModel(which.run, which.rob, which.outliers,
                       mGeneral, mComb, low.rob.filter, method.tau,
                       .raw.bin.es, .type.es, round.digits, nnt.cer)
    sendMessage(mRob, "rob", which.run)
    # If model failed, add to error model list
    if (mRob$has.error){
      error.model.list = append(error.model.list, "rob")
    }
  } else {
    mRob = NULL
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  8. Three-Level Model                                                     #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  if (method.tau != "REML" & "threelevel" %in% which.run){
    message("- ", crayon::green("[OK] "), 
            "3L-model tau(s) estimated using 'REML', since '",
            method.tau, "' is not applicable...")
  }
  
  if ("threelevel" %in% which.run){
    mThreeLevel = 
      fitThreeLevelModel(data, es.var, se.var, arm.var.1, arm.var.2,
                         measure.var, study.var, .raw.bin.es, .type.es, hakn,
                         method.tau.meta, method.tau.ci, method.tau,
                         dots, es.binary.raw.vars, round.digits,
                         nnt.cer, which.run, mGeneral, mCombined,
                         use.rve, i2.ci.boot, nsim.boot)
    sendMessage(mThreeLevel, "threelevel", which.run)
    # If model failed, add to error model list
    if (mThreeLevel$has.error){
      error.model.list = append(error.model.list, "threelevel")
    }
  } else {
    mThreeLevel = NULL
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  9. Three-Level CHE Model                                                 #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  if (method.tau != "REML" & "threelevel.che" %in% which.run){
    message("- ", crayon::green("[OK] "), 
            "Three-level CHE model tau(s) estimated using 'REML', since '",
            method.tau, "' is not applicable...")
  }
  
  if ("threelevel.che" %in% which.run){
    mCHEHasError = FALSE
    if (identical(vcov[1], "complex")){
      mCHE = 
        fitThreeLevelHACEModel(data, es.var, se.var, arm.var.1, arm.var.2,
                               measure.var, study.var, .raw.bin.es, .type.es, hakn,
                               method.tau.meta, method.tau.ci, method.tau,
                               dots, es.binary.raw.vars, round.digits,
                               nnt.cer, which.run, mGeneral, mCombined,
                               use.rve, rho.within.study, which.combine.var,
                               phi.within.study, n.var.arm1, 
                               n.var.arm2, w1.var, w2.var, time.var,
                               near.pd, i2.ci.boot, nsim.boot)
      sendMessage(mCHE, "threelevel.che", which.run)
      if (mCHE$has.error){
        message("- ", crayon::yellow("[!] "), 
                "model could not be fitted using",
                " vcov='complex'. Switching to 'simple'...")
        mCHEHasError = TRUE
      }
    } 
    if (identical(vcov[1], "simple") ||
        mCHEHasError){
      mCHE = 
        fitThreeLevelCHEModel(data, es.var, se.var, arm.var.1, arm.var.2,
                              measure.var, study.var, .raw.bin.es, .type.es, hakn,
                              method.tau.meta, method.tau.ci, method.tau,
                              dots, es.binary.raw.vars, round.digits,
                              nnt.cer, which.run, mGeneral, mCombined,
                              use.rve, rho.within.study, phi.within.study,
                              i2.ci.boot, nsim.boot)
      sendMessage(mCHE, "threelevel.che", which.run)
      # If model failed, add to error model list
      if (mCHE$has.error){
        error.model.list = append(error.model.list, "threelevel.che")
      }
    }
  } else {
    mCHE = NULL
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  RETURN                                                                   #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  
  # Combine everything
  rbind(mGeneral$res, mLowest$res, 
        mHighest$res, mOutliers$res,
        mInfluence$res, mRob$res, 
        mComb$res, mThreeLevel$res, mCHE$res) -> summary
  
  # Fill NAs 
  summary[summary == "[NA; NA]"] = "[-; -]"
  summary[is.na(summary)] = "-"
  if (identical(.type.es, "RR")){
    colnames(summary)[2:3] = c("rr", "rr.ci")
  }
  
  # Define ES type for all models
  mGeneral$m$.type.es = .type.es
  mComb$m$.type.es = .type.es
  mLowest$m$.type.es = .type.es
  mHighest$m$.type.es = .type.es
  mOutliers$m$.type.es = .type.es
  mInfluence$m$.type.es = .type.es
  mRob$m$.type.es = .type.es
  mThreeLevel$m$.type.es = .type.es
  mCHE$m$.type.es = .type.es
  
  # Return
  list(summary = summary,
       model.overall = mGeneral$m,
       model.combined = mComb$m,
       model.lowest = mLowest$m,
       model.highest = mHighest$m,
       model.outliers = mOutliers$m,
       model.influence = mInfluence$m,
       model.rob = mRob$m,
       model.threelevel = mThreeLevel$m,
       model.threelevel.var.comp = mThreeLevel$m$variance.components,
       model.threelevel.che = mCHE$m,
       model.threelevel.che.var.comp = mCHE$m$variance.components,
       influence.analysis = influenceRes,
       which.run = which.run,
       data = data.original,
       html = html,
       round.digits = round.digits,
       nnt.cer = nnt.cer,
       use.rve = use.rve,
       call = match.call(),
       .type.es = .type.es,
       .raw.bin.es = .raw.bin.es,
       error.model.list = error.model.list,
       i2.ci.boot = i2.ci.boot,
       nsim.boot = nsim.boot) -> returnlist
  class(returnlist) = c("runMetaAnalysis", "list")
  
  if ("lowest.highest" %in% which.run){
    message("- ", crayon::green("[OK] "), "Done!")
  }
  if (!identical(as.numeric(length(error.model.list)), 0) ||
      warn.end){
    warning("There were some issues during the ",
            "calculations. Please check the report above.", 
            call. = FALSE)
  }
  return(returnlist)
}


#' Print method for objects of class 'runMetaAnalysis'
#'
#' Print S3 method for objects of class \code{runMetaAnalysis}.
#'
#' @param x An object of class \code{runMetaAnalysis}.
#' @param ... Additional arguments.
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers 
#' \email{p.cuijpers@@vu.nl}
#'
#' @importFrom knitr kable
#' @importFrom dplyr as_tibble
#' @importFrom kableExtra kable_styling column_spec footnote
#' @importFrom crayon green blue magenta bold
#' @importFrom stringr str_sub
#'
#' @export
#' @method print runMetaAnalysis


print.runMetaAnalysis = function(x, ...){
  
  if (class(x)[2] == "correctPublicationBias"){
    
    # # # # # # # # # # # # # #
    # If Pubbias is corrected #
    # # # # # # # # # # # # # # 
    
    
    models = list("overall" = "Overall", 
                  "lowest.highest" = c("One ES/study (lowest)",
                                       "One ES/study (highest)"), 
                  "outliers" = "Outliers removed",
                  "influence" = "Influence Analysis", 
                  "rob" = "rob", 
                  "combined" = "Combined", 
                  "threelevel" = "Three-Level Model",
                  "threelevel.che" = "Three-Level Model (CHE)")
    
    run.models = unlist(models[x$which.run])
    if ("rob" %in% run.models){
      run.models[["rob"]] = 
        rownames(x$summary)[grepl("Only", rownames(x$summary))]
    }
    
    cat(crayon::blue$bold("Model results "))
    cat(crayon::blue$bold(
      "------------------------------------------------ \n"))
    dat = x$summary[run.models, 1:8]
    tbl = dplyr::as_tibble(cbind(Model = rownames(dat), dat))
    old = options(pillar.bold=TRUE)
    tbl.format = format(tbl)[-c(1,3)]
    tbl.format[-1] = lapply(tbl.format[-1], 
                            function(x) stringr::str_sub(x, 19))
    tbl.format[1] = stringr::str_sub(tbl.format[1], 3)
    cat(do.call(c, tbl.format), sep="\n")
    options(old)
    
    # Add publication bias corrected results
    cat("\n")
    cat(crayon::blue$bold(paste0("Publication bias correction ('", 
                                 x$correctPublicationBias$which.run,
                                 "' model) ")))
    cat(crayon::blue$bold("----------------------- \n"))
    dat.pb = x$correctPublicationBias$summary[,1:8]
    tbl.pb = dplyr::as_tibble(cbind(Model = rownames(dat.pb), dat.pb))
    old = options(pillar.bold=TRUE)
    tbl.format.pb = format(tbl.pb)[-c(1,3)]
    tbl.format.pb[-1] = lapply(tbl.format.pb[-1], 
                               function(x) stringr::str_sub(x, 19))
    tbl.format.pb[1] = stringr::str_sub(tbl.format.pb[1], 3)
    cat(do.call(c, tbl.format.pb), sep="\n")
    options(old)
    
    
    if ("threelevel" %in% x$which.run){
      cat("\n")
      cat(crayon::blue$bold("Variance components (three-level model) "))
      cat(crayon::blue$bold("---------------------- \n"))
      if (is.na(x$model.threelevel.var.comp)[1]){
        cat("-")
      } else {
        print(x$model.threelevel.var.comp)
      }
    }
    
    if ("threelevel.che" %in% x$which.run){
      cat("\n")
      cat(crayon::blue$bold("Variance components (three-level CHE model) "))
      cat(crayon::blue$bold("------------------ \n"))
      if (is.na(x$model.threelevel.che.var.comp)[1]){
        cat("-")
      } else {
        print(x$model.threelevel.che.var.comp)
      }
    }
    
    
    if (x$html == TRUE){
      
      x$summary = x$summary[run.models,]
      
      # Add publication bias
      pub.row = c(rep(" ", ncol(x$summary)-1),
                  paste0("Corrections were applied to the '",
                         x$correctPublicationBias$which.run, 
                         "' model."))
      rownames(x$correctPublicationBias$summary) = 
        paste("-", rownames(x$correctPublicationBias$summary))
      rbind(x$summary, "<i>Publication bias correction</i>" = pub.row,
            x$correctPublicationBias$summary) -> x$summary
      
      # Add footnote labels
      fn.rows = x$summary$excluded != "none"
      rownames(x$summary)[fn.rows] = paste0(rownames(x$summary[fn.rows,]), "<sup>",
                                            letters[1:nrow(x$summary[fn.rows,])], "</sup>")
      
      if (identical(x$summary$excluded[fn.rows], 
                    character(0))) {
        footnotes = "none"
      } else {
        footnotes = x$summary$excluded[fn.rows]
      }
      
      x$summary %>%
        {.$Analysis = rownames(.); rownames(.) = NULL; .} %>%
        dplyr::select(Analysis, dplyr::everything(), -excluded) %>%
        setColnames(c(".", "<i>k</i>", 
                      ifelse(identical(x$.type.es, "RR"), 
                             "<i>RR</i>",
                             "<i>g</i>"), 
                      "CI", "<i>p</i>",
                      "<i>I</i><sup>2</sup>",
                      "CI", "PI", "NNT")) %>%
        knitr::kable(escape = FALSE) %>%
        kableExtra::kable_styling(font_size = 8, full_width = FALSE) %>%
        kableExtra::column_spec(1, bold = TRUE, width_min = "13em") %>%
        kableExtra::footnote(general = "Excluded effect sizes/studies:",
                             alphabet = footnotes) %>%
        print()
    }
    
    
  } else {
    
    # # # # # # # # # # # # # #
    # Standard Result Class   #
    # # # # # # # # # # # # # # 
    
    models = list("overall" = "Overall", 
                  "lowest.highest" = c("One ES/study (lowest)",
                                       "One ES/study (highest)"), 
                  "outliers" = "Outliers removed",
                  "influence" = "Influence Analysis", 
                  "rob" = "rob", 
                  "combined" = "Combined", 
                  "threelevel" = "Three-Level Model",
                  "threelevel.che" = "Three-Level Model (CHE)")
    
    run.models = unlist(models[x$which.run])
    if ("rob" %in% run.models){
      run.models[["rob"]] = 
        rownames(x$summary)[grepl("Only", rownames(x$summary))]
    }
    
    cat(crayon::blue$bold("Model results "))
    cat(crayon::blue$bold(
      "------------------------------------------------ \n"))
    dat = x$summary[run.models, 1:8]
    tbl = dplyr::as_tibble(cbind(Model = rownames(dat), dat))
    old = options(pillar.bold=TRUE)
    tbl.format = format(tbl)[-c(1,3)]
    tbl.format[-1] = lapply(tbl.format[-1], 
                            function(x) stringr::str_sub(x, 19))
    tbl.format[1] = stringr::str_sub(tbl.format[1], 3)
    cat(do.call(c, tbl.format), sep="\n")
    options(old)
    
    if ("threelevel" %in% x$which.run){
      cat("\n")
      cat(crayon::blue$bold("Variance components (three-level model) "))
      cat(crayon::blue$bold("---------------------- \n"))
      if (is.na(x$model.threelevel.var.comp)[1]){
        cat("-")
      } else {
        print(x$model.threelevel.var.comp) 
      }
    }
    
    if ("threelevel.che" %in% x$which.run){
      cat("\n")
      cat(crayon::blue$bold("Variance components (three-level CHE model) "))
      cat(crayon::blue$bold("------------------ \n"))
      if (is.na(x$model.threelevel.che.var.comp)[1]){
        cat("-")
      } else {
        print(x$model.threelevel.che.var.comp) 
      }
    }
    
    if (x$html == TRUE){
      
      x$summary = x$summary[run.models,]
      # Add footnote labels
      fn.rows = x$summary$excluded != "none"
      rownames(x$summary)[fn.rows] = paste0(rownames(x$summary[fn.rows,]), "<sup>",
                                            letters[1:nrow(x$summary[fn.rows,])], "</sup>")
      
      if (identical(x$summary$excluded[fn.rows], 
                    character(0))) {
        footnotes = "none"
      } else {
        footnotes = x$summary$excluded[fn.rows]
      }
      
      x$summary %>%
        {.$Analysis = rownames(.); rownames(.) = NULL; .} %>%
        dplyr::select(Analysis, dplyr::everything(), -excluded) %>%
        setColnames(c(".", "<i>k</i>", 
                      ifelse(identical(x$.type.es, "RR"), 
                             "<i>RR</i>",
                             "<i>g</i>"), 
                      "CI", "<i>p</i>",
                      "<i>I</i><sup>2</sup>",
                      "CI", "PI", "NNT")) %>%
        knitr::kable(escape = FALSE) %>%
        kableExtra::kable_styling(font_size = 8, full_width = FALSE) %>%
        kableExtra::column_spec(1, bold = TRUE, width_min = "13em") %>%
        kableExtra::footnote(general = "Excluded effect sizes/studies:",
                             alphabet = footnotes) %>%
        print()
    }
    
  }
  
}


#' Print method for objects of class 'runMetaAnalysis'
#'
#' Print S3 method for objects of class \code{runMetaAnalysis}.
#'
#' @param x An object of class \code{runMetaAnalysis}.
#' @param ... Additional arguments.
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers 
#' \email{p.cuijpers@@vu.nl}
#'
#' @importFrom knitr kable
#' @importFrom dplyr as_tibble
#' @importFrom kableExtra kable_styling column_spec footnote
#' @importFrom crayon green blue magenta bold
#' @importFrom stringr str_sub
#'
#' @export
#' @method print runMetaAnalysis


print.runMetaAnalysis = function(x, ...){
  
  if (class(x)[2] == "correctPublicationBias"){
    
    # # # # # # # # # # # # # #
    # If Pubbias is corrected #
    # # # # # # # # # # # # # # 
    
    
    models = list("overall" = "Overall", 
                  "lowest.highest" = c("One ES/study (lowest)",
                                       "One ES/study (highest)"), 
                  "outliers" = "Outliers removed",
                  "influence" = "Influence Analysis", 
                  "rob" = "rob", 
                  "combined" = "Combined", 
                  "threelevel" = "Three-Level Model",
                  "threelevel.che" = "Three-Level Model (CHE)")
    
    run.models = unlist(models[x$which.run])
    if ("rob" %in% run.models){
      run.models[["rob"]] = 
        rownames(x$summary)[grepl("Only", rownames(x$summary))]
    }
    
    cat(crayon::blue$bold("Model results "))
    cat(crayon::blue$bold(
      "------------------------------------------------ \n"))
    dat = x$summary[run.models, 1:8]
    tbl = dplyr::as_tibble(cbind(Model = rownames(dat), dat))
    old = options(pillar.bold=TRUE)
    tbl.format = format(tbl)[-c(1,3)]
    tbl.format[-1] = lapply(tbl.format[-1], 
                            function(x) stringr::str_sub(x, 19))
    tbl.format[1] = stringr::str_sub(tbl.format[1], 3)
    cat(do.call(c, tbl.format), sep="\n")
    options(old)
    
    # Add publication bias corrected results
    cat("\n")
    cat(crayon::blue$bold(paste0("Publication bias correction ('", 
               x$correctPublicationBias$which.run,
               "' model) ")))
    cat(crayon::blue$bold("----------------------- \n"))
    dat.pb = x$correctPublicationBias$summary[,1:8]
    tbl.pb = dplyr::as_tibble(cbind(Model = rownames(dat.pb), dat.pb))
    old = options(pillar.bold=TRUE)
    tbl.format.pb = format(tbl.pb)[-c(1,3)]
    tbl.format.pb[-1] = lapply(tbl.format.pb[-1], 
                            function(x) stringr::str_sub(x, 19))
    tbl.format.pb[1] = stringr::str_sub(tbl.format.pb[1], 3)
    cat(do.call(c, tbl.format.pb), sep="\n")
    options(old)
    
    
    if ("threelevel" %in% x$which.run){
      cat("\n")
      cat(crayon::blue$bold("Variance components (three-level model) "))
      cat(crayon::blue$bold("---------------------- \n"))
      if (is.na(x$model.threelevel.var.comp)[1]){
        cat("-")
      } else {
        dat = x$model.threelevel.var.comp
        dat = within(dat,{tau2 = round(tau2, x$round.digits+1)})
        tbl = dplyr::as_tibble(cbind(Source = rownames(dat), dat))
        old = options(pillar.bold=TRUE)
        tbl.format = format(tbl)[-c(1,3)]
        tbl.format[-1] = lapply(tbl.format[-1], 
                                function(x) stringr::str_sub(x, 19))
        tbl.format[1] = stringr::str_sub(tbl.format[1], 3)
        cat(do.call(c, tbl.format), sep="\n")
        options(old)
      }
    }
    
    if ("threelevel.che" %in% x$which.run){
      cat("\n")
      cat(crayon::blue$bold("Variance components (three-level CHE model) "))
      cat(crayon::blue$bold("------------------ \n"))
      if (is.na(x$model.threelevel.che.var.comp)[1]){
        cat("-")
      } else {
        dat = x$model.threelevel.che.var.comp
        dat = within(dat,{tau2 = round(tau2, x$round.digits+1)})
        tbl = dplyr::as_tibble(cbind(Source = rownames(dat), dat))
        old = options(pillar.bold=TRUE)
        tbl.format = format(tbl)[-c(1,3)]
        tbl.format[-1] = lapply(tbl.format[-1], 
                                function(x) stringr::str_sub(x, 19))
        tbl.format[1] = stringr::str_sub(tbl.format[1], 3)
        cat(do.call(c, tbl.format), sep="\n")
        options(old)
      }
    }
    
    
    if (x$html == TRUE){
      
      x$summary = x$summary[run.models,]
      
      # Add publication bias
      pub.row = c(rep(" ", ncol(x$summary)-1),
                  paste0("Corrections were applied to the '",
                         x$correctPublicationBias$which.run, 
                         "' model."))
      rownames(x$correctPublicationBias$summary) = 
        paste("-", rownames(x$correctPublicationBias$summary))
      rbind(x$summary, "<i>Publication bias correction</i>" = pub.row,
            x$correctPublicationBias$summary) -> x$summary
      
      # Add footnote labels
      fn.rows = x$summary$excluded != "none"
      rownames(x$summary)[fn.rows] = paste0(rownames(x$summary[fn.rows,]), "<sup>",
                                            letters[1:nrow(x$summary[fn.rows,])], "</sup>")
      
      if (identical(x$summary$excluded[fn.rows], 
                    character(0))) {
        footnotes = "none"
      } else {
        footnotes = x$summary$excluded[fn.rows]
      }
      
      x$summary %>%
        {.$p = ifelse(.$p == "<0.001", "&lt;0.001", .$p);.} %>% 
        {.$Analysis = rownames(.); rownames(.) = NULL; .} %>%
        dplyr::select(Analysis, dplyr::everything(), -excluded) %>%
        setColnames(c(".", "<i>k</i>", 
                                 ifelse(identical(x$.type.es, "RR"), 
                                        "<i>RR</i>",
                                        "<i>g</i>"), 
                                 "CI", "<i>p</i>",
                                 "<i>I</i><sup>2</sup>",
                                 "CI", "PI", "NNT")) %>%
        knitr::kable(escape = FALSE) %>%
        kableExtra::kable_styling(font_size = 8, full_width = FALSE) %>%
        kableExtra::column_spec(1, bold = TRUE, width_min = "13em") %>%
        kableExtra::footnote(general = "Excluded effect sizes/studies:",
                             alphabet = footnotes) %>%
        print()
    }
    

  } else {
    
    # # # # # # # # # # # # # #
    # Standard Result Class   #
    # # # # # # # # # # # # # # 
    
    models = list("overall" = "Overall", 
                  "lowest.highest" = c("One ES/study (lowest)",
                                       "One ES/study (highest)"), 
                  "outliers" = "Outliers removed",
                  "influence" = "Influence Analysis", 
                  "rob" = "rob", 
                  "combined" = "Combined", 
                  "threelevel" = "Three-Level Model",
                  "threelevel.che" = "Three-Level Model (CHE)")
    
    run.models = unlist(models[x$which.run])
    if ("rob" %in% run.models){
      run.models[["rob"]] = 
        rownames(x$summary)[grepl("Only", rownames(x$summary))]
    }
    
    cat(crayon::blue$bold("Model results "))
    cat(crayon::blue$bold(
      "------------------------------------------------ \n"))
    dat = x$summary[run.models, 1:8]
    tbl = dplyr::as_tibble(cbind(Model = rownames(dat), dat))
    old = options(pillar.bold=TRUE)
    tbl.format = format(tbl)[-c(1,3)]
    tbl.format[-1] = lapply(tbl.format[-1], 
                            function(x) stringr::str_sub(x, 19))
    tbl.format[1] = stringr::str_sub(tbl.format[1], 3)
    cat(do.call(c, tbl.format), sep="\n")
    options(old)
    
    if ("threelevel" %in% x$which.run){
      cat("\n")
      cat(crayon::blue$bold("Variance components (three-level model) "))
      cat(crayon::blue$bold("---------------------- \n"))
      if (is.na(x$model.threelevel.var.comp)[1]){
        cat("-")
      } else {
        dat = x$model.threelevel.var.comp
        dat = within(dat,{tau2 = round(tau2, x$round.digits+1)})
        tbl = dplyr::as_tibble(cbind(Source = rownames(dat), dat))
        old = options(pillar.bold=TRUE)
        tbl.format = format(tbl)[-c(1,3)]
        tbl.format[-1] = lapply(tbl.format[-1], 
                                function(x) stringr::str_sub(x, 19))
        tbl.format[1] = stringr::str_sub(tbl.format[1], 3)
        cat(do.call(c, tbl.format), sep="\n")
        options(old)
      }
    }
    
    if ("threelevel.che" %in% x$which.run){
      cat("\n")
      cat(crayon::blue$bold("Variance components (three-level CHE model) "))
      cat(crayon::blue$bold("------------------ \n"))
      if (is.na(x$model.threelevel.che.var.comp)[1]){
        cat("-")
      } else {
        dat = x$model.threelevel.che.var.comp
        dat = within(dat,{tau2 = round(tau2, x$round.digits+1)})
        tbl = dplyr::as_tibble(cbind(Source = rownames(dat), dat))
        old = options(pillar.bold=TRUE)
        tbl.format = format(tbl)[-c(1,3)]
        tbl.format[-1] = lapply(tbl.format[-1], 
                                function(x) stringr::str_sub(x, 19))
        tbl.format[1] = stringr::str_sub(tbl.format[1], 3)
        cat(do.call(c, tbl.format), sep="\n")
        options(old)
      }
    }
    
    if (x$html == TRUE){
      
      x$summary = x$summary[run.models,]
      # Add footnote labels
      fn.rows = x$summary$excluded != "none"
      rownames(x$summary)[fn.rows] = paste0(rownames(x$summary[fn.rows,]), "<sup>",
                                            letters[1:nrow(x$summary[fn.rows,])], "</sup>")
      
      if (identical(x$summary$excluded[fn.rows], 
                    character(0))) {
        footnotes = "none"
      } else {
        footnotes = x$summary$excluded[fn.rows]
      }
      
      x$summary %>%
        {.$p = ifelse(.$p == "<0.001", "&lt;0.001", .$p);.} %>% 
        {.$Analysis = rownames(.); rownames(.) = NULL; .} %>%
        dplyr::select(Analysis, dplyr::everything(), -excluded) %>%
        setColnames(c(".", "<i>k</i>", 
                                 ifelse(identical(x$.type.es, "RR"), 
                                        "<i>RR</i>",
                                        "<i>g</i>"), 
                                 "CI", "<i>p</i>",
                                 "<i>I</i><sup>2</sup>",
                                 "CI", "PI", "NNT")) %>%
        knitr::kable(escape = FALSE) %>%
        kableExtra::kable_styling(font_size = 8, full_width = FALSE) %>%
        kableExtra::column_spec(1, bold = TRUE, width_min = "13em") %>%
        kableExtra::footnote(general = "Excluded effect sizes/studies:",
                             alphabet = footnotes) %>%
        print()
    }

  }

}

#' Plot method for objects of class 'runMetaAnalysis'
#'
#' Plot S3 method for objects of class \code{runMetaAnalysis}.
#'
#' @param x An object of class \code{runMetaAnalysis}.
#' @param which Model to be plotted. Can be one of \code{"overall"},
#' \code{"combined"}, \code{"lowest.highest"}, \code{"outliers"},
#' \code{"influence"}, \code{"threelevel"}, \code{"threelevel.che"},
#' \code{"baujat"}, \code{"loo-es"}, \code{"loo-i2"},
#' \code{"trimfill"}, \code{"limitmeta"} or \code{"selection"}.
#' @param eb Prints a forest plot with empirical Bayes point estimates and study-specific
#' prediction intervals as proposed by van Aert (2021). Defaults to `FALSE`.
#' @param eb.labels If `eb` is `TRUE`, should the empirical Bayes estimates 
#' and prediction intervals for each study be printed on the right side of the 
#' forest plot? Defaults to `FALSE`.
#' @param ... Additional arguments.
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @importFrom meta metagen forest.meta
#' @importFrom metafor forest.rma
#' @importFrom stringr str_replace_all
#' @importFrom purrr map
#' @importFrom dplyr mutate
#' @importFrom metasens funnel.limitmeta
#' @importFrom meta funnel.meta
#' @importFrom metafor plot.rma.uni.selmodel
#' @importFrom crayon green yellow cyan bold
#'
#' @export
#' @method plot runMetaAnalysis
#' @references 
#' van Aert, R. C., Schmid, C. H., Svensson, D., & Jackson, D. (2021). 
#' Study specific prediction intervals for random-effects meta-analysis: A tutorial.
#'  _Research Synthesis Methods, 12_(4), 429-447.

plot.runMetaAnalysis = function(x, which = NULL, eb = FALSE, 
                                eb.labels = FALSE, ...){

  if (!is.logical(eb[1])){
    stop("'eb' must be either TRUE or FALSE.")
  }
  
  if (!eb[1]){
    models = list("overall" = "model.overall",
                  "lowest.highest" = c("model.lowest", "model.highest"),
                  "outliers" = "model.outliers", 
                  "influence" = "model.influence",
                  "rob" = "model.rob", 
                  "combined" = "model.combined",
                  "threelevel" = "model.threelevel",
                  "che" = "model.threelevel.che")
    
    models.which.run = list("overall" = "overall", 
                            "combined" = "combined",
                            "lowest.highest" = "lowest.highest", 
                            "outliers" = "outliers",
                            "influence" = "influence", 
                            "rob" = "rob", 
                            "threelevel" = "threelevel",
                            "che" = "threelevel.che",
                            "threelevel.che" = "threelevel.che")
    
    if (!is.null(which)){
      if (!which %in% names(models.which.run[which])){
        if (!which %in% c("selection", "baujat",
                          "loo", "loo-es", "loo-i2", "trimfill",
                          "limitmeta", "summary")){
          stop("Model not available for plotting.")
        }
      } 
    }
    
    leftCols = c("studlab", "comparison.only", "instrument", "TE", "seTE")
    leftColsRR = c("studlab", "comparison.only", "instrument")
    leftLabs = c("Study", "Comparison", "Instrument", "g", "S.E.")
    leftLabsRR = c("Study", "Comparison", "Instrument")
    
    # print forest plot by default
    if (is.null(which)){
      if (models[[x$which.run[1]]][1] != "model.threelevel" &&
          models[[x$which.run[1]]][1] != "model.threelevel.che"){
        message("- ", crayon::green("[OK] "), 
                "Generating forest plot ('", x$which.run[1], "' model).")
        if (identical(x$.type.es, "RR")){
          meta::forest.meta(x[[models[[x$which.run[1]]][1]]],
                            layout = "JAMA", ...)
        } else {
          meta::forest.meta(x[[models[[x$which.run[1]]][1]]], 
                            layout = "JAMA", ...)
        }
      }
      if (models[[x$which.run[1]]][1] == "lowest.highest"){
        message("- ", crayon::green("[OK] "), "Generating forest plot ('", 
                "highest", "' model).")
        if (identical(x$.type.es, "RR")){
          meta::forest.meta(x$model.highest, 
                            layout = "JAMA", ...)
        } else {
          meta::forest.meta(x$model.highest, 
                            layout = "JAMA", ...)
        }
      }
      if (models[[x$which.run[1]]][1] == "model.threelevel"){
        if (identical(x$.type.es, "RR")){
          metafor::forest.rma(x$model.threelevel, 
                              transf = exp, ...)
        } else {
          metafor::forest.rma(x$model.threelevel, ...)
        }
      }
      if (models[[x$which.run[1]]][1] == "model.threelevel.che"){
        if (identical(x$.type.es, "RR")){
          metafor::forest.rma(x$model.threelevel.che, 
                              transf = exp, ...)
        } else {
          metafor::forest.rma(x$model.threelevel.che, ...)
        }
      }
      
    } else {
      
      if (which[1] == "overall"){
        message("- ", crayon::green("[OK] "), 
                "Generating forest plot ('overall' model).")
        if (identical(x$.type.es, "RR")){
          meta::forest.meta(x$model.overall, 
                            layout = "JAMA",  ...)
        } else {
          meta::forest.meta(x$model.overall, 
                            layout = "JAMA", ...)
        }
      }
      
      if (which[1] == "lowest.highest"){
        message("- ", crayon::green("[OK] "), 
                "Generating forest plot ('lowest' model).")
        if (identical(x$.type.es, "RR")){
          meta::forest.meta(x$model.lowest, 
                            layout = "JAMA", ...)
        } else {
          meta::forest.meta(x$model.lowest, 
                            layout = "JAMA", ...)
        }
        message("- ", crayon::green("[OK] "), 
                "Generating forest plot ('highest' model).")
        if (identical(x$.type.es, "RR")){
          meta::forest.meta(x$model.highest, 
                            layout = "JAMA", ...)
        } else {
          meta::forest.meta(x$model.highest, 
                            layout = "JAMA", ...)
        }
      }
      
      if (which[1] == "outliers"){
        message("- ", crayon::green("[OK] "), 
                "Generating forest plot ('outliers' model).")
        if (identical(x$.type.es, "RR")){
          meta::forest.meta(x$model.outliers, 
                            layout = "JAMA", ...)
        } else {
          meta::forest.meta(x$model.outliers, 
                            layout = "JAMA", ...)
        }
      }
      
      if (which[1] == "influence"){
        message("- ", crayon::green("[OK] "), 
                "Generating forest plot ('influence' model).")
        if (identical(x$.type.es, "RR")){
          meta::forest.meta(x$model.influence, 
                            layout = "JAMA", ...)
        } else {
          meta::forest.meta(x$model.influence, 
                            layout = "JAMA", ...)
        }
      }
      
      if (which[1] == "rob"){
        message("- ", crayon::green("[OK] "), 
                "Generating forest plot ('rob' model).")
        if (identical(x$.type.es, "RR")){
          meta::forest.meta(x$model.rob, 
                            layout = "JAMA", ...)
        } else {
          meta::forest.meta(x$model.rob, 
                            layout = "JAMA", ...)
        }
      }
      
      if (which[1] == "combined"){
        message("- ", crayon::green("[OK] "), 
                "Generating forest plot ('combined' model).")
        if (identical(x$.type.es, "RR")){
          meta::forest.meta(x$model.combined, 
                            layout = "JAMA", ...)
        } else {
          meta::forest.meta(x$model.combined, 
                            layout = "JAMA", ...)
        }
      }
      
      if (which[1] == "threelevel"){
        message("- ", crayon::green("[OK] "), 
                "Generating forest plot ('threelevel' model).")
        if (identical(x$.type.es, "RR")){
          metafor::forest.rma(x$model.threelevel, transf = exp, ...)
        } else {
          metafor::forest.rma(x$model.threelevel, ...)
        }
      }
      
      if (which[1] == "che"){
        message("- ", crayon::green("[OK] "), 
                "Generating forest plot ('threelevel.che' model).")
        if (identical(x$.type.es, "RR")){
          metafor::forest.rma(x$model.threelevel.che, transf = exp, ...)
        } else {
          metafor::forest.rma(x$model.threelevel.che, ...)
        }
      }
      
      if (which[1] == "threelevel.che"){
        message("- ", crayon::green("[OK] "), 
                "Generating forest plot ('threelevel.che' model).")
        if (identical(x$.type.es, "RR")){
          metafor::forest.rma(x$model.threelevel.che, transf = exp, ...)
        } else {
          metafor::forest.rma(x$model.threelevel.che, ...)
        }
      }
      
      if (which[1] == "baujat"){
        if (!"influence" %in% x$which.run){
          stop("Baujat plots are only available when influence analyses have",
               " been conducted.")
        }
        message("- ", crayon::green("[OK] "), 
                "Generating baujat plot.")
        plot(x$influence.analysis$BaujatPlot)
      }
      
      if (which[1] == "loo" | which[1] == "loo-es"){
        if (!"influence" %in% x$which.run){
          stop("L-O-O plots are only available when influence analyses have",
               " been conducted.")
        }
        message("- ", crayon::green("[OK] "), 
                "Generating leave-one-out forest plot.")
        suppressWarnings({plot(x$influence.analysis$ForestEffectSize)})
      }
      
      if (which[1] == "loo-i2"){
        if (!"influence" %in% x$which.run){
          stop("L-O-O plots are only available when influence analyses have",
               " been conducted.")
        }
        message("- ", crayon::green("[OK] "), 
                "Generating leave-one-out forest plot (sorted by I2).")
        suppressWarnings({plot(x$influence.analysis$ForestI2)})
      }
      
      if (which[1] == "trimfill"){
        if (is.null(x$correctPublicationBias)){
          stop("No publication bias analysis results found.",
               " Have you called correctPublicationBias()?")
        }
        message("- ", crayon::green("[OK] "), 
                "Generating funnel plot for the trim-and-fill analysis.")
        suppressWarnings({meta::funnel.meta(x$correctPublicationBias$model.trimfill)})
      }
      
      if (which[1] == "limitmeta"){
        if (is.null(x$correctPublicationBias)){
          stop("No publication bias analysis results found.",
               " Have you called correctPublicationBias()?")
        }
        message("- ", crayon::green("[OK] "), 
                "Generating funnel plot for the limit meta-analysis.")
        suppressWarnings({
          metasens::funnel.limitmeta(
            x$correctPublicationBias$model.limitmeta)})
      }
      
      if (which[1] == "selection"){
        if (is.null(x$correctPublicationBias)){
          stop("No publication bias analysis results found.",
               " Have you called correctPublicationBias()?")
        }
        if (inherits(x$correctPublicationBias$model.selection, "try-error")){
          stop("Selection model could not be calculated.", 
               " Plot cannot be generated.")
        } else {
          message("- ", crayon::green("[OK] "), 
                  "Generating selection model likelihood plot.")
          suppressWarnings({
            metafor::plot.rma.uni.selmodel(x$correctPublicationBias$model.selection,
                                           xlim = c(0, 0.2))})
        }
      }
      
      if (which[1] == "summary"){
        
        models = list("overall" = 1, "lowest.highest" = c(2,3), "outliers" = 4,
                      "influence" = 5, "rob" = 6, "combined" = 7, "threelevel" = 8,
                      "threelevel.che" = 9)
        
        x$summary = x$summary[unlist(models[x$which.run]),]
        
        if (x$.type.es == "RR"){
          stringr::str_replace_all(x$summary$rr.ci, ";|\\]|\\[", "") %>%
            strsplit(" ") %>% 
            purrr::map(~as.numeric(.)) %>% do.call(rbind,.) %>%
            {colnames(.) = c("lower", "upper");.} %>%
            cbind(model = rownames(x$summary), rr = x$summary$rr,
                  i2 = round(x$summary$i2,1) %>% format(1), .) %>%
            data.frame() %>%
            dplyr::mutate(rr = as.numeric(rr) %>% log(),
                          lower = as.numeric(lower) %>% log(),
                          upper = as.numeric(upper) %>% log()) %>%
            meta::metagen(TE = rr, 
                          lower = lower - 1e-50, 
                          upper = upper + 1e-50, 
                          studlab = model, sm = "RR",
                          data = .) %>%
            meta::forest.meta(
              col.square = "lightblue",
              rightcols = FALSE,
              overall.hetstat = FALSE,
              weight.study = "same",
              test.overall = FALSE, overall = FALSE,
              leftcols = c("studlab", "effect", "ci", "i2"),
              leftlabs = c(expression(bold(Model)), 
                           expression(bold(RR)),
                           expression(bold(CI)),
                           expression(bold(italic(I)^2)))) %>%
            suppressWarnings()
        } else {
          stringr::str_replace_all(x$summary$g.ci, ";|\\]|\\[", "") %>%
            strsplit(" ") %>% purrr::map(~as.numeric(.)) %>% do.call(rbind,.) %>%
            {colnames(.) = c("lower", "upper");.} %>%
            cbind(model = rownames(x$summary), g = x$summary$g,
                  i2 = round(x$summary$i2,1) %>% format(1), .) %>%
            data.frame() %>%
            dplyr::mutate(g = as.numeric(g),
                          lower = as.numeric(lower) - 1e-50,
                          upper = as.numeric(upper) + 1e-50) %>%
            meta::metagen(TE = g, lower = lower, upper = upper, studlab = model,
                          data = .) %>%
            meta::forest.meta(col.square = "lightblue",
                              rightcols = FALSE,
                              overall.hetstat = FALSE,
                              weight.study = "same",
                              test.overall = FALSE, overall = FALSE,
                              leftcols = c("studlab", "TE", "ci", "i2"),
                              leftlabs = c(expression(bold(Model)), expression(bold(g)),
                                           expression(bold(CI)),
                                           expression(bold(italic(I)^2)))) %>%
            suppressWarnings()
        }
      }
    }
  } else {
    dots = list(...)
    argsList = append(list(model = x, which = which, 
                           eb.labels = eb.labels), dots)
    do.call(forestBlup, argsList)
  }
}




#' Show details of 'runMetaAnalysis' class objects
#'
#' S3 method showing analysis settings for objects of class \code{runMetaAnalysis}.
#'
#' @param object An object of class \code{runMetaAnalysis}.
#' @param forest \code{logical}. Should a summary forest plot be returned? \code{TRUE} by default.
#' @param ... Additional arguments.
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @importFrom knitr kable
#' @importFrom meta metagen forest.meta
#' @importFrom dplyr as_tibble
#' @importFrom kableExtra kable_styling column_spec footnote
#' @importFrom stringr str_replace_all
#' @importFrom purrr map
#' @importFrom crayon green blue magenta
#'
#' @export
#' @method summary runMetaAnalysis


summary.runMetaAnalysis = function(object, forest = TRUE, ...){

  x = object
  cat("\n")
  cat("Analysis settings ")
  cat("---------------------------------------------------------- \n")
  cat("\n")

  if (x$model.overall$common == TRUE & x$model.overall$random == FALSE){

    cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Overall]")),
        "Effects pooled using inverse variance weighting. \n")
    
    if (isTRUE(x$.raw.bin.es) & identical(x$.type.es, "RR")){
      cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Overall]")),
          "Mantel-Haenszel method used to pool dichotomous outcome data. \n")
    }

    cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Overall]")),
        "Fixed ('common') effect model assumed. \n")

    if ("outliers" %in% x$which.run){
      cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Outliers removed]")),
          "ES removed as outliers if the CI did not overlap with pooled effect CI. \n")}

    if ("influence" %in% x$which.run){
      cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Influence Analysis]")),
          "Influential cases determined using diagnostics of Viechtbauer and Cheung (2010). \n")}

    if ("combined" %in% x$which.run){
      cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Combined]")),
          "'Combined' analysis: ES combined on", x$model.combined$which.combine,
          "level assuming a correlation of rho =", paste0(x$model.combined$rho, ". \n"))}

    if ("threelevel" %in% x$which.run || "threelevel.che" %in% x$which.run){
      cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Three-Level Model]")),
          "Three-level model estimated via restricted maximum likelihood, using the 'rma.mv'
  function in {metafor}. \n")
      if (x$model.threelevel$test == "t" | x$model.threelevel.che$test[1] == "t"){
        cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Three-Level Model]")),
            "Test statistics and CIs of the three-level model calculated based on
  a t-distribution. \n")
      } else {
        cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Three-Level Model]")),
            "Test statistics and CIs of the three-level model calculated based on
  a Wald-type normal approximation. \n")
      }
      if (x$use.rve[1] == TRUE){
        cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Robust Variance Estimation]")),
            "Robust variance estimation was used to guard the three-level model(s)
  from misspecification. This was done using functions of the {clubSandwich} package.\n")}
      }


    cat("\n")
    cat("\n")
    cat("Cite the following packages/methods: ---------------------------------------\n")
    cat("\n")
    cat(" -", crayon::bold(crayon::blue("{meta}:")), "Balduzzi S, R\u00FCcker G, Schwarzer G (2019),
          How to perform a meta-analysis with R: a practical tutorial,
          Evidence-Based Mental Health; 22: 153-160. \n")
    cat(" -", crayon::bold(crayon::blue("{metafor}:")), "Viechtbauer, W. (2010). Conducting meta-analyses in R
          with the metafor package. Journal of Statistical Software, 36(3), 1-48.
          https://doi.org/10.18637/jss.v036.i03. \n")
    cat(" -", crayon::bold(crayon::blue("{dmetar}:")), "Harrer, M., Cuijpers, P., Furukawa, T.A., & Ebert, D.D. (2021).
          Doing Meta-Analysis with R: A Hands-On Guide. Boca Raton, FL and London:
          Chapman & Hall/CRC Press. ISBN 978-0-367-61007-4. \n")
    cat(" -", crayon::bold(crayon::blue("R:")), "R Core Team (2021). R: A language and environment for statistical computing.
          R Foundation for Statistical Computing, Vienna, Austria.
          URL https://www.R-project.org/.\n")
    if ("threelevel" %in% x$which.run || "threelevel.che" %in% x$which.run){
      if (x$use.rve[1] == TRUE){
        cat(" -", crayon::bold(crayon::blue("{clubSandwich}:")), "Pustejovsky, J. (2021). clubSandwich:
          Cluster-Robust (Sandwich) Variance Estimators with Small-Sample Corrections.
          R package version 0.5.3. https://CRAN.R-project.org/package=clubSandwich \n")
        }
      }
    if (class(x)[2] == "correctPublicationBias"){
      cat(" -", crayon::bold(crayon::blue("{metasens}:")), "Schwarzer G, Carpenter J, R\u00FCcker G (2022). 
          metasens: Statistical Methods for Sensitivity Analysis in Meta-Analysis. 
          R package version 1.0-1, https://CRAN.R-project.org/package=metasens. \n")}
    if (isTRUE(x$.raw.bin.es) & identical(x$.type.es, "RR")){
      cat(" -", crayon::bold(crayon::blue("Mantel-Haenszel method:")), "Greenland S & Robins JM (1985). Estimation of a 
          common effect parameter from sparse follow-up data. Biometrics, 41, 55-68 \n")}
    if ("influence" %in% x$which.run){
      cat(" -", crayon::bold(crayon::blue("Influential cases:")), "Viechtbauer, W., & Cheung, M. W.-L. (2010). Outlier and influence
          diagnostics for meta-analysis. Research Synthesis Methods, 1(2), 112-125.
          https://doi.org/10.1002/jrsm.11 \n")}
    if ("threelevel.che" %in% x$which.run){
      cat(" -", crayon::bold(crayon::blue("Three-level CHE model:")), "Pustejovsky, J.E., Tipton, E. Meta-analysis with Robust
          Variance Estimation: Expanding the Range of Working Models. Prevention
          Science (2021). https://doi.org/10.1007/s11121-021-01246-3 \n")}
    if ("threelevel" %in% x$which.run || "threelevel.che" %in% x$which.run){
      if (x$use.rve[1] == TRUE){
        cat(" -", crayon::bold(crayon::blue("Robust variance estimation:")), "Pustejovsky, J.E., Tipton, E. Meta-analysis with Robust
          Variance Estimation: Expanding the Range of Working Models. Prevention
          Science (2021). https://doi.org/10.1007/s11121-021-01246-3 \n")} 
    }
    if (class(x)[2] == "correctPublicationBias"){
      cat(" -", crayon::bold(crayon::blue("Trim and fill method:")), "Duval S & Tweedie R (2000a): A nonparametric 'Trim and Fill' 
          method of accounting for publication bias in meta-analysis. Journal 
          of the American Statistical Association, 95, 89-98 \n")
      cat(" -", crayon::bold(crayon::blue("Limit meta-analysis:")), "R\u00FCcker G, Schwarzer G, Carpenter JR, Binder H, Schumacher M (2011): 
          Treatment-effect estimates adjusted for small-study effects via a 
          limit meta-analysis. Biostatistics, 12, 122-42 \n")
      cat(" -", crayon::bold(crayon::blue("Selection model:")), "Vevea, J. L., & Hedges, L. V. (1995). A general linear model 
          for estimating effect size in the presence of publication bias. 
          Psychometrika, 60(3), 419-435. https://doi.org/10.1007/BF02294384. \n")
    }


  } else {

    tau.methods = data.frame(method = c("REML", "DL", "PM", "ML", "HS", "SJ", "HE", "EB"),
                             name = c("restricted maximum-likelihood", "DerSimonian-Laird",
                                      "Paule-Mandel", "maximum likelihood", "Hunter and Schmidt",
                                      "Sidik-Jonkman", "Hedges", "empirical Bayes"),
                             citation =
                               c("Viechtbauer W. (2005): Bias and efficiency of meta-analytic variance
          estimators in the random-effects model.Journal of Educational and
          Behavioral Statistics, 30, 261-93",
          "DerSimonian R. & Laird N. (1986): Meta-analysis in clinical trials.
          Controlled Clinical Trials, 7, 177-88",
          "Paule RC & Mandel J (1982): Consensus values and weighting factors.
          Journal of Research of the National Bureau of Standards, 87, 377-85",
          "Viechtbauer W. (2005): Bias and efficiency of meta-analytic
          variance estimators in the random-effects model. Journal of Educational
          and Behavioral Statistics, 30, 261-93",
          "Hunter JE & Schmidt FL (2015): Methods of Meta-Analysis: Correcting
          Error and Bias in Research Findings (3rd edition). Thousand Oaks, CA: Sage",
          "Sidik K & Jonkman JN (2005): Simple heterogeneity variance estimation
          for meta-analysis. Journal of the Royal Statistical Society:
          Series C (Applied Statistics), 54, 367-84",
          "Hedges LV & Olkin I (1985): Statistical methods for meta-analysis.
          San Diego, CA: Academic Press",
          "Morris CN (1983): Parametric empirical Bayes inference: Theory
          and applications (with discussion). Journal of the American
          Statistical Association 78, 47-65"))


    cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Overall]")),
        "Random effects model assumed. \n")
    
    if (isTRUE(x$.raw.bin.es) & identical(x$.type.es, "RR")){
      cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Overall]")),
          "Mantel-Haenszel method used to pool dichotomous outcome data. \n")
    }

    cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Overall]")),
        "Heterogeneity variance (tau2) calculated using",
        tau.methods[tau.methods$method == x$model.overall$method.tau, "name"],
        "estimator. \n")

    if (x$model.overall$hakn == TRUE){
      cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Overall]")),
          "Test statistic and CI of the pooled effect calculated using the Knapp-Hartung adjustment. \n")
    } else {
      cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Overall]")),
          "Test statistic and CI of the pooled effect calculated based on
  a Wald-type normal approximation. \n")
    }

    if ("outliers" %in% x$which.run){
      cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Outliers removed]")),
          "ES removed as outliers if the CI did not overlap with pooled effect CI. \n")}

    if ("influence" %in% x$which.run){
      cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Influence Analysis]")),
          "Influential cases determined using diagnostics of Viechtbauer and Cheung (2010). \n")}

    if ("combined" %in% x$which.run){
      cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Combined]")),
          "'Combined' analysis: ES combined on", x$model.combined$which.combine,
          "level assuming a correlation of rho =", paste0(x$model.combined$rho, ". \n"))}

    if ("threelevel" %in% x$which.run || "threelevel.che" %in% x$which.run){
      cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Three-Level Model]")),
          "Three-level model estimated via restricted maximum likelihood, using the 'rma.mv'
  function in {metafor}. \n")

      if (x$model.threelevel$test[1] == "t" | x$model.threelevel.che$test[1] == "t"){
        cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Three-Level Model]")),
            "Test statistics and CIs of the three-level model calculated based on
  a t-distribution. \n")
      } else {
        cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Three-Level Model]")),
            "Test statistics and CIs of the three-level model calculated based on
  a Wald-type normal approximation. \n")
      }
      if (x$use.rve[1] == TRUE){
        cat(crayon::green("\u2713"), crayon::bold(crayon::magenta("[Robust Variance Estimation]")),
            "Robust variance estimation was used to guard the three-level model(s)
  from misspecification. This was done using functions of the {clubSandwich} package.\n")}  
    }



    cat("\n")
    cat("\n")
    cat("Cite the following packages/methods: ---------------------------------------\n")
    cat("\n")
    cat(" -", crayon::bold(crayon::blue("{meta}:")), "Balduzzi S, R\u00FCcker G, Schwarzer G (2019),
          How to perform a meta-analysis with R: a practical tutorial,
          Evidence-Based Mental Health; 22: 153-160. \n")
    cat(" -", crayon::bold(crayon::blue("{metafor}:")), "Viechtbauer, W. (2010). Conducting meta-analyses in R
          with the metafor package. Journal of Statistical Software, 36(3), 1-48.
          https://doi.org/10.18637/jss.v036.i03. \n")
    cat(" -", crayon::bold(crayon::blue("{dmetar}:")), "Harrer, M., Cuijpers, P., Furukawa, T.A., & Ebert, D.D. (2021).
          Doing Meta-Analysis with R: A Hands-On Guide. Boca Raton, FL and London:
          Chapman & Hall/CRC Press. ISBN 978-0-367-61007-4. \n")
    cat(" -", crayon::bold(crayon::blue("R:")), "R Core Team (2021). R: A language and environment for statistical computing.
          R Foundation for Statistical Computing, Vienna, Austria.
          URL https://www.R-project.org/.\n")
    if ("threelevel" %in% x$which.run || "threelevel.che" %in% x$which.run){
      if (x$use.rve[1] == TRUE){
        cat(" -", crayon::bold(crayon::blue("{clubSandwich}:")), "Pustejovsky, J. (2021). clubSandwich:
          Cluster-Robust (Sandwich) Variance Estimators with Small-Sample Corrections.
          R package version 0.5.3. https://CRAN.R-project.org/package=clubSandwich \n")} 
    }
    if (class(x)[2] == "correctPublicationBias"){
      cat(" -", crayon::bold(crayon::blue("{metasens}:")), "Schwarzer G, Carpenter J, R\u00FCcker G (2022). 
          metasens: Statistical Methods for Sensitivity Analysis in Meta-Analysis. 
          R package version 1.0-1, https://CRAN.R-project.org/package=metasens. \n")}
    if (isTRUE(x$.raw.bin.es) & identical(x$.type.es, "RR")){
      cat(" -", crayon::bold(crayon::blue("Mantel-Haenszel method:")), "Greenland S & Robins JM (1985). Estimation of a 
          common effect parameter from sparse follow-up data. Biometrics, 41, 55-68 \n")}
    if ("influence" %in% x$which.run){
      cat(" -", crayon::bold(crayon::blue("Influential cases:")), "Viechtbauer, W., & Cheung, M. W.-L. (2010). Outlier and influence
          diagnostics for meta-analysis. Research Synthesis Methods, 1(2), 112-125.
          https://doi.org/10.1002/jrsm.11 \n")}
    cat(" -", crayon::bold(crayon::blue("tau2 estimator:")),
        tau.methods[tau.methods$method == x$model.overall$method.tau, "citation"], "\n")
    if (x$model.overall$hakn == TRUE){
      cat(" -", crayon::bold(crayon::blue("Knapp-Hartung:")), "Hartung J, Knapp G (2001a): On tests of the overall treatment
          effect in meta-analysis with normally distributed responses.
          Statistics in Medicine, 20, 1771-82 \n")}
    if ("threelevel.che" %in% x$which.run){
      cat(" -", crayon::bold(crayon::blue("Three-level CHE model:")), "Pustejovsky, J.E., Tipton, E. Meta-analysis with Robust
          Variance Estimation: Expanding the Range of Working Models. Prevention
          Science (2021). https://doi.org/10.1007/s11121-021-01246-3 \n")}
    if ("threelevel" %in% x$which.run || "threelevel.che" %in% x$which.run) {
      if (x$use.rve[1] == TRUE){
        cat(" -", crayon::bold(crayon::blue("Robust variance estimation:")), "Pustejovsky, J.E., Tipton, E. Meta-analysis with Robust
          Variance Estimation: Expanding the Range of Working Models. Prevention
          Science (2021). https://doi.org/10.1007/s11121-021-01246-3 \n")} 
    }
    if (class(x)[2] == "correctPublicationBias"){
      cat(" -", crayon::bold(crayon::blue("Trim and fill method:")), "Duval S & Tweedie R (2000a): A nonparametric 'Trim and Fill' 
          method of accounting for publication bias in meta-analysis. Journal 
          of the American Statistical Association, 95, 89-98 \n")
      cat(" -", crayon::bold(crayon::blue("Limit meta-analysis:")), "R\u00FCcker G, Schwarzer G, Carpenter JR, Binder H, Schumacher M (2011): 
          Treatment-effect estimates adjusted for small-study effects via a 
          limit meta-analysis. Biostatistics, 12, 122-42 \n")
      cat(" -", crayon::bold(crayon::blue("Selection model:")), "Vevea, J. L., & Hedges, L. V. (1995). A general linear model 
          for estimating effect size in the presence of publication bias. 
          Psychometrika, 60(3), 419-435. https://doi.org/10.1007/BF02294384. \n")
    }
    
  }

  if (forest[1]){

    models = list("overall" = "Overall", 
                  "lowest.highest" = c("One ES/study (lowest)",
                                       "One ES/study (highest)"), 
                  "outliers" = "Outliers removed",
                  "influence" = "Influence Analysis", 
                  "rob" = "rob", 
                  "combined" = "Combined", 
                  "threelevel" = "Three-Level Model",
                  "threelevel.che" = "Three-Level Model (CHE)")
    
    run.models = unlist(models[x$which.run])
    if ("rob" %in% run.models){
      run.models[["rob"]] = 
        rownames(x$summary)[grepl("Only", rownames(x$summary))]
    }

    x$summary = x$summary[run.models,]

    if (x$.type.es == "RR"){
      stringr::str_replace_all(x$summary$rr.ci, ";|\\]|\\[", "") %>%
        strsplit(" ") %>% 
        purrr::map(~as.numeric(.)) %>% do.call(rbind,.) %>%
        {colnames(.) = c("lower", "upper");.} %>%
        cbind(model = rownames(x$summary), rr = x$summary$rr,
              i2 = round(x$summary$i2,1) %>% format(1), .) %>%
        data.frame() %>%
        dplyr::mutate(rr = as.numeric(rr) %>% log(),
                      lower = as.numeric(lower) %>% log(),
                      upper = as.numeric(upper) %>% log()) %>%
        meta::metagen(TE = rr, 
                      lower = lower - 1e-50, 
                      upper = upper + 1e-50, 
                      studlab = model, sm = "RR",
                      data = .) %>%
        meta::forest.meta(
          col.square = "lightblue",
          rightcols = FALSE,
          overall.hetstat = FALSE,
          weight.study = "same",
          test.overall = FALSE, overall = FALSE,
          leftcols = c("studlab", "effect", "ci", "i2"),
          leftlabs = c(expression(bold(Model)), 
                       expression(bold(RR)),
                       expression(bold(CI)),
                       expression(bold(italic(I)^2)))) %>%
        suppressWarnings()
    } else {
      stringr::str_replace_all(x$summary$g.ci, ";|\\]|\\[", "") %>%
        strsplit(" ") %>% purrr::map(~as.numeric(.)) %>% do.call(rbind,.) %>%
        {colnames(.) = c("lower", "upper");.} %>%
        cbind(model = rownames(x$summary), g = x$summary$g,
              i2 = round(x$summary$i2,1) %>% format(1), .) %>%
        data.frame() %>%
        dplyr::mutate(g = as.numeric(g),
                      lower = as.numeric(lower) - 1e-50,
                      upper = as.numeric(upper) + 1e-50) %>%
        meta::metagen(TE = g, lower = lower, upper = upper, studlab = model,
                      data = .) %>%
        meta::forest.meta(col.square = "lightblue",
                          rightcols = FALSE,
                          overall.hetstat = FALSE,
                          weight.study = "same",
                          test.overall = FALSE, overall = FALSE,
                          leftcols = c("studlab", "TE", "ci", "i2"),
                          leftlabs = c(expression(bold(Model)), expression(bold(g)),
                                       expression(bold(CI)),
                                       expression(bold(italic(I)^2)))) %>%
        suppressWarnings()
    }
  }

}



#' Best Linear Unbiased Predictions (BLUPs) for 'runMetaAnalysis' models.
#' 
#' Generates empirical Bayes (EB) estimates, also known as best linear unbiased 
#' predictions (BLUPs), by merging the fitted values obtained from fixed effects 
#' and estimated contributions of random effects. These estimates represent 
#' the study-specific true effect sizes or outcomes and are accompanied by 
#' standard errors and prediction interval bounds. Uses the [metafor::blup.rma.uni()] function
#' internally.
#' 
#' @param x An object of class \code{runMetaAnalysis}.
#' @param which Model for which estimates should be printed. Can be one of \code{"overall"},
#' \code{"combined"}, \code{"lowest.highest"}, \code{"outliers"},
#' \code{"influence"}, \code{"threelevel"}, or \code{"threelevel.che"}.
#' @param ... Additional arguments.
#' 
#' @importFrom crayon green yellow
#' @importFrom dplyr select ends_with filter group_map group_by
#' @importFrom stringr str_replace_all str_remove_all
#' @importFrom metafor rma.uni blup.rma.uni
#' @importFrom clubSandwich conf_int
#' @export
#' @method eb runMetaAnalysis
eb.runMetaAnalysis = function(x, which = NULL, ...){
  
  model = x
  
  # Check overall compatibility
  if (is.null(which[1])){
    which.run = model$which.run[1]
  } else {
    if (!which %in% model$which.run){
      stop("argument 'which' must be one of ",
           paste(model$which.run, collapse=", "), ".")
    } else {
      which.run = which
    }
  }
  
  
  if (which.run %in% c("trimfill", "limitmeta", "selection")){
    stop("Empirical Bayes estimates are not supported for '",
         which.run, "' models.")
  }
  
  
  # Switch to combined if three-level was used
  threeLevel = FALSE
  M.3l = NULL
  if (which.run %in% c("threelevel", "threelevel.che")){
    threeLevel = TRUE
    whichThreeLevel = which.run
    which.run = "combined"
    if (!which.run %in% model$which.run){
      stop("Model 'combined' not found. To generate an empirical",
           " Bayes estimates for a three-level model,",
           " make sure that 'which.run' includes 'combined'.")
    }
  }
  
  # Transform model to metafor
  models = list("overall" = "model.overall",
                "lowest.highest" = c("model.lowest", "model.highest"),
                "outliers" = "model.outliers",
                "influence" = "model.influence",
                "rob" = "model.rob",
                "combined" = "model.combined",
                "threelevel" = "model.threelevel",
                "che" = "model.threelevel.che",
                "threelevel.che" = "model.threelevel.che")
  
  if (identical(which.run[1], "lowest.highest")){
    
    # Generate "lowest" plot
    M = model[models[[which.run[1]]]][[1]]
    if (!is.null(M$exclude)){
      M$data = M$data[!M$exclude,]
    }
    dat = escalc(yi=.TE, sei=.seTE, data = M$data)
    res = metafor::rma(yi = .TE, sei = .seTE, dat = dat, slab = .studlab,
                       method = M$method.tau, test = ifelse(M$hakn, "knha", "z"))
    lowest = metafor::blup(res)
    
    # Generate "highest" plot
    M = model[models[[which.run[1]]]][[2]]
    if (!is.null(M$exclude)){
      M$data = M$data[!M$exclude,]
    }
    dat = escalc(yi=.TE, sei=.seTE, data = M$data)
    res = metafor::rma(yi = .TE, sei = .seTE, dat = dat, slab = .studlab,
                       method = M$method.tau, test = ifelse(M$hakn, "knha", "z"))
    highest = metafor::blup(res)
    
    message("- ", crayon::green("[OK] "), 
            "Calculated EB estimates/BLUPs ('", which.run,"' model).")
    
    return(list(lowest = lowest, 
                highest = highest))
    
    
  } else {
    
    M = model[models[[which.run[1]]]][[1]]
    if (!is.null(M$exclude)){
      M$data = M$data[!M$exclude,]
    }
    if (threeLevel) {
      warning("Empirical Bayes estimates are not available for three-level models.",
              " Values were calculated for the 'combined' model.")
    }
    dat = escalc(yi=.TE, sei=.seTE, data = M$data)
    res = metafor::rma(yi = .TE, sei = .seTE, dat = dat, slab = .studlab,
                       method = M$method.tau, test = ifelse(M$hakn, "knha", "z"))
    
    message("- ", crayon::green("[OK] "), 
            "Calculated EB estimates/BLUPs ('", which.run,"' model).")
    return(metafor::blup(res))
  }
}



#' Best Linear Unbiased Predictions (BLUPs) for 'runMetaAnalysis' models.
#' 
#' Generates empirical Bayes (EB) estimates, also known as best linear unbiased 
#' predictions (BLUPs), by merging the fitted values obtained from fixed effects 
#' and estimated contributions of random effects. These estimates represent 
#' the study-specific true effect sizes or outcomes and are accompanied by 
#' standard errors and prediction interval bounds. Uses the [metafor::blup.rma.uni()] function
#' internally.
#' 
#' @param x An object of class \code{runMetaAnalysis}.
#' @param which Model for which estimates should be printed. Can be one of \code{"overall"},
#' \code{"combined"}, \code{"lowest.highest"}, \code{"outliers"},
#' \code{"influence"}, \code{"threelevel"}, or \code{"threelevel.che"}.
#' @param ... Additional arguments.
#' 
#' @importFrom crayon green yellow
#' @importFrom dplyr select ends_with filter group_map group_by
#' @importFrom stringr str_replace_all str_remove_all
#' @importFrom metafor rma.uni blup.rma.uni
#' @importFrom clubSandwich conf_int
#' @export
#' @method blup runMetaAnalysis
blup.runMetaAnalysis = function(x, which = NULL, ...){
  
  model = x
  
  # Check overall compatibility
  if (is.null(which[1])){
    which.run = model$which.run[1]
  } else {
    if (!which %in% model$which.run){
      stop("argument 'which' must be one of ",
           paste(model$which.run, collapse=", "), ".")
    } else {
      which.run = which
    }
  }
  
  
  if (which.run %in% c("trimfill", "limitmeta", "selection")){
    stop("Empirical Bayes estimates are not supported for '",
         which.run, "' models.")
  }
  
  
  # Switch to combined if three-level was used
  threeLevel = FALSE
  M.3l = NULL
  if (which.run %in% c("threelevel", "threelevel.che")){
    threeLevel = TRUE
    whichThreeLevel = which.run
    which.run = "combined"
    if (!which.run %in% model$which.run){
      stop("Model 'combined' not found. To generate an empirical",
           " Bayes estimates for a three-level model,",
           " make sure that 'which.run' includes 'combined'.")
    }
  }
  
  # Transform model to metafor
  models = list("overall" = "model.overall",
                "lowest.highest" = c("model.lowest", "model.highest"),
                "outliers" = "model.outliers",
                "influence" = "model.influence",
                "rob" = "model.rob",
                "combined" = "model.combined",
                "threelevel" = "model.threelevel",
                "che" = "model.threelevel.che",
                "threelevel.che" = "model.threelevel.che")
  
  if (identical(which.run[1], "lowest.highest")){
    
    # Generate "lowest" plot
    M = model[models[[which.run[1]]]][[1]]
    if (!is.null(M$exclude)){
      M$data = M$data[!M$exclude,]
    }
    dat = escalc(yi=.TE, sei=.seTE, data = M$data)
    res = metafor::rma(yi = .TE, sei = .seTE, dat = dat, slab = .studlab,
                       method = M$method.tau, test = ifelse(M$hakn, "knha", "z"))
    lowest = metafor::blup(res)
    
    # Generate "highest" plot
    M = model[models[[which.run[1]]]][[2]]
    if (!is.null(M$exclude)){
      M$data = M$data[!M$exclude,]
    }
    dat = escalc(yi=.TE, sei=.seTE, data = M$data)
    res = metafor::rma(yi = .TE, sei = .seTE, dat = dat, slab = .studlab,
                       method = M$method.tau, test = ifelse(M$hakn, "knha", "z"))
    highest = metafor::blup(res)
    
    message("- ", crayon::green("[OK] "), 
            "Calculated EB estimates/BLUPs ('", which.run,"' model).")
    
    return(list(lowest = lowest, 
                highest = highest))
    
    
  } else {
    
    M = model[models[[which.run[1]]]][[1]]
    if (!is.null(M$exclude)){
      M$data = M$data[!M$exclude,]
    }
    if (threeLevel) {
      warning("Empirical Bayes estimates are not available for three-level models.",
              " Values were calculated for the 'combined' model.")
    }
    dat = escalc(yi=.TE, sei=.seTE, data = M$data)
    res = metafor::rma(yi = .TE, sei = .seTE, dat = dat, slab = .studlab,
                       method = M$method.tau, test = ifelse(M$hakn, "knha", "z"))
    
    message("- ", crayon::green("[OK] "), 
            "Calculated EB estimates/BLUPs ('", which.run,"' model).")
    return(metafor::blup(res))
  }
}



#' Profile Likelihood Plots for 'runMetaAnalysis' models.
#' 
#' Profiles the restricted log-likelihood of `threelevel` and `threelevel.che`
#' models, using the [metafor::profile.rma()] function. This functionality
#' can be used to check if the two heterogeneity variances (\mjeqn{\tau^2}{\tau^2} within 
#' and between studies) were identifiable and correctly estimated.
#' 
#' @param fitted An object of class \code{runMetaAnalysis}.
#' @param which Model for which estimates should be printed. Can be one of \code{"threelevel"} 
#' or \code{"threelevel.che"}.
#' @param ... Additional arguments.
#' 
#' @importFrom crayon green yellow
#' @importFrom dplyr select ends_with filter group_map group_by
#' @importFrom stringr str_replace_all str_remove_all
#' @importFrom metafor rma.uni profile.rma.mv
#' @export
#' @method profile runMetaAnalysis

profile.runMetaAnalysis = function(fitted, which = NULL, ...){
  
  x = fitted
  models = list("threelevel" = "model.threelevel",
                "threelevel.che" = "model.threelevel.che")
  
  if(is.null(which)){
    if (sum(x$which.run %in% c("threelevel", "threelevel.che")) == 0){
      stop("Either 'threelevel' or 'threelevel.che' needed in the 'runMetaAnalysis' model.")
    }
    which.run = x$which.run[x$which.run %in% c("threelevel", "threelevel.che")][1]
    message("- ", crayon::green("[OK] "), 
            "Generating profile likelihood plot ('", which.run,"' model).")
    metafor::profile.rma.mv(x[models[[which.run]]][[1]])
    
  } else {
    if (!which[1] %in% x$which.run){
      stop("Either 'threelevel' or 'threelevel.che' needed in the 'runMetaAnalysis' model.")
    }
    which.run = which
    message("- ", crayon::green("[OK] "), 
            "Generating profile likelihood plot ('", which.run,"' model).")
    metafor::profile.rma.mv(x[models[[which.run]]][[1]])
  }
}



