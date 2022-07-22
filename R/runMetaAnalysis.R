#' Run different types of meta-analyses
#'
#' This wrapper function allows to simultaneously pool effect sizes using
#' different meta-analytic approaches.
#'
#' @usage runMetaAnalysis(data,
#'                 which.run = c("overall", "combined",
#'                               "lowest.highest", "outliers",
#'                               "influence", "rob", "threelevel",
#'                               "threelevel.che"),
#'                 es.measure = c("g", "RR"),
#'                 es.type = c("precalculated", "raw"),
#'                 es.var = ifelse(identical(es.measure[1], "g"), 
#'                                 ".g", ".log_rr"),
#'                 se.var = ifelse(identical(es.measure[1], "g"), 
#'                                 ".g_se", ".log_rr_se"),
#'                 es.binary.raw.vars = 
#'                   c(".event_arm1", ".event_arm2",
#'                     ".totaln_arm1", ".totaln_arm2"),
#'                 method.tau = "REML",
#'                 hakn = TRUE,
#'                 study.var = "study",
#'                 extra.grouping.var = NULL,
#'                 arm.var.1 = "condition_arm1",
#'                 arm.var.2 = "condition_arm2",
#'                 measure.var = "instrument",
#'                 low.rob.filter = "rob > 2",
#'                 method.tau.ci = "Q-Profile",
#'                 round.digits = 2,
#'                 which.outliers = c("overall", "combined"),
#'                 which.influence = c("overall", "combined"),
#'                 which.rob = c("overall", "combined"),
#'                 nnt.cer = 0.2,
#'                 rho.within.study = 0.5,
#'                 use.rve = TRUE,
#'                 html = TRUE,
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
#' @param hakn \code{logical}. Should the Knapp-Hartung adjustment for effect size significance tests be used? Default is \code{TRUE}.
#' @param study.var \code{character}. The name of the variable in \code{data} in which the study IDs are stored.
#' @param extra.grouping.var \code{character}. Additional grouping variable within studies to be used for the \code{"combined"}
#' analysis. This is useful, for example, to \emph{not} pool the effects of different treatment arms within multi-arm trials.
#' \code{NULL} by default.
#' @param arm.var.1 \code{character}. The name of the variable in \code{data} in which the condition (e.g. "guided iCBT")
#' of the \emph{first} arm within a comparison are stored.
#' @param arm.var.2 \code{character}. The name of the variable in \code{data} in which the condition (e.g. "wlc")
#' of the \emph{second} arm within a comparison are stored.
#' @param measure.var \code{character}. The name of the variable in \code{data} 
#' in which the instrument used for the comparison is stored.
#' @param low.rob.filter \code{character}. A filtering statement by which to 
#' include studies for the "low RoB only" analysis. Please note that the name of 
#' the variable must be included as a column in \code{data}.
#' @param method.tau.ci \code{character}. A character string indicating which method is used to estimate the
#' confidence interval of tau/tau-squared. Either \code{"Q-Profile"} (default and recommended),
#' \code{"BJ"}, \code{"J"}, or \code{"PL"} can be abbreviated. See \code{\link[meta]{metagen}} for details.
#' @param round.digits \code{numeric}. Number of digits to round the (presented) results by. Default is \code{2}.
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
#' the variance-covariance matrices needed for the conditional and hierarchical effects three-level model. Default is \code{0.5}.
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
#' # For the combined analysis, we provide an extra variable to
#' # 'extra.grouping.var' so that different treatment arm 
#' # comparisons in multiarm trials are NOT pooled
#' data %>% 
#'   runMetaAnalysis(extra.grouping.var = "condition_arm1") %>% 
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
#' For more details see the [Get Started](https://tools.metapsy.org/articles/metapsytools) vignette.
#'
#' @import dplyr
#' @importFrom crayon green yellow cyan bold
#' @importFrom scales pvalue
#' @importFrom purrr map
#' @importFrom meta update.meta metagen
#' @importFrom metafor escalc aggregate.escalc rma.mv
#' @importFrom clubSandwich coef_test conf_int
#' @importFrom stats dffits model.matrix rnorm rstudent
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
                           hakn = TRUE,
                           study.var = "study",
                           extra.grouping.var = NULL,
                           arm.var.1 = "condition_arm1",
                           arm.var.2 = "condition_arm2",
                           measure.var = "instrument",
                           low.rob.filter = "rob > 2",
                           method.tau.ci = "Q-Profile",
                           round.digits = 2,
                           which.outliers = c("overall", "combined"),
                           which.influence = c("overall", "combined"),
                           which.rob = c("overall", "combined"),
                           nnt.cer = 0.2,
                           rho.within.study = 0.5,
                           use.rve = TRUE,
                           html = TRUE,
                           ...){

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

  warn.end = FALSE
  
  # Send message (beginning of analyses)
  message(crayon::cyan(crayon::bold("- Running meta-analyses...")))
  message("- ", crayon::green("[OK] "), 
          "Using ", ifelse(identical(.type.es, "RR"), 
                    ifelse(identical(es.type[1], "precalculated"), 
                      crayon::bold(crayon::magenta(
                        "risk ratio (pre-calculated)")), 
                      crayon::bold(crayon::magenta(
                        "risk ratio (raw event data)"))), 
                    crayon::bold(crayon::magenta("Hedges' g"))),
          " as effect size metric... ")

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  1. Overall Model                                                         #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  data.original = data
  round.digits = abs(round.digits)

  # Create comparison variable
  paste0(data[[study.var]], " (",
         data[[arm.var.1]], " vs. ",
         data[[arm.var.2]], "; ",
         data[[measure.var]], ")") -> data$comparison

  paste0(data[[arm.var.1]], " vs. ",
         data[[arm.var.2]]) -> data$comparison.only

  data[[measure.var]] -> data$instrument
  
  # Define multi.study
  multi.study = names(table(data[[study.var]])
                      [table(data[[study.var]]) > 1])

  # Select methods
  method.tau.ci = ifelse(identical(method.tau.ci,"Q-Profile"), 
                         "QP", method.tau.ci)
  method.tau.meta = ifelse(identical(method.tau[1], "FE"),
                           "REML", method.tau[1])
  
  if ("overall" %in% which.run){
    message("- ", crayon::green("[OK] "), 
            "Calculating overall effect size... ", 
            appendLF = FALSE);
  }

  if (!isTRUE(.raw.bin.es)){
    mGeneral = meta::metagen(TE = data[[es.var[1]]],
                seTE = data[[se.var[1]]],
                studlab = data[[study.var]],
                data = data,
                sm = .type.es,
                hakn = hakn,
                method.tau = method.tau.meta,
                method.tau.ci = method.tau.ci,
                prediction = TRUE,
                fixed = ifelse(method.tau == "FE", TRUE, FALSE),
                random = ifelse(method.tau == "FE", FALSE, TRUE),
                ...)
  } else {
    mGeneral = meta::metabin(event.e = data[[es.binary.raw.vars[1]]],
                n.e = data[[es.binary.raw.vars[3]]],
                event.c = data[[es.binary.raw.vars[2]]],
                n.c = data[[es.binary.raw.vars[4]]],
                studlab = data[[study.var]],
                data = data,
                sm = "RR",
                hakn = hakn,
                method.tau = method.tau.meta,
                method.tau.ci = method.tau.ci,
                prediction = TRUE,
                fixed = ifelse(method.tau == "FE", TRUE, FALSE),
                random = ifelse(method.tau == "FE", FALSE, TRUE),
                ...)
    with(meta::update.meta(mGeneral, sm="RD"), {
      ifelse(isTRUE(fixed) & isTRUE(random),
              abs(TE.fixed)^-1, abs(TE.random)^-1)
    }) -> nnt.raw.bin.es
  }


  mGeneralRes = with(mGeneral, {

    data.frame(
      k = k,
      g = ifelse(
        isTRUE(fixed) & isTRUE(random), 
        TE.fixed, TE.random) %>% 
        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
        round(round.digits),
      g.ci = paste0("[", 
                    ifelse(isTRUE(fixed) & isTRUE(random),
                    lower.fixed, lower.random) %>% 
                      ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                      round(round.digits), "; ",
                    ifelse(isTRUE(fixed) & isTRUE(random),
                    upper.fixed, upper.random) %>% 
                      ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                      round(round.digits), "]"),
      p = ifelse(isTRUE(fixed) & isTRUE(random), 
                 pval.fixed, pval.random) %>% 
                 scales::pvalue(),
      i2 = round(I2*100, 2),
      i2.ci = paste0("[", round(lower.I2*100, 2), "; ", 
                     round(upper.I2*100, 2), "]"),
      prediction.ci = paste0("[", 
                             round(
                               ifelse(identical(.type.es, "RR"), 
                                      exp(lower.predict), lower.predict), 
                               round.digits), "; ",
                             round(
                               ifelse(identical(.type.es, "RR"), 
                                      exp(upper.predict), upper.predict), 
                               round.digits), "]"),
      nnt = metapsyNNT(
        ifelse(isTRUE(fixed) & isTRUE(random), 
               TE.fixed, TE.random), nnt.cer) %>%
        ifelse(identical(.type.es, "RR"), NA, .) %>% 
        ifelse(isTRUE(.raw.bin.es), 
               nnt.raw.bin.es, .) %>% 
        round(round.digits) %>% abs(),
      excluded = "none"
    )
  })
  rownames(mGeneralRes) = "Overall"

  if ("overall" %in% which.run){
    message(crayon::green("DONE"))
  }


  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  2. Lowest Only                                                           #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  if ("lowest.highest" %in% which.run){
    
    message("- ", crayon::green("[OK] "), 
              "Calculating effect size using only lowest effect... ",
              appendLF = FALSE)
    
    multi.study = names(table(data[[study.var]])
                        [table(data[[study.var]]) > 1])
  
    if (length(multi.study) > 0){
      if (length(multi.study) > 0){
        data$.TE = mGeneral$TE
        data %>%
          split(.[[study.var]]) %>%
          purrr::map(function(x){
            tiebreaker = rnorm(nrow(x), sd=1e-10)
            x$.TE + tiebreaker == min(x$.TE + tiebreaker)}) %>%
          do.call(c,.) -> lowest
        data$.TE = NULL
      } else {
        lowest = NULL
      }
  
      mLowest = meta::update.meta(mGeneral, exclude = !lowest, 
                                  id = NULL)
      with(update.meta(mLowest, sm="RD"), {
        ifelse(isTRUE(fixed) & isTRUE(random),
               abs(TE.fixed)^-1, abs(TE.random)^-1)
      }) -> nnt.raw.bin.es
      
      mLowestRes = with(mLowest,{
        data.frame(
          k = k,
          g = ifelse(
            isTRUE(fixed) & isTRUE(random), 
            TE.fixed, TE.random) %>% 
            ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
            round(round.digits),
          g.ci = paste0("[", 
                        ifelse(isTRUE(fixed) & isTRUE(random),
                               lower.fixed, lower.random) %>% 
                          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                          round(round.digits), "; ",
                        ifelse(isTRUE(fixed) & isTRUE(random),
                               upper.fixed, upper.random) %>% 
                          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                          round(round.digits), "]"),
          p = ifelse(isTRUE(fixed) & isTRUE(random), 
                     pval.fixed, pval.random) %>% 
            scales::pvalue(),
          i2 = round(I2*100, 2),
          i2.ci = paste0("[", round(lower.I2*100, 2), "; ", 
                         round(upper.I2*100, 2), "]"),
          prediction.ci = paste0("[", 
                                 round(
                                   ifelse(identical(.type.es, "RR"), 
                                          exp(lower.predict), lower.predict), 
                                   round.digits), "; ",
                                 round(
                                   ifelse(identical(.type.es, "RR"), 
                                          exp(upper.predict), upper.predict), 
                                   round.digits), "]"),
          nnt = metapsyNNT(
            ifelse(isTRUE(fixed) & isTRUE(random), 
                   TE.fixed, TE.random), nnt.cer) %>%
            ifelse(identical(.type.es, "RR"), NA, .) %>% 
            ifelse(isTRUE(.raw.bin.es), 
                   nnt.raw.bin.es, .) %>% 
            round(round.digits) %>% abs(),
          excluded = "none")
      })
    } else {
      mLowest = mGeneral
      mLowestRes = mGeneralRes
    }
  
    if (is.null(mLowest$exclude[1])){
      mLowestRes$excluded = "none"
    } else {
      mLowestRes$excluded = paste(data$comparison[mLowest$exclude],
                                  collapse = "; ")
    }
    rownames(mLowestRes) = "One ES/study (lowest)"
  
    if ("lowest.highest" %in% which.run){
      message(crayon::green("DONE"))
    }
  
  } else {
    mLowest = NULL
    mLowestRes = NULL
  }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  3. Highest Only                                                          #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  if ("lowest.highest" %in% which.run){
    
    message("- ", crayon::green("[OK] "), 
            "Calculating effect size using only highest effect... ",
            appendLF = FALSE)
  
  
    multi.study = names(table(data[[study.var]])
                        [table(data[[study.var]]) > 1])
  
    if (length(multi.study) > 0){
      if (length(multi.study) > 0){
        data$.TE = mGeneral$TE
        data %>%
          split(.[[study.var]]) %>%
          purrr::map(function(x){
            tiebreaker = rnorm(nrow(x), sd=1e-10)
            x$.TE + tiebreaker == max(x$.TE + tiebreaker)}) %>%
          do.call(c,.) -> highest
        data$.TE = NULL
      } else {
        highest = NULL
      }
  
      mHighest = meta::update.meta(mGeneral, exclude = !highest,
                                   id = NULL)
      with(update.meta(mHighest, sm="RD",
                       id = NULL), {
        ifelse(isTRUE(fixed) & isTRUE(random),
               abs(TE.fixed)^-1, abs(TE.random)^-1)
      }) -> nnt.raw.bin.es
      
      mHighestRes = with(mHighest,{
        data.frame(
          k = k,
          g = ifelse(
            isTRUE(fixed) & isTRUE(random), 
            TE.fixed, TE.random) %>% 
            ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
            round(round.digits),
          g.ci = paste0("[", 
                        ifelse(isTRUE(fixed) & isTRUE(random),
                               lower.fixed, lower.random) %>% 
                          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                          round(round.digits), "; ",
                        ifelse(isTRUE(fixed) & isTRUE(random),
                               upper.fixed, upper.random) %>% 
                          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                          round(round.digits), "]"),
          p = ifelse(isTRUE(fixed) & isTRUE(random), 
                     pval.fixed, pval.random) %>% 
            scales::pvalue(),
          i2 = round(I2*100, 2),
          i2.ci = paste0("[", round(lower.I2*100, 2), "; ", 
                         round(upper.I2*100, 2), "]"),
          prediction.ci = paste0("[", 
                                 round(
                                   ifelse(identical(.type.es, "RR"), 
                                          exp(lower.predict), lower.predict), 
                                   round.digits), "; ",
                                 round(
                                   ifelse(identical(.type.es, "RR"), 
                                          exp(upper.predict), upper.predict), 
                                   round.digits), "]"),
          nnt = metapsyNNT(
            ifelse(isTRUE(fixed) & isTRUE(random), 
                   TE.fixed, TE.random), nnt.cer) %>%
            ifelse(identical(.type.es, "RR"), NA, .) %>% 
            ifelse(isTRUE(.raw.bin.es), 
                   nnt.raw.bin.es, .) %>% 
            round(round.digits) %>% abs(),
          excluded = "none")
      })
    } else {
      mHighest = mGeneral
      mHighestRes = mGeneralRes
    }
  
    if (is.null(mHighest$exclude[1])){
      mHighestRes$excluded = "none"
    } else {
      mHighestRes$excluded = 
        paste(data$comparison[mHighest$exclude],
              collapse = "; ")}
    rownames(mHighestRes) = "One ES/study (highest)"
  
    if ("lowest.highest" %in% which.run){
      message(crayon::green("DONE"))
    }
    
  } else {
    mHighest = NULL
    mHighestRes = NULL
  }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  4. Combined Effect Sizes                                                 #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  
  if ("combined" %in% which.run){
    message("- ", crayon::green("[OK] "), 
            "Calculating effect size using combined effects (rho=", 
            rho.within.study, ")... ",
            appendLF = FALSE)

  
    if (!is.null(extra.grouping.var)){
      data$study.var.comb = paste(data[[study.var]], 
                                  data[[extra.grouping.var]])
      study.var.comb = "study.var.comb"
    } else {
      data$study.var.comb = data[[study.var]]
      study.var.comb = study.var
    }
  
    # Define study ID variable
    svc = data[[study.var.comb]]
    multi.study.comb = names(table(data[[study.var.comb]])
                             [table(data[[study.var.comb]]) > 1])
    if (identical(.type.es, "g")){
      data = metafor::escalc("SMD", yi = data[[es.var[1]]],
                             sei = data[[se.var[1]]], data = data)
    } else {
      if (isTRUE(.raw.bin.es)){
       data = metafor::escalc("RR", 
                        ai = data[[es.binary.raw.vars[1]]], 
                        ci = data[[es.binary.raw.vars[2]]], 
                        n1i = data[[es.binary.raw.vars[3]]], 
                        n2i = data[[es.binary.raw.vars[4]]], 
                        data = data)
      } else {
        data = metafor::escalc("RR", yi = data[[es.var[1]]],
                               sei = data[[se.var[1]]], data = data)
      }
    }
    
    rho.within.study = abs(rho.within.study[1])
    if (rho.within.study[1] > 1){
      stop("'rho.within.study' must be in [-1,1].")
    }
    if (rho.within.study[1] >.99){
      message("- ", crayon::yellow("[!] "), 
          "'rho.within.study' is very close to 1.",
          " This can lead to non-positive-definite variance-covariance matrices.",
          " If the V matrix is not positive-definite, assume a lower value...")
      warn.end = TRUE}
    if (rho.within.study[1] <.01){
      message("- ", crayon::yellow("[!] "), 
          "'rho.within.study' is very close to 0.",
          " This can lead to non-positive-definite variance-covariance matrices.",
          " If the V matrix is not positive-definite, assume a higher value...")
      warn.end = TRUE
    }
    
    data.comb = metafor::aggregate.escalc(data, cluster = svc,
                                          rho = rho.within.study)
  
    # Get studlabs of studies with multiple arms
    table(data.comb[[study.var]], 
          data.comb[[study.var.comb]]) %>%
      {names(rowSums(.)[rowSums(.) > 1])} -> multi.entries
    
    mComb = meta::metagen(TE = data.comb$yi,
                          seTE = sqrt(data.comb$vi),
                          studlab = with(data.comb,{
                            ifelse(study %in% multi.entries,
                                   study.var.comb, study)}),
                          data = data.comb,
                          sm = .type.es,
                          hakn = hakn,
                          method.tau = method.tau.meta,
                          method.tau.ci = method.tau.ci,
                          prediction = TRUE,
                          fixed = ifelse(method.tau == "FE", TRUE, FALSE),
                          random = ifelse(method.tau == "FE", FALSE, TRUE),
                          ...)
    
    if (isTRUE(.raw.bin.es)){
      meta::metabin(
       event.e = data.comb[[es.binary.raw.vars[1]]],
       n.e = data.comb[[es.binary.raw.vars[3]]],
       event.c = data.comb[[es.binary.raw.vars[2]]],
       n.c = data.comb[[es.binary.raw.vars[4]]],
       studlab = data.comb[[study.var]],
       data = data.comb,
       sm = "RD",
       hakn = hakn,
       method.tau = method.tau.meta,
       method.tau.ci = method.tau.ci,
       prediction = TRUE,
       fixed = ifelse(method.tau == "FE", TRUE, FALSE),
       random = ifelse(method.tau == "FE", FALSE, TRUE)) %>% 
        with(., {
          ifelse(isTRUE(fixed) & isTRUE(random),
                 abs(TE.fixed)^-1, abs(TE.random)^-1)
        }) -> nnt.raw.bin.es
    }
    mComb$data$.event.e = data.comb[[es.binary.raw.vars[1]]]
    mComb$data$.n.e = data.comb[[es.binary.raw.vars[3]]]
    mComb$data$.event.c = data.comb[[es.binary.raw.vars[2]]]
    mComb$data$.n.c = data.comb[[es.binary.raw.vars[4]]]
  
    
    mCombRes = with(mComb,{
  
      data.frame(
        k = k,
        g = ifelse(
          isTRUE(fixed) & isTRUE(random), 
          TE.fixed, TE.random) %>% 
          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
          round(round.digits),
        g.ci = paste0("[", 
                      ifelse(isTRUE(fixed) & isTRUE(random),
                             lower.fixed, lower.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "; ",
                      ifelse(isTRUE(fixed) & isTRUE(random),
                             upper.fixed, upper.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "]"),
        p = ifelse(isTRUE(fixed) & isTRUE(random), 
                   pval.fixed, pval.random) %>% 
          scales::pvalue(),
        i2 = round(I2*100, 2),
        i2.ci = paste0("[", round(lower.I2*100, 2), "; ", 
                       round(upper.I2*100, 2), "]"),
        prediction.ci = paste0("[", 
                               round(
                                 ifelse(identical(.type.es, "RR"), 
                                        exp(lower.predict), lower.predict), 
                                 round.digits), "; ",
                               round(
                                 ifelse(identical(.type.es, "RR"), 
                                        exp(upper.predict), upper.predict), 
                                 round.digits), "]"),
        nnt = metapsyNNT(
          ifelse(isTRUE(fixed) & isTRUE(random), 
                 TE.fixed, TE.random), nnt.cer) %>%
          ifelse(identical(.type.es, "RR"), NA, .) %>% 
          ifelse(isTRUE(.raw.bin.es), 
                 nnt.raw.bin.es, .) %>% 
          round(round.digits) %>% abs(),
        excluded = "none")
    })
    rownames(mCombRes) = "Combined"
    mCombRes$excluded = paste("combined:", 
                              ifelse(length(multi.study.comb) > 0,
                                     paste(multi.study.comb, collapse = "; "), 
                                     "none"))
  
    # Add rho to mComb model
    mComb$rho = rho.within.study
    class(data) = "data.frame"
    if ("combined" %in% which.run){
      message(crayon::green("DONE"))
    }
    
  } else {
    mComb = NULL
    mCombRes = NULL
  }
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  5. Outliers Removed                                                      #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  if ("outliers" %in% which.run){
    message("- ", crayon::green("[OK] "), 
            "Calculating effect size with outliers removed... ",
            appendLF = FALSE)
  
    if (!("overall" %in% which.run) && 
        ("combined" %in% which.run)){
      which.outliers[1] = "combined"
    }
  
    if (which.outliers[1] %in% c("general", "overall")){
      if (is.null(mGeneral)) {
        stop("If which.outliers = 'overall', the 'overall' model must be",
             " included in which.run.", call. = FALSE)
      }
      m.for.outliers = mGeneral
    }

    if (which.outliers[1] == "combined"){
      if (is.null(mComb)) {
        stop("If which.outliers = 'combined', the 'combined' model must be",
             " included in which.run.", call. = FALSE)
      }
      m.for.outliers = mComb
    }

    if (!(which.outliers[1] %in% c("general", "combined", "overall"))){
      stop("'which.outliers' must either be 'overall' or 'combined'.")
    }

    if (method.tau == "FE"){
      mOutliers = metapsyFindOutliers(m.for.outliers)$m.fixed
    } else {
      mOutliers = metapsyFindOutliers(m.for.outliers)$m.random
    }
  
    if (isTRUE(.raw.bin.es)){
      with(update.meta(mOutliers, sm="RD",
                       id = NULL), {
        ifelse(isTRUE(fixed) & isTRUE(random),
               abs(TE.fixed)^-1, abs(TE.random)^-1)
      }) -> nnt.raw.bin.es
    }

    mOutliersRes = with(mOutliers,{

      data.frame(
        k = k,
        g = ifelse(
          isTRUE(fixed) & isTRUE(random), 
          TE.fixed, TE.random) %>% 
          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
          round(round.digits),
        g.ci = paste0("[", 
                      ifelse(isTRUE(fixed) & isTRUE(random),
                             lower.fixed, lower.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "; ",
                      ifelse(isTRUE(fixed) & isTRUE(random),
                             upper.fixed, upper.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "]"),
        p = ifelse(isTRUE(fixed) & isTRUE(random), 
                   pval.fixed, pval.random) %>% 
          scales::pvalue(),
        i2 = round(I2*100, 2),
        i2.ci = paste0("[", round(lower.I2*100, 2), "; ", 
                       round(upper.I2*100, 2), "]"),
        prediction.ci = paste0("[", 
                               round(
                                 ifelse(identical(.type.es, "RR"), 
                                        exp(lower.predict), lower.predict), 
                                 round.digits), "; ",
                               round(
                                 ifelse(identical(.type.es, "RR"), 
                                        exp(upper.predict), upper.predict), 
                                 round.digits), "]"),
        nnt = metapsyNNT(
          ifelse(isTRUE(fixed) & isTRUE(random), 
                 TE.fixed, TE.random), nnt.cer) %>%
          ifelse(identical(.type.es, "RR"), NA, .) %>% 
          ifelse(isTRUE(.raw.bin.es), 
            nnt.raw.bin.es, .) %>% 
          round(round.digits) %>% abs(),
        excluded = "none")
    })
    rownames(mOutliersRes) = "Outliers removed"

    if (length(mOutliers$data$comparison
               [mOutliers$exclude]) != 0){
      if (identical(which.outliers[1], "combined")){
        mOutliersRes$excluded = paste(mOutliers$studlab
                                      [mOutliers$exclude], collapse = "; ")
      } else {
        mOutliersRes$excluded = paste(mOutliers$data$comparison
                                      [mOutliers$exclude], collapse = "; ")
      }
    } else {
      mOutliersRes$excluded = "no outliers detected"
    }

    if ("outliers" %in% which.run){
      message(crayon::green("DONE"))
    }
    
  } else {
    mOutliers = NULL
    mOutliersRes = NULL
  }


  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  6. Influential Cases                                                     #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  if ("influence" %in% which.run){
    message("- ", crayon::green("[OK] "), 
            "Calculating effect size with influential cases removed... ",
            appendLF = FALSE)
    
    if (!("overall" %in% which.run) && 
        ("combined" %in% which.run)){
      which.influence[1] = "combined"
    }
    
    if (which.influence[1] %in% c("general", "overall")){
      if (is.null(mGeneral)) {
        stop("If which.influence = 'overall', the 'overall' model must be",
             " included in which.run.", call. = FALSE)
      }
      m.for.influence = mGeneral
    }
  
    if (which.influence[1] == "combined"){
      if (is.null(mComb)) {
        stop("If which.influence = 'combined', the 'combined' model must be",
             " included in which.run.", call. = FALSE)
      }
      m.for.influence = mComb
    }
  
    if (!(which.influence[1] %in% c("general", "combined", "overall"))){
      stop("'which.influence' must either be 'overall' or 'combined'.")
    }
  
    influenceRes = 
      metapsyInfluenceAnalysis(
        m.for.influence,
        random = ifelse(method.tau == "FE",
                        FALSE, TRUE))
    mInfluence = 
      meta::update.meta(
        m.for.influence,
        exclude = influenceRes$Data$is.infl == "yes",
        id = NULL)
    
    if (isTRUE(.raw.bin.es)){
      with(update.meta(mInfluence, sm="RD",
                       id = NULL), {
        ifelse(isTRUE(fixed) & isTRUE(random),
               abs(TE.fixed)^-1, abs(TE.random)^-1)
      }) -> nnt.raw.bin.es
    }
    
    mInfluenceRes = with(mInfluence,{
  
      data.frame(
        k = k,
        g = ifelse(
          isTRUE(fixed) & isTRUE(random), 
          TE.fixed, TE.random) %>% 
          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
          round(round.digits),
        g.ci = paste0("[", 
                      ifelse(isTRUE(fixed) & isTRUE(random),
                             lower.fixed, lower.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "; ",
                      ifelse(isTRUE(fixed) & isTRUE(random),
                             upper.fixed, upper.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "]"),
        p = ifelse(isTRUE(fixed) & isTRUE(random), 
                   pval.fixed, pval.random) %>% 
          scales::pvalue(),
        i2 = round(I2*100, 2),
        i2.ci = paste0("[", round(lower.I2*100, 2), "; ", 
                       round(upper.I2*100, 2), "]"),
        prediction.ci = paste0("[", 
                               round(
                                 ifelse(identical(.type.es, "RR"), 
                                        exp(lower.predict), lower.predict), 
                                 round.digits), "; ",
                               round(
                                 ifelse(identical(.type.es, "RR"), 
                                        exp(upper.predict), upper.predict), 
                                 round.digits), "]"),
        nnt = metapsyNNT(
          ifelse(isTRUE(fixed) & isTRUE(random), 
                 TE.fixed, TE.random), nnt.cer) %>%
          ifelse(identical(.type.es, "RR"), NA, .) %>% 
          ifelse(isTRUE(.raw.bin.es), 
                 nnt.raw.bin.es, .) %>% 
            round(round.digits) %>% abs(),
        excluded = "none")
    })
    rownames(mInfluenceRes) = "Influence Analysis"
    
    if (identical(which.influence[1], "combined")){
      mInfluenceRes$excluded = 
        paste("removed as influential cases:",
              ifelse(sum(influenceRes$Data$is.infl == "yes") > 0,
                     paste(mInfluence$studlab[influenceRes$Data$is.infl == "yes"],
                           collapse = "; "), "none"))
    } else {
      mInfluenceRes$excluded = 
        paste("removed as influential cases:",
              ifelse(sum(influenceRes$Data$is.infl == "yes") > 0,
                     paste(mInfluence$data$comparison[influenceRes$Data$is.infl == "yes"],
                           collapse = "; "), "none"))
    }
    
    if ("influence" %in% which.run){
      message(crayon::green("DONE"))
    }
    
  } else {
    mInfluence = NULL
    mInfluenceRes = NULL
    influenceRes = NULL
  }



  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  7. (Low) risk of bias                                                    #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  if ("rob" %in% which.run){
    message("- ", crayon::green("[OK] "), 
            "Calculating effect size using only low RoB information... ",
            appendLF = FALSE)
  
    if (!("overall" %in% which.run) && 
        ("combined" %in% which.run)){
      which.rob[1] = "combined"
    }
    
    if (which.rob[1] %in% c("general", "overall")){
      if (is.null(mGeneral)) {
        stop("If which.influence = 'overall', the 'overall' model must be",
             " included in which.run.", call. = FALSE)
      }
      m.for.rob = mGeneral
      data.for.rob = mGeneral$data
    }
  
    if (which.rob[1] == "combined"){
      if (is.null(mComb)) {
        stop("If which.influence = 'overall', the 'overall' model must be",
             " included in which.run.", call. = FALSE)
      }
      m.for.rob = mComb
      data.for.rob = mComb$data
    }
  
    if (!(which.rob[1] %in% c("general", "combined", "overall"))){
      stop("'which.rob' must either be 'general' or 'combined'.")
    }
  
    if ("rob" %in% which.run){
      robVar = strsplit(low.rob.filter, " ")[[1]][1]
      m.for.rob$data[[robVar]] = 
        as.numeric(m.for.rob$data[[robVar]])
      robFilter = paste0("data.for.rob$", low.rob.filter, 
                         " & !is.na(data.for.rob$", robVar, ")")
      robMask = eval(parse(text = robFilter))
      if (sum(robMask) == 0){
        message("\n- ", crayon::yellow("[!] "), 
          "No low risk of bias studies detected! Switching to 'general'...")
        robMask = rep(TRUE, nrow(mGeneral$data))
        mRob = meta::update.meta(mGeneral, exclude = !robMask,
                                 id = NULL)
        which.run[!which.run == "rob"] -> which.run
        warn.end = TRUE
      } else {
        mRob = meta::update.meta(m.for.rob, exclude = !robMask,
                                 id = NULL)
      }
    } else {
      robMask = rep(TRUE, nrow(mGeneral$data))
      mRob = meta::update.meta(mGeneral, exclude = !robMask,
                               id = NULL)
    }
    
    if (isTRUE(.raw.bin.es)) {
     with(update.meta(mRob, sm="RD",
                      id = NULL), {
       ifelse(isTRUE(fixed) & isTRUE(random),
              abs(TE.fixed)^-1, abs(TE.random)^-1)
     }) -> nnt.raw.bin.es
    }
  
    mRobRes = with(mRob,{
  
      data.frame(
        k = k,
        g = ifelse(
          isTRUE(fixed) & isTRUE(random), 
          TE.fixed, TE.random) %>% 
          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
          round(round.digits),
        g.ci = paste0("[", 
                      ifelse(isTRUE(fixed) & isTRUE(random),
                             lower.fixed, lower.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "; ",
                      ifelse(isTRUE(fixed) & isTRUE(random),
                             upper.fixed, upper.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "]"),
        p = ifelse(isTRUE(fixed) & isTRUE(random), 
                   pval.fixed, pval.random) %>% 
          scales::pvalue(),
        i2 = round(I2*100, 2),
        i2.ci = paste0("[", round(lower.I2*100, 2), "; ", 
                       round(upper.I2*100, 2), "]"),
        prediction.ci = paste0("[", 
                               round(
                                 ifelse(identical(.type.es, "RR"), 
                                        exp(lower.predict), lower.predict), 
                                 round.digits), "; ",
                               round(
                                 ifelse(identical(.type.es, "RR"), 
                                        exp(upper.predict), upper.predict), 
                                 round.digits), "]"),
        nnt = metapsyNNT(
          ifelse(isTRUE(fixed) & isTRUE(random), 
                 TE.fixed, TE.random), nnt.cer) %>%
          ifelse(identical(.type.es, "RR"), NA, .) %>% 
          ifelse(isTRUE(.raw.bin.es), 
                 nnt.raw.bin.es, .) %>% 
          round(round.digits) %>% abs(),
        excluded = "none")
    })
    rownames(mRobRes) = paste("Only", low.rob.filter)
  
    if (length(mRob$data$comparison[!robMask]) != 0){
      if (identical(which.outliers[1], "combined")){
        mRobRes$excluded = paste(mRob$studlab[!robMask], collapse = "; ")
      } else {
        mRobRes$excluded = 
          paste(mRob$data$comparison[!robMask], collapse = "; ")
      }
    } else {
      mRobRes$excluded = paste0("no studies removed; ", 
                                low.rob.filter, " applies to all studies")
    }
    
    if ("rob" %in% which.run){
      message(crayon::green("DONE"))
    }
  } else {
    mRob = NULL
    mRobRes = NULL
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
    

    data$es.id = 1:nrow(data)
    formula = paste0("~ 1 | ", colnames(data[study.var]), " / es.id")
  
    mThreeLevel = metafor::rma.mv(
                   yi = mGeneral[["TE"]],
                   V = mGeneral[["seTE"]]^2,
                   slab = mGeneral[["studlab"]],
                   data = data,
                   random = as.formula(formula),
                   test = ifelse(hakn == TRUE, "t", "z"),
                   method = "REML",
                   ...)
    
    model.threelevel.legacy = list(
      slab = data[[study.var]],
      data = data,
      es.var = es.var[1],
      se.var = se.var[1],
      yi = mGeneral[["TE"]],
      V = mGeneral[["seTE"]]*mGeneral[["seTE"]],
      formula.rnd = as.formula(formula))
    mThreeLevel$legacy = model.threelevel.legacy
    
    if (isTRUE(.raw.bin.es)){
      # For RR analyses: re-run analyses using g
      # This is needed for NNTs
      meta::metaprop(
        event = data[[es.binary.raw.vars[2]]],
        n = data[[es.binary.raw.vars[4]]],
        fixed = ifelse(method.tau == "FE", 
                       TRUE, FALSE),
        random = ifelse(method.tau == "FE", 
                        FALSE, TRUE)) %>%
        {ifelse(isTRUE(.$fixed) & isTRUE(.$random), 
                .$TE.fixed, .$TE.random)} %>% 
        {exp(.)/(1+exp(.))} -> cer          
      
      apply(data, 1, function(x){
        esc::esc_2x2(grp1yes = as.numeric(x[[es.binary.raw.vars[1]]]), 
                     grp1no = as.numeric(x[[es.binary.raw.vars[3]]]) - 
                       as.numeric(x[[es.binary.raw.vars[1]]]),
                     grp2yes = as.numeric(x[[es.binary.raw.vars[2]]]),
                     grp2no = as.numeric(x[[es.binary.raw.vars[4]]]) -
                       as.numeric(x[[es.binary.raw.vars[2]]]),
                     es.type = "g") %>% 
          {data.frame(es = .$es, se = .$se)}
      }) %>% 
        do.call(rbind, .) -> data.g
      
      metafor::rma.mv(
        yi = data.g$es, V = data.g$se^2, 
        data = data, random = as.formula(formula),
        test = ifelse(hakn == TRUE, "t", "z"),
        method = "REML") -> mThreeLevel.g
      
      mThreeLevel.g %>% 
        {.[["b"]][1]} %>% 
        {ifelse(.==0, Inf, 
                metapsyNNT(., cer))} -> nnt.g
    } else {
      nnt.g = NA
    }
  
    # Calculate total I2
    W = diag(1/(mGeneral[["seTE"]]^2))
    X = model.matrix(mThreeLevel)
    P = W - W %*% X %*% 
      solve(t(X) %*% W %*% X) %*% 
      t(X) %*% W
    with(mThreeLevel, {
      100 * sum(sigma2) / 
        (sum(sigma2) + (k-p)/sum(diag(P)))
    }) -> mThreeLevel$I2
    
    # Calculate I2 per level
    with(mThreeLevel, {
      (100 * sigma2 / 
       (sum(sigma2) + (k-p)/sum(diag(P))))
    }) -> I2.bw
    mThreeLevel$I2.between.studies = I2.bw[1]
    mThreeLevel$I2.within.studies = I2.bw[2]
  
    # Get tau and I2
    with(mThreeLevel, 
         {data.frame(tau2 = c(sigma2, sum(sigma2)),
                     i2 = c(I2.between.studies, 
                            I2.within.studies,
                            I2))}) -> mThreeLevel$variance.components
    
    mThreeLevel$variance.components$tau2 = 
      round(mThreeLevel$variance.components$tau2, 4)
    mThreeLevel$variance.components$i2 = 
      round(mThreeLevel$variance.components$i2, 1)
    rownames(mThreeLevel$variance.components) = 
      c("Between Studies", "Within Studies", "Total")
  
    if (use.rve[1] == FALSE){
      
      if ("threelevel" %in% which.run){
        message("- ", crayon::green("[OK] "), 
                "Calculating effect size using three-level MA model... ",
                appendLF = FALSE)
      }
      
      mThreeLevelRes = with(mThreeLevel, {
        data.frame(k = k.all,
                   g = as.numeric(b[,1]) %>%
                     ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                     round(round.digits),
                   g.ci = paste0("[", 
                                 ci.lb %>% 
                                   ifelse(identical(.type.es, "RR"), 
                                          exp(.), .) %>% 
                                   round(round.digits), "; ",
                                 ci.ub %>% 
                                   ifelse(identical(.type.es, "RR"), 
                                          exp(.), .) %>% 
                                   round(round.digits), "]"),
                   p = pval %>% scales::pvalue(),
                   i2 = round(I2, 1),
                   i2.ci = "-",
                   prediction.ci = paste0("[", 
                         round(predict(mThreeLevel)$pi.lb %>% 
                                 ifelse(identical(.type.es, "RR"), exp(.), .), 
                               round.digits), "; ",
                         round(predict(mThreeLevel)$pi.ub %>% 
                                 ifelse(identical(.type.es, "RR"), exp(.), .), 
                               round.digits), "]"),
                   nnt = ifelse(identical(.type.es, "RR"), 
                                NA, metapsyNNT(as.numeric(b[,1]), nnt.cer)) %>% 
                         ifelse(isTRUE(.raw.bin.es), nnt.g, .) %>% 
                         round(round.digits) %>% abs(),
                   excluded = "none")
      })
      rownames(mThreeLevelRes) = "Three-Level Model"
      mThreeLevelRes$excluded = 
        paste("Number of clusters/studies:", 
              mThreeLevel$s.nlevels[1])
  
      if ("threelevel" %in% which.run){
        message(crayon::green("DONE"))
      }
  
      if (length(multi.study) == 0 & 
          "threelevel" %in% which.run){
        message("- ", crayon::yellow("[!] "), 
                "All included ES seem to be independent.",
                " A three-level model is not adequate and ",
                " tau/I2 estimates are not trustworthy!")
        warn.end = TRUE
        which.run[!which.run == "threelevel"] -> which.run
      }
  
    } else {
      
      if ("threelevel" %in% which.run){
        message("- ", crayon::green("[OK] "),
          "Calculating effect size using three-level MA model... ",
          appendLF = FALSE)
        message(crayon::green("DONE"))
        message("- ", crayon::green("[OK] "),
          "Robust variance estimation (RVE) used for three-level MA model... ",
          appendLF = FALSE)
      }
  
      # Get results: RVE used
      crit = qt(0.025, mThreeLevel$ddf[[1]], lower.tail = FALSE)
      tau2 = sum(mThreeLevel$sigma2)
      SE = clubSandwich::conf_int(mThreeLevel, vcov = "CR2")[["SE"]]
      pi.lb.rve = 
        clubSandwich::conf_int(mThreeLevel, "CR2")[["beta"]] -
        crit * sqrt(tau2 + (SE^2))
      pi.ub.rve = 
        clubSandwich::conf_int(mThreeLevel, "CR2")[["beta"]] +
        crit * sqrt(tau2 + (SE^2))
      if (isTRUE(.raw.bin.es)){
        as.numeric(clubSandwich::conf_int(
          mThreeLevel.g, "CR2")[["beta"]]) %>% 
          {ifelse(.==0, Inf, metapsyNNT(., cer))} -> nnt.g
      }
  
      mThreeLevelRes.RVE = with(mThreeLevel, {
        data.frame(k = k.all,
                   g = as.numeric(
                     clubSandwich::conf_int(
                       mThreeLevel, "CR2")[["beta"]]) %>%
                     ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                     round(round.digits),
                   g.ci = paste0(
                     "[",
                     clubSandwich::conf_int(
                       mThreeLevel, "CR2")[["CI_L"]] %>%
                      ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                      round(round.digits), "; ",
                     clubSandwich::conf_int(
                       mThreeLevel, "CR2")[["CI_U"]] %>%
                      ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                      round(round.digits), "]"),
                   p = clubSandwich::coef_test(
                     mThreeLevel, "CR2")[["p_Satt"]] %>%
                     scales::pvalue(),
                   i2 = round(I2, 1),
                   i2.ci = "-",
                   prediction.ci = paste0(
                     "[", 
                     round(pi.lb.rve %>% 
                           ifelse(identical(.type.es, "RR"), exp(.), .), 
                          round.digits), "; ",
                     round(pi.ub.rve %>% 
                           ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "]"),
                   nnt = ifelse(identical(.type.es, "RR"), 
                          NA, metapsyNNT(as.numeric(
                            clubSandwich::conf_int(
                              mThreeLevel, "CR2")[["beta"]]), nnt.cer)) %>% 
                     ifelse(isTRUE(.raw.bin.es), nnt.g, .) %>% 
                     round(round.digits) %>% abs(),
                   excluded = "none")
      })
      rownames(mThreeLevelRes.RVE) = "Three-Level Model"
      mThreeLevelRes.RVE$excluded = 
        paste0("Number of clusters/studies: ", 
               mThreeLevel$s.nlevels[1], 
               "; robust variance estimation (RVE) used.")
      mThreeLevelRes = mThreeLevelRes.RVE
  
      if ("threelevel" %in% which.run){
        message(crayon::green("DONE"))
      }
  
      if (length(multi.study) == 0 & "threelevel" %in% which.run){
        message("- ", crayon::yellow("[!] "), 
                "All included ES seem to be independent.",
                " A three-level model is not adequate and ",
                " tau/I2 estimates are not trustworthy!")
        warn.end = TRUE
        which.run[!which.run == "threelevel"] -> which.run
      }
    }
  } else {
    mThreeLevel = NULL
    mThreeLevelRes = NULL
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

    data$es.id = 1:nrow(data)
    formula.fixed = as.formula("TE ~ 1")
    formula.rnd = as.formula(paste0("~ 1 | ", 
                     colnames(data[study.var]), "/ es.id"))
    data.che = data
    data.che$TE = mGeneral[["TE"]]
  
    Vmat = clubSandwich::impute_covariance_matrix(
              mGeneral[["seTE"]]^2,
              cluster = data[[study.var]],
              r = rho.within.study,
              smooth_vi = TRUE)
  
    mCHE = metafor::rma.mv(
      formula.fixed,
      V = Vmat,
      slab = data[[study.var]],
      data = data.che,
      random = formula.rnd,
      test = ifelse(hakn == TRUE, "t", "z"),
      method = "REML", 
      sparse = TRUE,
      ...)
  
    model.threelevel.che.legacy = list(
      slab = data.che[[study.var]],
      data = data.che,
      formula.rnd = formula.rnd,
      formula.fixed = formula.fixed,
      Vmat = Vmat)
    mCHE$legacy = model.threelevel.che.legacy
    
    if (isTRUE(.raw.bin.es)){
      # For RR analyses: re-run analyses using g
      # This is needed for NNTs
      meta::metaprop(
        event = data.che[[es.binary.raw.vars[2]]],
        n = data.che[[es.binary.raw.vars[4]]],
        fixed = ifelse(method.tau == "FE", 
                       TRUE, FALSE),
        random = ifelse(method.tau == "FE", 
                        FALSE, TRUE)) %>%
        {ifelse(isTRUE(.$fixed) & isTRUE(.$random), 
                .$TE.fixed, .$TE.random)} %>% 
        {exp(.)/(1+exp(.))} -> cer          
      
      apply(data.che, 1, function(x){
        esc::esc_2x2(grp1yes = as.numeric(x[[es.binary.raw.vars[1]]]), 
                     grp1no = as.numeric(x[[es.binary.raw.vars[3]]]) - 
                       as.numeric(x[[es.binary.raw.vars[1]]]),
                     grp2yes = as.numeric(x[[es.binary.raw.vars[2]]]),
                     grp2no = as.numeric(x[[es.binary.raw.vars[4]]]) -
                       as.numeric(x[[es.binary.raw.vars[2]]]),
                     es.type = "g") %>% 
          {data.frame(es = .$es, se = .$se)}
      }) %>% 
        do.call(rbind, .) -> data.g
      
      Vmat.g = clubSandwich::impute_covariance_matrix(
        data.g$se^2,
        cluster = data.che[[study.var]],
        r = rho.within.study,
        smooth_vi = TRUE)
      
      data.che$TE = data.g$es
      metafor::rma.mv(
        formula.fixed,
        V = Vmat.g,
        slab = data.che[[study.var]],
        data = data.che,
        random = formula.rnd,
        test = ifelse(hakn == TRUE, "t", "z"),
        method = "REML", 
        sparse = TRUE) -> mCHE.g
      
      mCHE.g %>% 
        {.[["b"]][1]} %>% 
        {ifelse(.==0, Inf, 
                metapsyNNT(., cer))} -> nnt.g
    } else {
      nnt.g = NA
    }
  
    # Calculate total I2
    W = diag(1/(data[[se.var]]^2))
    X = model.matrix(mCHE)
    P = W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    mCHE$I2 = 100 * sum(mCHE$sigma2) /
      (sum(mCHE$sigma2) + (mCHE$k-mCHE$p)/sum(diag(P)))
  
    # Calculate I2 per level
    with(mCHE, {
      (100 * sigma2 / 
         (sum(sigma2) + (k-p)/sum(diag(P))))
    }) -> I2.bw
    mCHE$I2.between.studies = I2.bw[1]
    mCHE$I2.within.studies = I2.bw[2]
  
    # Get tau and I2
    with(mCHE, 
         {data.frame(tau2 = c(sigma2, sum(sigma2)),
                     i2 = c(I2.between.studies, 
                            I2.within.studies,
                            I2))}) -> mCHE$variance.components
    
    mCHE$variance.components$tau2 = 
      round(mCHE$variance.components$tau2, 4)
    mCHE$variance.components$i2 = 
      round(mCHE$variance.components$i2, 1)
    rownames(mCHE$variance.components) = 
      c("Between Studies", "Within Studies", "Total")
  
    # Get results: no RVE
    if (use.rve[1] == FALSE){
      if ("threelevel.che" %in% which.run){
        message("- ", crayon::green("[OK] "), 
                "Calculating effect size using three-level CHE model... ",
                appendLF = FALSE)
      }
      mCHERes = with(mCHE, {
        data.frame(k = k.all,
                   g = as.numeric(b[,1]) %>%
                     ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                     round(round.digits),
                   g.ci = paste0("[", 
                                 ci.lb %>% 
                                   ifelse(identical(.type.es, "RR"), 
                                          exp(.), .) %>% 
                                   round(round.digits), "; ",
                                 ci.ub %>% 
                                   ifelse(identical(.type.es, "RR"), 
                                          exp(.), .) %>% 
                                   round(round.digits), "]"),
                   p = pval %>% scales::pvalue(),
                   i2 = round(I2, 1),
                   i2.ci = "-",
                   prediction.ci = paste0(
                     "[", 
                     round(predict(mCHE)$pi.lb %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "; ",
                     round(predict(mCHE)$pi.ub %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "]"),
                   nnt = ifelse(identical(.type.es, "RR"), 
                                NA, metapsyNNT(as.numeric(b[,1]), nnt.cer)) %>% 
                     ifelse(isTRUE(.raw.bin.es), nnt.g, .) %>% 
                     round(round.digits) %>% abs(),
                   excluded = "none")
      })
      rownames(mCHERes) = "Three-Level Model (CHE)"
      mCHERes$excluded = paste("Number of clusters/studies:", 
                               mCHE$s.nlevels[1])
  
      if ("threelevel.che" %in% which.run){
        message(crayon::green("DONE"))
      }
  
      if (length(multi.study) == 0 & "threelevel.che" %in% which.run){
        message("- ", crayon::yellow("[!] "), 
                "All included ES seem to be independent.",
                " A three-level CHE model is not adequate and ",
                "tau/I2 estimates are not trustworthy!")
        warn.end = TRUE
        which.run[!which.run == "threelevel.che"] -> which.run
      }
    } else {
      
      if ("threelevel.che" %in% which.run){
        message("- ", crayon::green("[OK] "), 
          "Calculating effect size using three-level CHE model (rho=", 
                rho.within.study, ")... ",
          appendLF = FALSE)
        message(crayon::green("DONE"))
        message("- ", crayon::green("[OK] "), 
          "Robust variance estimation (RVE) used for three-level CHE model... ",
          appendLF = FALSE)
      }
  
      # Get results: RVE used
      crit = qt(0.025, mCHE$ddf[[1]], lower.tail = FALSE)
      tau2 = sum(mCHE$sigma2)
      SE = clubSandwich::conf_int(mCHE, vcov = "CR2")[["SE"]]
      pi.lb.rve = clubSandwich::conf_int(mCHE, "CR2")[["beta"]] -
        crit * sqrt(tau2 + (SE^2))
      pi.ub.rve = clubSandwich::conf_int(mCHE, "CR2")[["beta"]] +
        crit * sqrt(tau2 + (SE^2))
      if (isTRUE(.raw.bin.es)){
        as.numeric(clubSandwich::conf_int(
          mCHE.g, "CR2")[["beta"]]) %>% 
          {ifelse(.==0, Inf, metapsyNNT(., cer))} -> nnt.g
      }
  
      mCHERes.RVE = with(mCHE, {
        data.frame(k = k.all,
                   g = as.numeric(
                     clubSandwich::conf_int(
                       mCHE, "CR2")[["beta"]]) %>%
                     ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                     round(round.digits),
                   g.ci = paste0(
                     "[",
                     clubSandwich::conf_int(
                       mCHE, "CR2")[["CI_L"]] %>%
                       ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                       round(round.digits), "; ",
                     clubSandwich::conf_int(
                       mCHE, "CR2")[["CI_U"]] %>%
                       ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                       round(round.digits), "]"),
                   p = clubSandwich::coef_test(
                     mCHE, "CR2")[["p_Satt"]] %>%
                     scales::pvalue(),
                   i2 = round(I2, 1),
                   i2.ci = "-",
                   prediction.ci = paste0(
                     "[", 
                     round(pi.lb.rve %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "; ",
                     round(pi.ub.rve %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "]"),
                   nnt = ifelse(identical(.type.es, "RR"), 
                                NA, metapsyNNT(as.numeric(
                                  clubSandwich::conf_int(
                                    mCHE, "CR2")[["beta"]]), nnt.cer)) %>% 
                     ifelse(isTRUE(.raw.bin.es), nnt.g, .) %>% 
                     round(round.digits) %>% abs(),
                   excluded = "none")
      })
      rownames(mCHERes.RVE) = "Three-Level Model (CHE)"
      mCHERes.RVE$excluded = 
        paste0("Number of clusters/studies: ", mCHE$s.nlevels[1],
               "; robust variance estimation (RVE) used.")
      mCHERes = mCHERes.RVE
  
      if ("threelevel.che" %in% which.run){
        message(crayon::green("DONE"))
      }
  
      if (length(multi.study) == 0 & "threelevel.che" %in% which.run){
        message("- ", crayon::yellow("[!] "), 
                "All included ES seem to be independent.",
                " A three-level CHE model is not adequate and ",
                "tau/I2 estimates are not trustworthy!")
        warn.end = TRUE
        which.run[!which.run == "threelevel.che"] -> which.run
      }
    }
  } else {
    mCHE = NULL
    mCHERes = NULL
  }


  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  RETURN                                                                   #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


  # Combine everything
  rbind(mGeneralRes, mLowestRes, 
        mHighestRes, mOutliersRes,
        mInfluenceRes, mRobRes, 
        mCombRes, mThreeLevelRes, mCHERes) -> summary
  
  # Fill NAs 
  summary[summary == "[NA; NA]"] = "[-; -]"
  summary[is.na(summary)] = "-"
  if (identical(.type.es, "RR")){
    colnames(summary)[2:3] = c("rr", "rr.ci")
  }

  # Define ES type for all models
  mGeneral$.type.es = .type.es
  mComb$.type.es = .type.es
  mLowest$.type.es = .type.es
  mHighest$.type.es = .type.es
  mOutliers$.type.es = .type.es
  mInfluence$.type.es = .type.es
  mRob$.type.es = .type.es
  mThreeLevel$.type.es = .type.es
  mCHE$.type.es = .type.es
  
  # Return
  list(summary = summary,
       model.overall = mGeneral,
       model.combined = mComb,
       model.lowest = mLowest,
       model.highest = mHighest,
       model.outliers = mOutliers,
       model.influence = mInfluence,
       model.rob = mRob,
       model.threelevel = mThreeLevel,
       model.threelevel.var.comp = mThreeLevel$variance.components,
       model.threelevel.che = mCHE,
       model.threelevel.che.var.comp = mCHE$variance.components,
       influence.analysis = influenceRes,
       which.run = which.run,
       data = data.original,
       html = html,
       round.digits = round.digits,
       nnt.cer = nnt.cer,
       use.rve = use.rve,
       .type.es = .type.es,
       .raw.bin.es = .raw.bin.es) -> returnlist
  

  class(returnlist) = c("runMetaAnalysis", "list")

  if ("lowest.highest" %in% which.run){
    message("- ", crayon::green("[OK] "), "Done!")
  }

  if (warn.end == TRUE){
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
#' @importFrom magrittr set_colnames
#' @importFrom dplyr as_tibble
#' @importFrom kableExtra kable_styling column_spec footnote
#' @importFrom crayon green blue magenta
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
    
    
    cat("Model results ")
    cat("------------------------------------------------ \n")
    dat = x$summary[run.models,1:8]
    print(dplyr::as_tibble(cbind(model = rownames(dat), dat)))
    
    # Add publication bias corrected results
    cat("\n")
    cat(paste0("Publication bias correction ('", 
               x$correctPublicationBias$which.run,
               "' model) "))
    cat("----------------------- \n")
    dat.pb = x$correctPublicationBias$summary[,1:8]
    print(dplyr::as_tibble(cbind(model = rownames(dat.pb), dat.pb)))
    
    
    if ("threelevel" %in% x$which.run){
      cat("\n")
      cat("Variance components (three-level model) ")
      cat("---------------------- \n")
      print(x$model.threelevel.var.comp)
    }
    
    if ("threelevel.che" %in% x$which.run){
      cat("\n")
      cat("Variance components (three-level CHE model) ")
      cat("------------------ \n")
      print(x$model.threelevel.che.var.comp)
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
        magrittr::set_colnames(c(".", "<i>k</i>", 
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
    
    cat("Model results ")
    cat("------------------------------------------------ \n")
    dat = x$summary[run.models, 1:8]
    print(dplyr::as_tibble(cbind(model = rownames(dat), dat)))
    
    if ("threelevel" %in% x$which.run){
      cat("\n")
      cat("Variance components (three-level model) ")
      cat("---------------------- \n")
      print(x$model.threelevel.var.comp)
    }
    
    if ("threelevel.che" %in% x$which.run){
      cat("\n")
      cat("Variance components (three-level CHE model) ")
      cat("------------------ \n")
      print(x$model.threelevel.che.var.comp)
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
        magrittr::set_colnames(c(".", "<i>k</i>", 
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
#' \code{"influence"}, \code{"baujat"}, \code{"loo-es"}, \code{"loo-i2"},
#' \code{"trimfill"}, \code{"limitmeta"} or \code{"selection"}.
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

plot.runMetaAnalysis = function(x, which = NULL, ...){

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
#' @importFrom magrittr set_colnames
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

  if (x$model.overall$fixed == TRUE & x$model.overall$random == FALSE){

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
          "'Combined' analysis: ES combined on study level assuming a correlation of rho =", paste0(x$model.combined$rho, ". \n"))}

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
          "'Combined' analysis: ES combined on study level assuming a correlation of rho =", paste0(x$model.combined$rho, ". \n"))}

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








