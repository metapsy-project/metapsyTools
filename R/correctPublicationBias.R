#' Correct the effect size for publication bias/small-study effects.
#'
#' This function allows to add effect sizes estimates corrected for publication bias/
#' small-study effects to results of the \code{runMetaAnalysis} function.
#'
#' @usage correctPublicationBias(model, 
#'                        which.run = model$which.run[1],
#'                        lower.is.better = TRUE,
#'                        selmodel.steps = c(0.025, 0.05),
#'                        ...)
#'
#' @param model An object of class \code{runMetaAnalysis}, created by the \code{runMetaAnalysis} function.
#' @param which.run The model in \code{model} that should be used for the publication bias analyses. 
#' Uses the default analysis in \code{model} if no value is specified by the user. Possible values are
#' \code{"overall"}, \code{"combined"}, \code{"lowest"}, \code{"highest"}, \code{"outliers"},
#' \code{"influence"} and \code{"rob"}.
#' @param lower.is.better Do lower values indicate better outcomes (i.e. higher effects)? Default is \code{TRUE}.
#' @param selmodel.steps Thresholds to be assumed for the step function in the selection model. Must be a vector
#' of numbers referring to the cut-points in the selection models. If two-sided testing is assumed for the
#' included studies, the cut-point must be doubled to obtain the assumed \emph{p}-value 
#' (e.g. \code{selmodel.steps = c(0.03, 0.05)} means that \emph{p}=0.06 and \emph{p}=0.10 are assumed as selection thresholds).
#' The default is \code{c(0.025, 0.05)}.
#' @param ... Additional arguments. See \link[meta]{trimfill.default} and \link[metasens]{limitmeta}.
#'
#' @return Returns an object of class \code{"runMetaAnalysis"} and \code{"correctPublicationBias"}. 
#' This object includes all original objects included in \code{model}, but adds a \code{list} object
#' with the name \code{correctPublicationBias}. This list object includes all three fitted publication
#' bias analysis models, as well as the generated results.
#'
#'
#' @examples
#' \dontrun{
#' # Run meta-analysis
#' data("psyCtrSubsetWide")
#' psyCtrSubsetWide %>%
#'   checkDataFormat() %>%
#'   checkConflicts() %>%
#'   expandMultiarmTrials() %>%
#'   calculateEffectSizes() %>%
#'   filterPoolingData(year > 2010) %>%
#'   runMetaAnalysis() -> res
#' 
#' # Correct for small-study-effects/publication bias
#' res %>% correctPublicationBias()
#' 
#' # Use additional arguments to control settings of the trim-and-fill
#' # and limit meta-analysis
#' correctPublicationBias(res,
#'                        which.run = "combined",
#'                        type = "R",
#'                        method.adjust = "mulim") 
#' 
#' # Generate plots
#' correctPublicationBias(res) %>% plot("trimfill")
#' correctPublicationBias(res) %>% plot("limitmeta")
#' correctPublicationBias(res) %>% plot("selection")
#' 
#' # Returned object is of class "runMetaAnalysis"; therefore,
#' # all S3 methods are available:
#' res %>% 
#'   correctPublicationBias() %>% 
#'   subgroupAnalysis(country)
#' }
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{runMetaAnalysis}}
#'
#' @details The \code{runMetaAnalysis} function is a wrapper running three meta-analytic 
#' methods to control the pooled effect size for publication bias and/or small-study effects:
#'
#' \itemize{
#'   \item \code{"trimfill"}. Applies Duval and Tweedie's (2000a, 2000b) trim-and-fill algorithm, using the 
#'   \link[meta]{trimfill} method in the \code{meta} package.
#'   \item \code{"limitmeta"}. Runs a limit meta-analysis as described in Rücker et al. (2011),
#'   using the implementation in the \code{limitmeta} package.
#'   \item \code{"selection"}. Runs a step function selection model using the \link[metafor]{selmodel}
#'   function in \code{metafor}. For details see e.g. Vevea and Hedges (1995).
#' }
#' 
#' @references 
#' Duval S & Tweedie R (2000a): A nonparametric "Trim and Fill" method of accounting for 
#' publication bias in meta-analysis. \emph{Journal of the American Statistical Association, 95}, 89–98
#' 
#' Duval S & Tweedie R (2000b): Trim and Fill: A simple funnel-plot-based method of 
#' testing and adjusting for publication bias in meta-analysis. \emph{Biometrics, 56}, 455–63
#' 
#' Rücker G, Schwarzer G, Carpenter JR, Binder H, Schumacher M (2011): Treatment-effect estimates adjusted for 
#' small-study effects via a limit meta-analysis. \emph{Biostatistics, 12}, 122–42
#' 
#' Vevea, J. L., & Hedges, L. V. (1995). A general linear model for estimating effect size in the 
#' presence of publication bias. \emph{Psychometrika, 60}(3), 419–435.  ⁠https://doi.org/10.1007/BF02294384⁠
#' 
#' @import dplyr numDeriv
#' @importFrom crayon green yellow cyan bold
#' @importFrom scales pvalue
#' @importFrom meta trimfill
#' @importFrom metasens limitmeta
#' @importFrom metafor selmodel rma.uni
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn argsAnywhere
#' @export correctPublicationBias


correctPublicationBias = function(model, 
                                  which.run = model$which.run[1],
                                  lower.is.better = TRUE,
                                  selmodel.steps = c(0.025, 0.05),
                                  ...){
  
  # Check class
  if (class(model)[1] != "runMetaAnalysis"){
    stop("Input must be of class 'runMetaAnalysis'. Did you apply 'runMetaAnalysis' first?")
  }
  
  
  # Throw error if which.run is not included
  if (which.run[1] == "combined" & model$model.combined$k == 1){
    stop("'which.run' is set to 'combined', but there is only k=1 study/ES.")}
  if (which.run[1] == "lowest" & model$model.lowest$k == 1){
    stop("'which.run' is set to 'lowest', but there is only k=1 study/ES.")}
  if (which.run[1] == "highest" & model$model.highest$k == 1){
    stop("'which.run' is set to 'highest', but there is only k=1 study/ES.")}
  if (which.run[1] == "influence" & model$model.influence$k == 1){
    stop("'which.run' is set to 'influence', but there is only k=1 study/ES.")}
  if (which.run[1] == "rob" & model$model.rob$k == 1){
    stop("'which.run' is set to 'rob', but there is only k=1 study/ES.")}
  if (!(which.run[1] %in% c("overall", "combined",
                          "outliers", "rob", "highest",
                          "lowest", "influence"))){
    stop("'which.run' must be 'overall', 'combined', 'outliers', 'rob',",
         " 'highest', 'lowest' or 'influence'.")
  }
  
  # Throw error if multilevel model is used
  if (which.run %in% c("threelevel", "threelevel.che")){
    stop("Adjustment cannot be performed for three-level models. ",
         "Change the model provided in 'which.run'")
  }
  
  # Get model type
  model.type = paste0("model.", which.run)
  M = model[[model.type]]
  message("- ", crayon::green("[OK] "), "'", model.type, 
          "' adjusted for small-study effects/publication bias.")
  
  # Get three dots
  dots = list(...)
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  1. Trim and Fill                                                         #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  tf.args = dots[names(dots) %in% names(as.list(argsAnywhere("trimfill.default")))]
  mTrimFill = do.call(meta::trimfill, c(x = list(M), tf.args))
  
  mTrimFillRes = with(mTrimFill,{
    
    data.frame(
      k = k,
      g = ifelse(fixed == TRUE & random == FALSE, TE.fixed, TE.random) %>% round(model$round.digits),
      g.ci = paste0("[", ifelse(fixed == TRUE & random == FALSE,
                                lower.fixed, lower.random) %>% round(model$round.digits),"; ",
                    ifelse(fixed == TRUE & random == FALSE,
                           upper.fixed, upper.random) %>% round(model$round.digits), "]"),
      p = ifelse(fixed == TRUE & random == FALSE, pval.fixed, pval.random) %>% scales::pvalue(),
      i2 = round(I2*100, 2),
      i2.ci = paste0("[", round(lower.I2*100, 2), "; ", round(upper.I2*100, 2), "]"),
      prediction.ci = paste0("[", round(lower.predict, model$round.digits), "; ",
                             round(upper.predict, model$round.digits), "]"),
      nnt = metapsyNNT(ifelse(fixed == TRUE & random == FALSE, TE.fixed, TE.random), model$nnt.cer) %>%
        round(model$round.digits) %>% abs(),
      excluded = "none"
    )
  })
  
  mTrimFillRes$excluded = paste(sum(mTrimFill[["trimfill"]]), 
                                "studies added")
  
  rownames(mTrimFillRes) = "Trim-and-fill method"
  
  message("- ", crayon::green("[OK] "),
          "Trim-and-fill method applied using the '", mTrimFill["type"], "' estimator.")
  
  
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  2. Limit Meta-Analysis                                                   #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  lm.args = dots[names(dots) %in% names(as.list(args(metasens::limitmeta)))]
  mLimit = do.call(metasens::limitmeta, c(x = list(M), lm.args))
  
  mLimitRes = with(mLimit,{
    
    sd.pi = (qt(0.975, k-1) * sqrt(seTE.adjust^2 + tau^2))
    lower.predict = TE.adjust - sd.pi
    upper.predict = TE.adjust + sd.pi
    
    data.frame(
      k = k,
      g = TE.adjust %>% round(model$round.digits),
      g.ci = paste0("[", lower.adjust %>% round(model$round.digits),"; ",
                    upper.adjust %>% round(model$round.digits), "]"),
      p = pval.adjust %>% scales::pvalue(),
      i2 = round(G.squared*100, 2),
      i2.ci = "-",
      prediction.ci = paste0("[", round(lower.predict, model$round.digits), "; ",
                             round(upper.predict, model$round.digits), "]"),
      nnt = metapsyNNT(TE.adjust, model$nnt.cer) %>%
        round(model$round.digits) %>% abs(),
      excluded = paste("For the limit meta-analysis, the value under I-squared refers",
                       "to the G-squared heterogeneity statistic.")
    )
  })
  
  rownames(mLimitRes) = "Limit meta-analysis"
  
  message("- ", crayon::green("[OK] "),
          "Shrunken estimate of the limit meta-analysis based on the '", 
          mLimit["method.adjust"], "' estimator.")
  
  
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  3. Three-Parameter Selection Model                                       #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  # Turn into rma.uni object
  M.metafor = with(M, {
    test = ifelse(hakn, "t", "z")
    
    if (is.null(exclude)){
      metafor::rma.uni(yi = TE, sei = seTE, slab = studlab,
                       method = method.tau, test = test)
    } else {
      metafor::rma.uni(yi = TE, sei = seTE, slab = studlab,
                       method = method.tau, test = test, 
                       subset = !exclude)
    }
  }) 
  
  mSelmodel = metafor::selmodel(M.metafor,
                                type = "stepfun",
                                alternative = 
                                  ifelse(lower.is.better, "less", "greater"),
                                steps = selmodel.steps) 
  
  # Extract information from overall model
  r.digits = model$round.digits
  m.nnt.cer = model$nnt.cer
  
  mSelmodelRes = with(mSelmodel,{
    
    data.frame(
      k = k,
      g = b %>% round(r.digits),
      g.ci = paste0(
        "[", ci.lb %>% round(r.digits),"; ",
        ci.ub %>% round(r.digits), "]"),
      p = pval %>% scales::pvalue(),
      i2 = round(M.metafor$I2, 2),
      i2.ci = "-",
      prediction.ci = paste0(
        "[", round(predict(mSelmodel)[["pi.lb"]], r.digits), "; ",
        round(predict(mSelmodel)[["pi.ub"]], r.digits), "]"),
      nnt = metapsyNNT(b[1], m.nnt.cer) %>%
        round(r.digits) %>% abs(),
      excluded = paste0(
        "Step-function selection model with cutpoints p=",
        paste(selmodel.steps*2, collapse = ", "), ". ",
        "The selection model parameter test was ",
        ifelse(LRTp < 0.05, "significant: ", "not significant: "),
        "\u03C7\u00B2","=", round(LRT, 3), " (p=", scales::pvalue(LRTp), "). ",
        "The model was fitted using maximum likelihood estimation."
      )
    )
  })
  
  rownames(mSelmodelRes) = "Selection model"
  
  message("- ", crayon::green("[OK] "),
          "Estimated effect using step function selection model (p=",
          paste(selmodel.steps*2, collapse = ", "), ").")
  
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  RETURN                                                                   #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  # Combine all results tables
  summary = rbind(mTrimFillRes, mLimitRes, mSelmodelRes)
  
  # Prepare for return
  returnlist = model
  list(summary = summary,
       model.trimfill = mTrimFill,
       model.limitmeta = mLimit,
       model.selection = mSelmodel,
       which.run = which.run) -> returnlist[["correctPublicationBias"]]
  
  class(returnlist) = c("runMetaAnalysis", "correctPublicationBias", "list")
  
  return(returnlist)
  
  
}

