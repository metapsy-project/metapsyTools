#' Correct the effect size for publication bias/small-study effects
#'
#' This function allows to add effect sizes estimates corrected for publication bias/
#' small-study effects to results of the \code{runMetaAnalysis} function.
#'
#' @usage correctPublicationBias(model, 
#'                        which.run = model$which.run[1],
#'                        lower.is.better = NULL,
#'                        selmodelSteps = 0.05,
#'                        ...)
#'
#' @param model An object of class \code{runMetaAnalysis}, created by the \code{runMetaAnalysis} function.
#' @param which.run The model in \code{model} that should be used for the publication bias analyses. 
#' Uses the first model in \code{model} if no value is specified by the user. Possible values are
#' \code{"overall"}, \code{"combined"}, \code{"lowest"}, \code{"highest"}, \code{"outliers"},
#' \code{"influence"} and \code{"rob"}.
#' @param lower.is.better Only used for the selection model. Should significant \emph{negative} (\code{TRUE}) or
#' significant \emph{positive} (\code{FALSE}) effects be assumed to be favored for publication? This refers to the sign of 
#' the effect size in the data, not to the scale of the outcome measure. If set to \code{NULL} (default), the direction 
#' is determined by the sign of the pooled effect in the selected model. A warning is issued if the 
#' provided value contradicts the sign of the pooled effect.
#' @param selmodelSteps Thresholds to be assumed for the step function in the selection model. Must be a vector
#' of numbers referring to the cut-points in the selection models. If two-sided testing is assumed for the
#' included studies, the cut-point must be doubled to obtain the assumed \emph{p}-value 
#' (e.g. \code{selmodelSteps = c(0.03, 0.05)} means that \emph{p}=0.06 and \emph{p}=0.10 are assumed as selection thresholds).
#' The default is \code{0.05}.
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
# Run meta-analysis
#' data("depressionPsyCtr")
#' depressionPsyCtr %>%
#'   checkDataFormat() %>%
#'   checkConflicts() %>%
#'   calculateEffectSizes() %>%
#'   filterPoolingData(condition_arm1 %in% c("cbt", "pst")) %>%
#'   runMetaAnalysis() -> res
#' 
#' # Correct for small-study-effects/publication bias
#' res %>% correctPublicationBias()
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
#' res <- res %>% correctPublicationBias()
#' metaRegression(res$model.overall, ~country)
#' }
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{runMetaAnalysis}}
#'
#' @details The \code{correctPublicationBias} function is a wrapper running three meta-analytic 
#' methods to control the pooled effect size for publication bias and/or small-study effects:
#'
#' \itemize{
#'   \item \code{"trimfill"}. Applies Duval and Tweedie's (2000a, 2000b) trim-and-fill algorithm, using the 
#'   \link[meta]{trimfill} method in the \code{meta} package.
#'   \item \code{"limitmeta"}. Runs a limit meta-analysis as described in Rücker et al. (2011),
#'   using the \link[metasens]{limitmeta} function in the \code{metasens} package.
#'   \item \code{"selection"}. Runs a step function selection model using the \link[metafor]{selmodel}
#'   function in \code{metafor}. For details see e.g. Vevea and Hedges (1995). The cut-points in 
#'   \code{selmodelSteps} refer to one-sided \emph{p}-values. Selection models are not estimated for 
#'   event rates (\code{es.measure = "EER"} or \code{"CER"}), for which selection on
#'   significance is not a meaningful assumption.
#' }
#' 
#' If the selected model is a fixed-effect ("common effect") model (\code{method.tau = "FE"}), the trim-and-fill and
#' selection model are also fitted assuming a fixed effect. The limit meta-analysis is 
#' based on a random-effects model by definition.
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
#' @importFrom metafor selmodel rma.uni confint.rma.uni.selmodel predict.rma
#' @importFrom stats dffits model.matrix rnorm rstudent quantile plogis qt sd
#' @importFrom utils combn argsAnywhere
#' @export correctPublicationBias


correctPublicationBias = function(model, 
                                  which.run = model$which.run[1],
                                  lower.is.better = NULL,
                                  selmodelSteps = 0.05,
                                  ...){
  
  # Check class
  if (class(model)[1] != "runMetaAnalysis"){
    stop("Input must be of class 'runMetaAnalysis'. Did you apply 'runMetaAnalysis' first?")
  }
  
  which.run = which.run[1]
  
  # 'lowest.highest' are two models; one of them must be selected
  if (identical(which.run, "lowest.highest")){
    stop("'which.run' is 'lowest.highest', which includes two models.",
         " Please set 'which.run' to 'lowest' or 'highest'.")
  }
  
  # Throw an error if which.run does not contain required model
  if (!(which.run %in% c("overall", "combined",
                         "outliers", "rob", "highest",
                         "lowest", "influence"))){
    stop("'which.run' must be 'overall', 'combined', 'outliers', 'rob',",
         " 'highest', 'lowest' or 'influence'.")
  }
  
  models.fitted = model$which.run
  if ("lowest.highest" %in% models.fitted){
    models.fitted = c(models.fitted, "lowest", "highest")
  }
  if (!(which.run %in% models.fitted)){
    stop("'", which.run, "' model has not been fitted using runMetaAnalysis.")
  }
  
  if (which.run %in% unlist(model$error.model.list)){
    stop("There was an error fitting the '", which.run, "' model using runMetaAnalysis;",
         " the correctPublicationBias function cannot be applied.")
  }
  
  # Get model type
  model.type = paste0("model.", which.run)
  M = model[[model.type]]
  if (is.null(M)){
    stop("The '", which.run, "' model is not available in 'model'.")
  }
  
  # Throw error if which.run only includes one study
  if (M$k == 1){
    stop("'which.run' is set to '", which.run, "', but there is only k=1 study/ES.")
  }
  
  message("- ", crayon::green("[OK] "), "'", model.type, 
          "' adjusted for small-study effects/publication bias.")
  
  # Get three dots
  dots = list(...)
  
  # Get steering variables
  .raw.bin.es = model[[".raw.bin.es"]]
  .type.es = model[[".type.es"]]
  nntCer = model[["nntCer"]]
  round.digits = model[["round.digits"]]
  M.orig = M
  is.common = isTRUE(M$common) && !isTRUE(M$random)
  is.ratio = identical(.type.es, "RR") || identical(.type.es, "ROM")
  is.event.rate = identical(.type.es, "EER") || identical(.type.es, "CER")
  
  # Functions to back-transform and present the estimates
  back.transform = function(x){
    if (is.ratio) { exp(x) } else if (is.event.rate) { plogis(x) } else { x }
  }
  format.es = function(x){
    if (is.null(x)) return(NA)
    round(back.transform(as.numeric(x)), round.digits)
  }
  format.ci = function(lower, upper){
    paste0("[", format.es(lower), "; ", format.es(upper), "]")
  }
  format.p = function(p){
    if (is.event.rate) { "-" } else { scales::pvalue(p) }
  }
  # NNTs are based on the same effect metric as in the original model
  calc.nnt = function(te){
    if (is.event.rate || identical(.type.es, "ROM")) return(NA)
    if (identical(.type.es, "RR")){
      return(metapsyNNT(RR = exp(as.numeric(te)), CER = nntCer))
    }
    metapsyNNT(abs(as.numeric(te)), nntCer)
  }
  failed.row = function(note){
    data.frame(k = "-", g = "-", g.ci = "[-; -]", p = "-", i2 = "-", i2.ci = "-",
               prediction.ci = "[-; -]", nnt = "-", excluded = note)
  }
  
  # Risk difference model; used to calculate NNTs from raw event data
  if (isTRUE(.raw.bin.es)) {
    M.rd = updateMeta(M, sm = "RD")
    if (sum(M.rd$exclude) == 0){
      M.rd$exclude = NULL
    }
  }
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  1. Trim and Fill                                                         #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  M.tf = M
  if (sum(M.tf$exclude) == 0){
    M.tf$exclude = NULL
  }
  tf.args = dots[names(dots) %in% 
                   names(as.list(argsAnywhere("trimfill.default")))]
  # Pool using the same type of model as in 'model'
  if (!"common" %in% names(tf.args)) tf.args$common = is.common
  if (!"random" %in% names(tf.args)) tf.args$random = !is.common
  
  mTrimFill = tryCatch2(
    do.call(meta::trimfill, c(x = list(M.tf), tf.args)))
  
  if (mTrimFill$has.error){
    message("- ", crayon::yellow("[!] "),
            "Trim-and-fill: ", mTrimFill$error$message)
    mTrimFillRes = failed.row("Trim-and-fill method could not be applied.")
    
  } else {
    
    if (isTRUE(.raw.bin.es)) {
      with(do.call(meta::trimfill, 
                   c(x = list(M.rd), tf.args)), {
                     ifelse(isTRUE(common) & !isTRUE(random),
                            abs(TE.common)^-1, abs(TE.random)^-1)
                   }) -> nnt.raw.bin.es
    }
    
    mTrimFillRes = with(mTrimFill$value,{
      use.common = isTRUE(common) & !isTRUE(random)
      data.frame(
        k = k,
        g = format.es(ifelse(use.common, TE.common, TE.random)),
        g.ci = format.ci(ifelse(use.common, lower.common, lower.random),
                         ifelse(use.common, upper.common, upper.random)),
        p = format.p(ifelse(use.common, pval.common, pval.random)),
        i2 = round(I2*100, 2),
        i2.ci = paste0("[", round(lower.I2*100, 2), "; ", 
                       round(upper.I2*100, 2), "]"),
        prediction.ci = format.ci(lower.predict, upper.predict),
        nnt = calc.nnt(ifelse(use.common, TE.common, TE.random)) %>% 
          ifelse(isTRUE(.raw.bin.es), nnt.raw.bin.es, .) %>% 
          round(round.digits) %>% abs(),
        excluded = "none")
    })
    
    mTrimFillRes$excluded = paste(sum(mTrimFill$value[["trimfill"]]), 
                                  "studies added.")
    
    message("- ", crayon::green("[OK] "),
            "Trim-and-fill method applied using the '", 
            mTrimFill$value["type"], "' estimator.")
    
    if (identical(.type.es, "RR") && isTRUE(.raw.bin.es) && mTrimFillRes$g == 1){
      mTrimFillRes$nnt = Inf
    }
  }
  rownames(mTrimFillRes) = "Trim-and-fill method"
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  2. Limit Meta-Analysis                                                   #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  lm.args = dots[names(dots) %in% names(as.list(args(metasens::limitmeta)))]
  
  mLimit = tryCatch2(
    do.call(metasens::limitmeta, c(x = list(M), lm.args)))
  
  if (mLimit$has.error){
    message("- ", crayon::yellow("[!] "),
            "Limit meta-analysis: ", mLimit$error$message)
    mLimitRes = failed.row("Limit meta-analysis method could not be applied.")
    
  } else {
    
    # Use bootstrapping if selected
    if (isTRUE(model$i2.ci.boot)){
      message(
        crayon::cyan(
          crayon::bold(
            "- Parametric bootstrap for g2 confidence intervals (limit meta-analysis) ...")))
      
      keep = if (is.null(M$exclude)) rep(TRUE, length(M$TE)) else !M$exclude
      boot.se = M$seTE[keep]
      boot.mu = ifelse(is.common, M$TE.common, M$TE.random)
      boot.tau2 = ifelse(is.common, 0, M$tau2)
      boot.method.tau = M$method.tau
      
      # Sample
      g2.vec = rep(NA_real_, model$nsim.boot)
      tau2.vec = rep(NA_real_, model$nsim.boot)
      counter = 0
      for (i in seq_len(model$nsim.boot)){
        tmp = try({
          sim.te = rnorm(length(boot.se), boot.mu, sqrt(boot.tau2 + boot.se^2))
          M.boot = meta::metagen(TE = sim.te, seTE = boot.se, 
                                 method.tau = boot.method.tau)
          res.boot = do.call(metasens::limitmeta, c(x = list(M.boot), lm.args))
          c("G.squared" = res.boot$G.squared, "tau2" = res.boot$tau^2)
        }, silent = TRUE)
        counter = counter + 1
        if (!inherits(tmp, "try-error")) {
          g2.vec[i] = tmp[["G.squared"]]
          tau2.vec[i] = tmp[["tau2"]]
        }
        if (counter %in% seq(0, model$nsim.boot, model$nsim.boot/100)){
          cat(crayon::green(
            paste0((counter/model$nsim.boot)*100, "% completed | ")))
        }
      }
      cat(crayon::green("DONE \n"))
      
      mLimit$value$G.squared.ci = quantile(g2.vec, c(.025, .975), na.rm = TRUE)
      mLimit$value$se.tau2 = sd(tau2.vec, na.rm = TRUE)
    }
    
    if (isTRUE(.raw.bin.es)) {
      with(do.call(metasens::limitmeta, 
                   c(x = list(M.rd), lm.args)), {
                     abs(TE.adjust)^-1}) -> nnt.raw.bin.es
    }
    
    mLimitRes = with(mLimit$value,{
      sd.pi = (qt(0.975, k-1) * sqrt(seTE.adjust^2 + tau^2))
      data.frame(
        k = k,
        g = format.es(TE.adjust),
        g.ci = format.ci(lower.adjust, upper.adjust),
        p = format.p(pval.adjust),
        i2 = round(G.squared*100, 2),
        i2.ci = ifelse(isTRUE(model$i2.ci.boot),
                       paste0(
                         "[", 
                         round(G.squared.ci[1]*100, round.digits), "; ",
                         round(G.squared.ci[2]*100, round.digits), "]"), "-"),
        prediction.ci = format.ci(TE.adjust - sd.pi, TE.adjust + sd.pi),
        nnt = calc.nnt(TE.adjust) %>% 
          ifelse(isTRUE(.raw.bin.es), nnt.raw.bin.es, .) %>% 
          round(round.digits) %>% abs(),
        excluded = paste(
          "For the limit meta-analysis, the value under I-squared refers",
          "to the G-squared heterogeneity statistic.")
      )
    })
    
    if (identical(.type.es, "RR") && isTRUE(.raw.bin.es) && mLimitRes$g == 1){
      mLimitRes$nnt = Inf
    }
    
    message("- ", crayon::green("[OK] "),
            "Shrunken estimate of the limit meta-analysis based on the '", 
            mLimit$value["method.adjust"], "' estimator.")
  }
  rownames(mLimitRes) = "Limit meta-analysis"
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  3. Three-Parameter Selection Model                                       #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  mSelmodel = NULL
  r.digits = model$round.digits
  m.nntCer = model$nntCer
  
  # A selection model assumes selection on the significance of an effect relative
  # to the null value. Event rates have no meaningful null value (a logit of 0 would 
  # imply an event rate of 50%), so the model is not estimated for them.
  if (is.event.rate){
    
    message("- ", crayon::yellow("[!] "),
            "Selection model not applicable to event rates; it was not estimated.")
    mSelmodelRes = failed.row(
      paste("Selection models assume selection on significance;",
            "this is not meaningful for pooled event rates."))
    
  } else {
    
    # Turn into rma.uni object
    test = ifelse(identical(M$method.random.ci[1], "HK"), "t", "z")
    method.sm = ifelse(is.common, "FE", M$method.tau)
    has.exclude = !is.null(M$exclude) && sum(M$exclude) > 0
    keep = if (has.exclude) !M$exclude else rep(TRUE, length(M$TE))
    M.metafor = metafor::rma.uni(yi = M$TE[keep], sei = M$seTE[keep], 
                                 slab = M$studlab[keep],
                                 method = method.sm, test = test)
    
    # Direction of the favored effects: by default, significant effects with the 
    # same sign as the pooled effect are assumed to be more likely to be published.
    pooled.es = ifelse(is.common, M$TE.common, M$TE.random)
    if (is.null(lower.is.better[1])){
      alternative = ifelse(pooled.es < 0, "less", "greater")
      message("- ", crayon::green("[OK] "),
              "Selection model: effects ", 
              ifelse(pooled.es < 0, "below", "above"), " the null value",
              " assumed to be favored (based on the sign of the pooled effect).")
    } else {
      alternative = ifelse(isTRUE(lower.is.better[1]), "less", "greater")
      if ((pooled.es > 0 && alternative == "less") || 
          (pooled.es < 0 && alternative == "greater")){
        warning("Selection model: 'lower.is.better = ", isTRUE(lower.is.better[1]), 
                "' assumes that ", ifelse(alternative == "less", "negative", "positive"), 
                " effects are favored, but the pooled effect is ",
                ifelse(alternative == "less", "positive", "negative"), 
                ". If effects in your data are coded differently, the selection model",
                " is estimated in the wrong direction and can be severely biased.",
                " Check the coding of your effect sizes, or use the default",
                " 'lower.is.better = NULL'.", call. = FALSE)
      }
    }
    
    mSelmodel = try({
      metafor::selmodel(M.metafor, type = "stepfun", alternative = alternative,
                        steps = selmodelSteps)}, silent = TRUE)
    
    if (inherits(mSelmodel, "try-error")){
      message("- ", crayon::yellow("[!] "),
              "Selection model: problems detected! ",
              as.character(attr(mSelmodel, "condition")$message))
      mSelmodelRes = failed.row("Selection model could not be calculated.")
      
    } else {
      
      # Re-run using Hedges g and calculate CER for raw binary ES
      if (isTRUE(.raw.bin.es)){
        meta::metaprop(
          event = round(M$data$.event.c[keep]),
          n = round(M$data$.n.c[keep]),
          common = is.common,
          random = !is.common) %>%
          {ifelse(is.common, .$TE.common, .$TE.random)} %>% 
          plogis() -> cer
        
        with(M$data[keep,], {
          data.frame(.event.e, .event.c, .n.e, .n.c)
        }) %>% 
          apply(., 1, function(x){
            esc::esc_2x2(grp1yes = as.numeric(x[".event.e"]), 
                         grp1no = as.numeric(x[".n.e"]) - as.numeric(x[".event.e"]),
                         grp2yes = as.numeric(x[".event.c"]),
                         grp2no = as.numeric(x[".n.c"]) - as.numeric(x[".event.c"]),
                         es.type = "g") %>% 
              {data.frame(es = .$es, se = .$se)}
          }) %>% 
          do.call(rbind, .) -> data.g
        
        metafor::rma.uni(yi = es, sei = se, data = data.g,
                         method = method.sm, test = test) %>% 
          metafor::selmodel(type = "stepfun", alternative = alternative,
                            steps = selmodelSteps) %>% 
          {metapsyNNT(abs(as.numeric(.$b)), CER = cer)} -> nnt.g
      }
      
      # Calculate profile-likelihood CIs (not available for fixed-effect models)
      tau2.sm = tryCatch(
        metafor::confint.rma.uni.selmodel(mSelmodel)[[1]]$random["tau^2", 1:3],
        error = function(e) c(mSelmodel$tau2, NA, NA))
      i2.sm = i2.gen(tau2.sm, mSelmodel$k, mSelmodel$vi)
      pred.sm = predict.rma(mSelmodel)
      
      mSelmodelRes = with(mSelmodel,{
        data.frame(
          k = k,
          g = format.es(as.numeric(b)),
          g.ci = format.ci(ci.lb, ci.ub),
          p = format.p(pval),
          i2 = round(i2.sm[1], 2),
          i2.ci = paste0("[", round(i2.sm[2], r.digits), "; ",
                         round(i2.sm[3], r.digits), "]"),
          prediction.ci = format.ci(pred.sm[["pi.lb"]], pred.sm[["pi.ub"]]),
          nnt = calc.nnt(as.numeric(b)) %>% 
            ifelse(isTRUE(.raw.bin.es), nnt.g, .) %>% 
            round(round.digits) %>% abs(),
          excluded = paste0(
            "Step-function selection model with cutpoints p=",
            paste(selmodelSteps*2, collapse = ", "), 
            " (two-sided; ", alternative, " than the null value favored). ",
            "The selection model parameter test was ",
            ifelse(LRTp < 0.05, "significant: ", "not significant: "),
            "\u03C7\u00B2","=", round(LRT, 3), 
            " (p=", scales::pvalue(LRTp), "). ",
            "The model was fitted using maximum likelihood estimation."
          )
        )
      })
      
      message("- ", crayon::green("[OK] "),
              "Estimated effect using step function selection model (p=",
              paste(selmodelSteps*2, collapse = ", "), ").")
    }
  }
  
  rownames(mSelmodelRes) = "Selection model"
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  RETURN                                                                   #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  # Combine all results tables
  summary = rbind(mTrimFillRes, mLimitRes, mSelmodelRes)
  # Fill NAs 
  summary[summary == "[NA; NA]"] = "[-; -]"
  summary[is.na(summary)] = "-"
  if (identical(.type.es, "RR")){
    colnames(summary)[2:3] = c("rr", "rr.ci")
  }
  if (identical(.type.es, "ROM")){
    colnames(summary)[2:3] = c("rom", "rom.ci")
    summary$nnt = "-"
  }
  if (identical(.type.es, "EER")){
    colnames(summary)[2:3] = c("eer", "eer.ci")
  }
  if (identical(.type.es, "CER")){
    colnames(summary)[2:3] = c("cer", "cer.ci")
  }
  
  # Prepare for return
  returnlist = model
  list(summary = summary,
       model.trimfill = mTrimFill$value,
       model.limitmeta = mLimit$value,
       model.selection = mSelmodel,
       which.run = which.run,
       call = match.call()) -> returnlist[["correctPublicationBias"]]
  
  class(returnlist) = c("runMetaAnalysis", "correctPublicationBias", "list")
  
  return(returnlist)
  
}
