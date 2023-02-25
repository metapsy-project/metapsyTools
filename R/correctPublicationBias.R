#' Correct the effect size for publication bias/small-study effects
#'
#' This function allows to add effect sizes estimates corrected for publication bias/
#' small-study effects to results of the \code{runMetaAnalysis} function.
#'
#' @usage correctPublicationBias(model, 
#'                        which.run = model$which.run[1],
#'                        lower.is.better = TRUE,
#'                        selmodel.steps = 0.05,
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
#' res %>% 
#'   correctPublicationBias() %>% 
#'   metaRegression(country)
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
#' @importFrom metafor selmodel rma.uni confint.rma.uni.selmodel
#' @importFrom stats dffits model.matrix rnorm rstudent quantile
#' @importFrom utils combn argsAnywhere
#' @export correctPublicationBias


correctPublicationBias = function(model, 
                                  which.run = model$which.run[1],
                                  lower.is.better = TRUE,
                                  selmodel.steps = 0.05,
                                  ...){
  
  # Check class
  if (class(model)[1] != "runMetaAnalysis"){
    stop("Input must be of class 'runMetaAnalysis'. Did you apply 'runMetaAnalysis' first?")
  }
  
  # Throw an error if which.run does not contain required model
  if (!(which.run[1] %in% c("overall", "combined",
                            "outliers", "rob", "highest",
                            "lowest", "influence"))){
    stop("'which.run' must be 'overall', 'combined', 'outliers', 'rob',",
         " 'highest', 'lowest' or 'influence'.")
  }
  
  if (!(which.run[1] %in% model$which.run)){
    stop("'", which.run[1], "' model has not been fitted using runMetaAnalysis.")
  }
  
  if (which.run[1] %in% unlist(model$error.model.list)){
    stop("There was an error fitting the '", which.run[1], "' model using runMetaAnalysis;",
         " the correctPublicationBias function cannot be applied.")
  }
  
  # Throw error if which.run is not included
  if (which.run[1] == "combined"){
    if (model$model.combined$k == 1){
      stop("'which.run' is set to 'combined', but there is only k=1 study/ES.")
    }
  }
  if (which.run[1] == "lowest"){
    if (model$model.lowest$k == 1){
      stop("'which.run' is set to 'lowest', but there is only k=1 study/ES.") 
    }
  }
  if (which.run[1] == "highest"){
    if (model$model.highest$k == 1){
      stop("'which.run' is set to 'highest', but there is only k=1 study/ES.")
    }
  }
  if (which.run[1] == "influence" ){
    if (model$model.influence$k == 1){
      stop("'which.run' is set to 'influence', but there is only k=1 study/ES.")
    }
  }
  if (which.run[1] == "rob" ){
    if (model$model.rob$k == 1){
      stop("'which.run' is set to 'rob', but there is only k=1 study/ES.")
    }
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
  
  # Get steering variables
  .raw.bin.es = model[[".raw.bin.es"]]
  .type.es = model[[".type.es"]]
  nnt.cer = model[["nnt.cer"]]
  round.digits = model[["round.digits"]]
  M.orig = M
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  1. Trim and Fill                                                         #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  if (sum(M$exclude) == 0){
    M$exclude = NULL
  }
  tf.args = dots[names(dots) %in% 
                   names(as.list(argsAnywhere("trimfill.default")))]
  
  mTrimFill = tryCatch2(
    do.call(meta::trimfill, c(x = list(M), tf.args)))
  
  if (is.null(mTrimFill$value) & 
      is.null(mTrimFill$warning)){
    message("- ", crayon::yellow("[!] "),
            "Trim-and-fill: ", mTrimFill$error$message)
    
    mTrimFillRes = data.frame(
      k = "-",
      g = "-",
      g.ci = paste0("[", "-","; ", "-", "]"),
      p = "-",
      i2 = "-",
      i2.ci = "-",
      prediction.ci = paste0(
        "[", "-", "; ", "-", "]"),
      nnt = "-",
      excluded = "Trim-and-fill method could not be applied.")
    rownames(mTrimFillRes) = "Trim-and-fill method"
    
  } else {
    
    if (isTRUE(.raw.bin.es)) {
      M.rd = meta::update.meta(M, sm="RD")
      M.orig.rd = meta::update.meta(M.orig, sm="RD")
      if (sum(M.rd$exclude) == 0){
        M.rd$exclude = NULL
      }
      with(do.call(meta::trimfill, 
                   c(x = list(M.rd), tf.args)), {
                     ifelse(isTRUE(common) & isTRUE(random),
                            abs(TE.common)^-1, abs(TE.random)^-1)
                   }) -> nnt.raw.bin.es
    }
    
    mTrimFillRes = with(mTrimFill$value,{
      data.frame(
        k = k,
        g = ifelse(
          isTRUE(common) & isTRUE(random), 
          TE.common, TE.random) %>% 
          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
          round(round.digits),
        g.ci = paste0("[", 
                      ifelse(isTRUE(common) & isTRUE(random),
                             lower.common, lower.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "; ",
                      ifelse(isTRUE(common) & isTRUE(random),
                             upper.common, upper.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "]"),
        p = ifelse(isTRUE(common) & isTRUE(random), 
                   pval.common, pval.random) %>% 
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
          ifelse(isTRUE(common) & isTRUE(random), 
                 TE.common, TE.random), nnt.cer) %>%
          ifelse(identical(.type.es, "RR"), NA, .) %>% 
          ifelse(isTRUE(.raw.bin.es), 
                 nnt.raw.bin.es, .) %>% 
          round(round.digits) %>% abs(),
        excluded = "none")
    })
    
    mTrimFillRes$excluded = paste(sum(mTrimFill$value[["trimfill"]]), 
                                  "studies added.")
    
    rownames(mTrimFillRes) = "Trim-and-fill method"
    
    message("- ", crayon::green("[OK] "),
            "Trim-and-fill method applied using the '", 
            mTrimFill$value["type"], "' estimator.")
    
    if (isTRUE(.raw.bin.es)){
      M.rd = M.orig.rd
    }
    
    if (identical(.type.es, "RR")){
      if (mTrimFillRes$g == 1 & isTRUE(.raw.bin.es)){
        mTrimFillRes$nnt = Inf
      }
    }
  }
  
  M = M.orig
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  2. Limit Meta-Analysis                                                   #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  lm.args = dots[names(dots) %in% names(as.list(args(metasens::limitmeta)))]
  
  mLimit = tryCatch2(
    do.call(metasens::limitmeta, c(x = list(M), lm.args)))
  
  # Use bootstrapping if selected
  if (model$i2.ci.boot){
    message(
      crayon::cyan(
        crayon::bold(
          "- Parametric bootstrap for g2 confidence intervals (limit meta-analysis) ...")))
    
    dat = list(dat = M$data, lm.args = lm.args, call = M$call)
    # Bootstrap worker
    boot.func = function(data.boot) {
      data.bt = data.boot$data
      lm.args = data.boot$lm.args
      call = data.boot$call
      call$TE = data.bt$.g
      M.update = eval(call)
      res = do.call(metasens::limitmeta, c(x = list(M.update), lm.args))
      return(c("G.squared" = res$G.squared, "tau2" = res$tau^2))
    }
    # Generate data
    data.gen = function(dat, mle=list(mu=M$TE.random, tau2=M$tau2)) {
      dat.bt = dat$dat
      dat.bt$.g = rnorm(nrow(dat.bt), mle$mu, sqrt(mle$tau2 + dat.bt$.g_se^2))
      dat.bt$.g_se = dat$.g_se
      dat.ret = list(data = dat.bt, lm.args = dat$lm.args, call = dat$call)
      return(dat.ret)
    }
    # Sample
    g2.vec = numeric(model$nsim.boot)
    tau2.vec = numeric(model$nsim.boot)
    counter = 0
    for (i in 1:model$nsim.boot){
      tmp = try({boot.func(data.gen(dat))}, silent=TRUE)
      if (inherits(tmp, "try-error")) {
        g2.vec[i] = NA
        tau2.vec[i] = NA
        counter = counter + 1
      } else {
        g2.vec[i] = tmp[["G.squared"]]
        tau2.vec[i] = tmp[["tau2"]]
        counter = counter + 1
        if (counter %in% seq(0, model$nsim.boot, model$nsim.boot/100)){
          cat(crayon::green(
            paste0((counter/model$nsim.boot)*100, "% completed | ")))
        }
        if (identical(counter, model$nsim.boot)){
          cat(crayon::green("DONE \n"))
        }
      }
    }
    
    sav.g2 = g2.vec[!is.na(g2.vec)]
    sav.tau2 = tau2.vec[!is.na(tau2.vec)]
    mLimit$value$G.squared.ci = quantile(sav.g2, c(.025, .975))
    mLimit$value$se.tau2 = sd(sav.tau2)
  }
  
  if (is.null(mLimit$value) & 
      is.null(mLimit$warning)){
    message("- ", crayon::yellow("[!] "),
            "Limit meta-analysis: ", mLimit$error$message)
    
    mLimitRes = data.frame(
      k = "-",
      g = "-",
      g.ci = paste0("[", "-","; ", "-", "]"),
      p = "-",
      i2 = "-",
      i2.ci = "-",
      prediction.ci = paste0(
        "[", "-", "; ", "-", "]"),
      nnt = "-",
      excluded = "Limit meta-analysis method could not be applied.")
    rownames(mLimitRes) = "Limit meta-analysis"
    
  } else {
    
    if (isTRUE(.raw.bin.es)) {
      with(do.call(metasens::limitmeta, 
                   c(x = list(M.rd), lm.args)), {
                     abs(TE.adjust)^-1}) -> nnt.raw.bin.es
    }
    
    mLimitRes = with(mLimit$value,{
      
      sd.pi = (qt(0.975, k-1) * sqrt(seTE.adjust^2 + tau^2))
      lower.predict = TE.adjust - sd.pi
      upper.predict = TE.adjust + sd.pi
      
      data.frame(
        k = k,
        g = TE.adjust %>% 
          ifelse(identical(.type.es, "RR"), 
                 exp(.), .) %>% 
          round(round.digits),
        g.ci = paste0("[", 
                      lower.adjust %>% 
                        ifelse(identical(.type.es, "RR"), 
                               exp(.), .) %>% 
                        round(round.digits),"; ",
                      upper.adjust %>% 
                        ifelse(identical(.type.es, "RR"), 
                               exp(.), .) %>% 
                        round(round.digits), "]"),
        p = pval.adjust %>% scales::pvalue(),
        i2 = round(G.squared*100, 2),
        i2.ci = ifelse(model$i2.ci.boot,
                       paste0(
                         "[", 
                         round(G.squared.ci[1]*100, round.digits), "; ",
                         round(G.squared.ci[2]*100, round.digits), "]"), "-"),
        prediction.ci = paste0(
          "[", 
          round(lower.predict %>% 
                  ifelse(identical(.type.es, "RR"), 
                         exp(.), .), 
                round.digits), "; ",
          round(upper.predict %>% 
                  ifelse(identical(.type.es, "RR"), 
                         exp(.), .), 
                round.digits), "]"),
        nnt = metapsyNNT(
          TE.adjust, nnt.cer) %>%
          ifelse(identical(.type.es, "RR"), NA, .) %>% 
          ifelse(isTRUE(.raw.bin.es), 
                 nnt.raw.bin.es, .) %>% 
          round(round.digits) %>% abs(),
        excluded = paste(
          "For the limit meta-analysis, the value under I-squared refers",
          "to the G-squared heterogeneity statistic.")
      )
    })
    
    if (identical(.type.es, "RR")){
      if (mLimitRes$g == 1 & isTRUE(.raw.bin.es)){
        mLimitRes$nnt = Inf
      }
    }
    rownames(mLimitRes) = "Limit meta-analysis"
    message("- ", crayon::green("[OK] "),
            "Shrunken estimate of the limit meta-analysis based on the '", 
            mLimit$value["method.adjust"], "' estimator.")
    }
  
  
  
  
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
  
  mSelmodel = try({metafor::selmodel(M.metafor,
                   type = "stepfun",
                   alternative = 
                     ifelse(lower.is.better, "less", "greater"),
                   steps = selmodel.steps)}, silent = TRUE)
  .selmodel.error = FALSE
  
  if (inherits(mSelmodel, "try-error")){
    message("- ", crayon::yellow("[!] "),
            "Selection model: problems detected!")
    mSelmodel = try(metafor::selmodel(M.metafor,
                  type = "stepfun",
                  alternative = 
                    ifelse(lower.is.better, "less", "greater"),
                  steps = selmodel.steps, 
                  verbose = TRUE))
    .selmodel.error = TRUE
  }
  
  # Extract information from overall model
  r.digits = model$round.digits
  m.nnt.cer = model$nnt.cer
  
  if (isTRUE(.selmodel.error)){
    mSelmodelRes = data.frame(
                    k = "-",
                    g = "-",
                    g.ci = paste0("[", "-","; ", "-", "]"),
                    p = "-",
                    i2 = "-",
                    i2.ci = "-",
                    prediction.ci = paste0(
                      "[", "-", "; ", "-", "]"),
                    nnt = "-",
                    excluded = "Selection model could not be calculated.")
  } else {
    # Re-run using Hedges g and calculate CER for raw binary ES
    if (isTRUE(.raw.bin.es)){
      with(M, {
        meta::metaprop(
          event = round(data$.event.c),
          n = round(data$.n.c),
          common = ifelse(method.tau == "FE", 
                         TRUE, FALSE),
          random = ifelse(method.tau == "FE", 
                          FALSE, TRUE)) %>%
          {ifelse(isTRUE(common) & isTRUE(random), 
                  .$TE.common, .$TE.random)} %>% 
          {exp(.)/(1+exp(.))}
      }) -> cer    
      
      with(M$data, {
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
      
      test = ifelse(M$hakn, "t", "z")
      exclude = M$exclude
      
      if (is.null(exclude)){
        metafor::rma.uni(yi = es, sei = se, data = data.g,
                         method = M$method.tau, test = test) %>% 
          metafor::selmodel(type = "stepfun",
                            alternative = 
                              ifelse(lower.is.better, "less", "greater"),
                            steps = selmodel.steps) %>% 
          {metapsyNNT(as.numeric(.$b), 
                      CER = cer)} -> nnt.g
      } else {
        metafor::rma.uni(yi = es, sei = se, 
                         data = data.g,
                         method = M$method.tau, test = test, 
                         subset = !exclude) %>% 
          metafor::selmodel(type = "stepfun",
                            alternative = 
                              ifelse(lower.is.better, "less", "greater"),
                            steps = selmodel.steps) %>% 
          {metapsyNNT(as.numeric(.$b), 
                      CER = cer)} -> nnt.g
      }
    }

    # Calculate profile-likelihood CIs
    tau2.sm = metafor::confint.rma.uni.selmodel(mSelmodel)[[1]]$random["tau^2",1:3]
    mSelmodel$i2 = i2.gen(tau2.sm, mSelmodel$k, mSelmodel$vi)[1]
    mSelmodel$i2.lower = i2.gen(tau2.sm, mSelmodel$k, mSelmodel$vi)[2]
    mSelmodel$i2.upper = i2.gen(tau2.sm, mSelmodel$k, mSelmodel$vi)[3]
    
    mSelmodelRes = with(mSelmodel,{
      data.frame(
        k = k,
        g = b %>% 
          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
          round(r.digits),
        g.ci = paste0(
          "[", 
          ci.lb %>% 
            ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
            round(r.digits),"; ",
          ci.ub %>% 
            ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
            round(r.digits), "]"),
        p = pval %>% scales::pvalue(),
        i2 = round(i2, 2),
        i2.ci = paste0(
          "[", 
          round(i2.lower, r.digits), "; ",
          round(i2.upper, r.digits), "]"),
        prediction.ci = paste0(
          "[", 
          round(predict(mSelmodel)[["pi.lb"]] %>% 
                  ifelse(identical(.type.es, "RR"), exp(.), .), 
                r.digits), "; ",
          round(predict(mSelmodel)[["pi.ub"]] %>% 
                  ifelse(identical(.type.es, "RR"), exp(.), .), 
                r.digits), "]"),
        nnt = ifelse(identical(.type.es, "RR"), 
                     NA, metapsyNNT(as.numeric(b[,1]), m.nnt.cer)) %>% 
                     ifelse(isTRUE(.raw.bin.es), nnt.g, .) %>% 
                     round(round.digits) %>% abs(),
        excluded = paste0(
          "Step-function selection model with cutpoints p=",
          paste(selmodel.steps*2, collapse = ", "), ". ",
          "The selection model parameter test was ",
          ifelse(LRTp < 0.05, "significant: ", "not significant: "),
          "\u03C7\u00B2","=", round(LRT, 3), 
          " (p=", scales::pvalue(LRTp), "). ",
          "The model was fitted using maximum likelihood estimation."
        )
      )
    })
  }
  
  rownames(mSelmodelRes) = "Selection model"
  
  if (!isTRUE(.selmodel.error)){
    message("- ", crayon::green("[OK] "),
            "Estimated effect using step function selection model (p=",
            paste(selmodel.steps*2, collapse = ", "), ").")
  } 
  
  
  
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

