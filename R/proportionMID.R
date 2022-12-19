#' Calculate the proportion of true effect sizes above a meaningful threshold
#'
#' Based on results of the [runMetaAnalysis()], this function allows to
#' estimate the proportion of true effect sizes that exceed a user-defined
#' meaningful (e.g. clinically relevant) threshold.
#'
#' @usage proportionMID(model, 
#'               mid = NULL, 
#'               which = "all", 
#'               test = "smaller", 
#'               plot = FALSE)
#'
#' @param model A class \code{runMetaAnalysis} object, created by the [runMetaAnalysis()] function.
#' @param mid A `numeric` value, indicating a clinically relevant effect threshold (e.g. a minimally
#' important difference; \mjeqn{MID}{MID}; Cuijpers et al., [2014](https://onlinelibrary.wiley.com/doi/full/10.1002/da.22249)) 
#' that should be used to estimate the
#' proportion of true effect sizes that exceed this cut-off. If the outcome measure used in `model` is 
#' Hedges' \mjeqn{g}{g}, the provided value should also be a standardized mean difference.
#' If the outcome measure of `model` is a risk ratio, the treshold should also be provided
#' as an (untransformed) risk ratio.
#' @param which The model in \code{model} that should be used to estimate the proportions.
#' Defaults to `"all"`, which means that proportions are calculated for all models.
#' Alternatively, possible values are \code{"overall"}, \code{"combined"}, \code{"lowest"}, 
#' \code{"highest"}, \code{"outliers"}, \code{"influence"} and \code{"rob"}, if these
#' models are available in the `model` object. If 
#' [correctPublicationBias()] has been run, `"trimfill"`, `limitmeta` and `selection`
#' are also possible options. It is also possible to concatenate model names, meaning
#' that proportions are calculated for all the supplied models.  
#' @param test By default, the function estimates the proportion of true effects _below_
#' the provided threshold in `mid` (`test="smaller"`). Alternatively, one can specify
#' `test="bigger"`. This will calculate the proportion of true effects _above_ the treshold.
#' @param plot Should a density plot illustrating the proportions be returned? Defaults
#' to `FALSE`. Please note that an S3 `plot` method is available for outputs of this
#' function even when `plot=FALSE` (see "Details").
#'
#' @return Returns an object of class \code{"proportionMID"}. An S3 `plot` method
#' is defined for this object class, which allows to create a density plot illustrating
#' the estimated proportions, using the model-based estimate of the pooled effect size
#' and between-study heterogeneity \mjeqn{\tau^2}{\tau^2}. 
#'
#'
#' @examples
#' \dontrun{
#' # Run meta-analysis; then estimate the proportion
#  # of true effect sizes that exceed a MID of -0.24
#' depressionPsyCtr %>% 
#'   filterPoolingData(condition_arm1 == "cbt") %>% 
#'   runMetaAnalysis()  -> x
#' 
#' proportionMID(x, mid = -0.24)
#' proportionMID(x, mid = -0.24, "outliers") %>% plot()
#' 
#' 
#' # If bootstrap CIs are requested in runMetaAnalysis,
#' # calculation of CIs around p is possible for all available
#' # models.
#' depressionPsyCtr %>% 
#'   filterPoolingData(condition_arm1 == "cbt") %>% 
#'   runMetaAnalysis(i2.ci.boot = TRUE, nsim.boot = 1000) %>% 
#'   correctPublicationBias() -> x
#' 
#' proportionMID(x, mid = -0.33)
#' 
#' # Run meta-analysis based on RRs; then estimate proportion
#' # of true effects bigger than the defined threshold
#' depressionPsyCtr %>% 
#'   filterPoolingData(
#'     condition_arm1 %in% c("cbt", "pst", 
#'                           "dyn", "3rd wave")) %>% 
#'   runMetaAnalysis(es.measure = "RR") %>% 
#'   proportionMID(mid = 1.13, test = "bigger")
#' }
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{runMetaAnalysis}}, \code{\link{correctPublicationBias}}
#'
#' @details The `proportionMID` function implements an approach to estimate
#' the proportion of true effect sizes exceeding a (scientifically or clinically)
#' relevant threshold, as proposed by Mathur & VanderWeele ([2019](https://www.tandfonline.com/doi/10.1080/01621459.2018.1529598)). These estimated 
#' proportions have been suggested as a useful metric to determine
#' the impact that between-study heterogeneity in a meta-analysis has on the 
#' "real-life" interpretation of results.
#' 
#' If, for example, a pooled effect is significant, high between-study heterogeneity
#' can still mean that a substantial proportion of true effects in the studies
#' population are practicially irrelevant, or even negative. Conversely
#' overall non-significant effects, in the face of large heterogeneity, can still
#' mean that a substantial proportion of studies have non-negligible _true_ effects.
#' 
#' As recommended by Mathur & VanderWeele ([2019](https://www.tandfonline.com/doi/10.1080/01621459.2018.1529598)), 
#' the `proportionMID` function
#' also automatically calculates the proportion of true effects exceeding the 
#' _"inverse"_ of the user-defined effect (e.g., if `mid=-0.24`, by default,
#' the function also estimates the proportion of true effects that is larger
#' than \mjeqn{g}{g}=0.24; note the changed sign). This can be used to check for e.g. clinically
#' relevant negative effects.
#' 
#' When the `plot` method is used, or when `plot` is set to `TRUE` in the function,
#' a plot showing the assumed distribution of true effects based on the estimated
#' meta-analytic model is created. Notably, is is assumed that the random-effects
#' distribution of true effect sizes is approximately normal. This simplifying assumption
#' is required for this (and many other meta-analytic methods; 
#' Jackson & White, [2018](https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201800071)) 
#' to hold.
#' 
#' Confidence intervals provided by the functions are calculated using the
#' asymptotic closed-form solution derived using the Delta method in Mathur &
#' VanderWeele ([2019](https://www.tandfonline.com/doi/10.1080/01621459.2018.1529598)). 
#' Following their recommendations, a warning is printed
#' when \mjeqn{p}{p}<0.15 or \mjeqn{p}{p}>0.85, since in this case the asymptotic
#' CIs should be interpreted cautiously; CIs based on boostrapping would be preferable
#' in this scenario and can be calculated using the [MetaUtility::prop_stronger()]
#' function.
#' 
#' @references 
#' Cuijpers, P., Turner, E. H., Koole, S. L., Van Dijke, A., & Smit, F. (2014). 
#' What is the threshold for a clinically relevant effect? The case of 
#' major depressive disorders. _Depression and Anxiety, 31_(5), 374-378.
#' 
#' Jackson, D., & White, I. R. (2018). When should meta-analysis avoid making 
#' hidden normality assumptions?. _Biometrical Journal, 60_(6), 1040-1058.
#' 
#' Mathur, M. B., & VanderWeele, T. J. (2019). Sensitivity analysis for unmeasured 
#' confounding in meta-analyses. _Journal of the American Statistical Association_.
#' 
#' @import dplyr numDeriv
#' @importFrom stats dffits model.matrix rnorm rstudent quantile sd qnorm dnorm
#' @importFrom utils combn argsAnywhere
#' @export proportionMID

proportionMID = function(model, mid = NULL, which = "all", 
                         test = "smaller", plot = FALSE){
  
  x = model
  
  if (is.null(mid)){
    stop("A value of 'mid' must be specified.")
  }
  
  if (!(test=="smaller" | test=="bigger")){
    stop("'test' must be either 'smaller' (default) or 'bigger'.")
  }
  
  if (!(which[1] == "all" | 
      any(which %in% c("overall", "combined",
                   "lowest.highest", "outliers",
                   "influence", "rob", "threelevel",
                   "threelevel.che", "trimfill", 
                   "limitmeta", "selection")))){
    stop("'which' must be either 'all' or one of",
         " 'overall', 'combined', 'lowest.highest',",
         " 'outliers', 'influence', 'rob',",
         " 'threelevel', 'threelevel.che', 'trimfill',",
         " 'limitmeta', 'selection'.")}
  
  # Extract relevant information
  if (identical(which[1],"all")){
    run = paste0("model.", x$which.run)
    if (!is.null(x$correctPublicationBias)){
      x$model.trimfill = x$correctPublicationBias$model.trimfill
      x$model.limitmeta = x$correctPublicationBias$model.limitmeta
      x$model.selection = x$correctPublicationBias$model.selection
      run = c(run, c("model.trimfill", "model.limitmeta", "model.selection"))
    }
  } else {
    run = paste0("model.", which)
    if (!is.null(x$correctPublicationBias)){
      x$model.trimfill = x$correctPublicationBias$model.trimfill
      x$model.limitmeta = x$correctPublicationBias$model.limitmeta
      x$model.selection = x$correctPublicationBias$model.selection
    }
  }
  if ("model.lowest.highest" %in% run){
    run = c(run[-which(run %in% "model.lowest.highest")],
            "model.lowest", "model.highest")
  }
  
  lapply(x[run], function(z){
    if (identical(class(z)[1], "metagen") | 
        identical(class(z)[1], "metabin")){
      with(z, {
        c(es = ifelse(method.tau=="FE", TE.common, TE.random),
          se = ifelse(method.tau=="FE", seTE.common, seTE.random),
          tau2 = tau2, se.tau2 = se.tau2)
      }) -> ret
    }
    if (identical(class(z)[1], "rma.mv")){
      if (!is.numeric(z$se.sigma2)){
        z$se.sigma2 = NA
      }
      with(z, {
        c(es = b[[1,1]], se = se[1],
          tau2 = sigma2[1], 
          se.tau2 = ifelse(is.na(se.sigma2[1]), NA, se.sigma2[1,1]))
      }) -> ret
    }
    if (identical(class(z)[1], "limitmeta")){
      if (!is.numeric(z$se.tau2)){
        z$se.tau2 = NA
      }
      with(z, {
        c(es = TE.adjust, se = seTE.adjust,
          tau2 = tau^2, 
          se.tau2 = ifelse(is.na(se.tau2), NA, se.tau2))
      }) -> ret
    }
    if (identical(class(z)[1], "rma.uni.selmodel")){
      if (!is.numeric(z$se.tau2)){
        z$se.tau2 = NA
      }
      with(z, {
        c(es = b[[1,1]], se = se[1],
          tau2 = tau2, 
          se.tau2 = se.tau2)
      }) -> ret
    }
    return(ret)
  }) %>% do.call(rbind, .) -> dat.calc
  
  # If es.measure is "RR", log-transform MID
  if (identical(x$.type.es[1], "RR")){
    mid = log(mid)
  }
  
  # Calculate q-hat & asymptomatic CIs
  if (identical(test[1], "smaller")){
    with(as.data.frame(dat.calc), {
      data.frame(
        "p.below" = pnorm((mid - es) / sqrt(tau2)),
        "p.above.rev" = 1 - pnorm((-mid - es) / sqrt(tau2)))
    }) -> sw
  }
  if (identical(test[1], "bigger")){
    with(as.data.frame(dat.calc), {
      data.frame(
        "p.below" = 1- pnorm((mid - es) / sqrt(tau2)),
        "p.above.rev" = pnorm((-mid - es) / sqrt(tau2)))
    }) -> sw
  }
  
  sw %>% 
    {cbind(dat.calc, .)} %>% 
    {SE = sqrt((.$se^2/.$tau2) + 
                 ((.$se.tau2^2*(mid-.$es)^2)/(4*.$tau2^3))) *
      dnorm(((mid-.$es)/sqrt(.$tau2)))
    SE.rev = sqrt((.$se^2/.$tau2) + 
                    ((.$se.tau2^2*(-mid-.$es)^2)/(4*.$tau2^3))) *
      dnorm(((-mid-.$es)/sqrt(.$tau2)))
    .$p.below.lo = .$p.below + qnorm(0.025)*SE
    .$p.below.hi = .$p.below - qnorm(0.025)*SE
    .$p.above.rev.lo = .$p.above.rev + qnorm(0.025)*SE.rev
    .$p.above.rev.hi = .$p.above.rev - qnorm(0.025)*SE.rev; 
    .} -> res.unrounded
  
  res.unrounded %>% 
    {within(.,{
      p.below.lo = ifelse(p.below.lo < 0, 0, p.below.lo) 
      p.below.hi = ifelse(p.below.hi > 1, 1, p.below.hi) 
      p.above.rev.lo = ifelse(p.above.rev.lo < 0, 0, p.above.rev.lo)
      p.above.rev.hi = ifelse(p.above.rev.hi > 1, 1, p.above.rev.hi)
      tau = round(sqrt(tau2), x$round.digits+1)
      p.below = sprintf(paste0("%.", x$round.digits,"f"), p.below*100)
      p.below.ci = paste0(
        "[", sprintf(paste0("%.", x$round.digits,"f"), p.below.lo*100),
        "; ", sprintf(paste0("%.", x$round.digits,"f"), p.below.hi*100), "]")
      p.above.rev = sprintf(paste0("%.", x$round.digits,"f"), p.above.rev*100)
      p.above.rev.ci = paste0(
        "[", sprintf(paste0("%.", x$round.digits,"f"), p.above.rev.lo*100),
        "; ", sprintf(paste0("%.", x$round.digits,"f"), p.above.rev.hi*100), "]")
    })} %>% 
    {.[,c("es", "tau", "p.below", "p.below.ci",
          "p.above.rev", "p.above.rev.ci")]
    } -> Qhat
  
  # Convert if es.measure was RR
  if (identical(x$.type.es[1], "RR")){
    Qhat$es = round(exp(Qhat$es), x$round.digits)
    mid = exp(mid)
    mid.neg = round(exp(-log(mid)),2)
    es.name = "rr"
  } else {
    Qhat$es = round(Qhat$es, x$round.digits)
    es.name = "g"
    mid.neg = -mid
  }
  
  # Define colnames
  Qhat %>% 
    {colnames(.) = 
      c(es.name, "tau",
        paste0("p(", ifelse(x$.type.es=="RR", "rr", "g"), 
               ifelse(test=="smaller", "<", ">"), mid, ")"),
        "95% ci", paste0("p(", ifelse(x$.type.es=="RR", "rr", "g"), 
               ifelse(test=="smaller", ">", "<"), mid.neg, ")"), "95% ci");
    .} -> Qhat
  
  
  # Create Model list
  model.list = c(
    "model.overall" = "Overall", 
    "model.combined" = "Combined",
    "model.lowest" = "One ES/study (lowest)",
    "model.highest" = "One ES/study (highest)",
    "model.outliers" = "Outliers removed",
    "model.influence" = "Influence Analysis",
    "model.rob" = ifelse(is.null(x$model.rob), 
                         "model.rob", x$model.rob$title2),
    "model.threelevel" = "Three-Level Model",
    "model.threelevel.che" = "Three-Level Model (CHE)",
    "model.trimfill" = "Trim-and-fill method",
    "model.limitmeta" = "Limit meta-analysis",
    "model.selection" = "Selection model")
  
  Qhat[names(model.list),] %>% 
    {.[rownames(.) %in% names(model.list),]} %>% 
    {rownames(.) = model.list[rownames(.)];.} -> Qhat
  
  # Remove NA fields
  Qhat[Qhat=="[NA; NA]"] = "-"
  
  # Generate plot data
  res.unrounded %>% 
    with({
      data.frame(es = es, tau = sqrt(tau2),
                 lo=es-(sqrt(tau2)*3), hi=es+(sqrt(tau2)*3),
                 mid = ifelse(x$.type.es=="g", mid, log(mid)), 
                 mid.neg = ifelse(x$.type.es=="g", mid.neg, log(mid.neg)))
    }) %>% 
    {rownames(.) = rownames(res.unrounded);.} -> plot.dat
  
  plot.dat[names(model.list),] %>% 
    {.[rownames(.) %in% names(model.list),]} %>% 
    {rownames(.) = model.list[rownames(.)];.} -> plot.dat
  
  # Generate bell curve
  if (identical(x$.type.es[1], "RR")){xlab = "log(RR)"
  } else { xlab = "Hedges' g" }
  
  if (plot){
    
    par(mfrow=c(1,1))
    if (nrow(plot.dat) > 1){
      par(mfrow=c(1,2))
    }
    if (nrow(plot.dat) > 2){
      par(mfrow=c(2,2))
    } 
    if (nrow(plot.dat) > 4){
      par(mfrow=c(3,2))
    }
    
    for (i in 1:nrow(plot.dat)){
      X = seq(min(plot.dat), max(plot.dat),length=2000)
      y = dnorm(X, mean=plot.dat$es[i], sd=plot.dat$tau[i])
      plot(X,y,type="l", ylab="", yaxt="n",
           xlab = xlab)
      
      if (identical(test[1], "smaller")){
        x.mid = seq(abs(plot.dat$mid[i]), max(plot.dat), length=2000)
        y.mid = dnorm(x.mid, mean=plot.dat$es[i], sd=plot.dat$tau[i])
        polygon(c(abs(plot.dat$mid[i]),x.mid,max(plot.dat)),
                c(0,y.mid,0), col = adjustcolor("lightblue", alpha.f=0.5))
        
        x.mid.neg = seq(min(plot.dat), -abs(plot.dat$mid.neg[i]), length=2000)
        y.mid.neg = dnorm(x.mid.neg, mean=plot.dat$es[i], sd=plot.dat$tau[i])
        polygon(c(min(plot.dat),x.mid.neg,-abs(plot.dat$mid.neg[i])),
                c(0,y.mid.neg,0), col = adjustcolor("lightyellow", alpha.f=0.5))
        abline(v = plot.dat$mid[i], col="gray", lty=2)
        abline(v = plot.dat$mid.neg[i], col="gray", lty=2)
        abline(v = plot.dat$es[i], col="black", lwd=2)
      } else {
        x.mid = seq(plot.dat$mid[i], max(plot.dat), length=2000)
        y.mid = dnorm(x.mid, mean=plot.dat$es[i], sd=plot.dat$tau[i])
        polygon(c(plot.dat$mid[i],x.mid,max(plot.dat)),
                c(0,y.mid,0), col = adjustcolor("lightblue", alpha.f=0.5))
        
        x.mid.neg = seq(min(plot.dat), plot.dat$mid.neg[i], length=2000)
        y.mid.neg = dnorm(x.mid.neg, mean=plot.dat$es[i], sd=plot.dat$tau[i])
        polygon(c(min(plot.dat),x.mid.neg,plot.dat$mid.neg[i]),
                c(0,y.mid.neg,0), col = adjustcolor("lightyellow", alpha.f=0.5))
        abline(v = plot.dat$mid[i], col="gray", lty=2)
        abline(v = plot.dat$mid.neg[i], col="gray", lty=2)
        abline(v = plot.dat$es[i], col="black", lwd=2)
      }
      title(rownames(plot.dat)[i]) 
    }
  }

  # Return 
  ret.obj = list(qhat = Qhat, xlab = xlab, 
                 plot.dat = plot.dat, 
                 test = test)
  class(ret.obj) = c("proportionMID")
  return(ret.obj)
}



#' Print method for objects of class 'proportionMID'
#'
#' Print S3 method for objects of class \code{proportionMID}.
#'
#' @param x An object of class \code{proportionMID}.
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
#' @importFrom stringr str_sub str_replace_all
#'
#' @export
#' @method print proportionMID

print.proportionMID = function(x, ...){
  res = x
  cat(crayon::blue$bold("Proportion of meaningful true effects "))
  cat(crayon::blue$bold(
    "------------------------------------------------ \n"))
  cat(crayon::green$bold(paste0("\u2192 Assuming an effect of ",
                                ifelse(res$xlab=="Hedges' g", "g", "rr"), "=", 
                                ifelse(res$xlab=="Hedges' g", 
                                       res$plot.dat$mid[1], 
                                       round(exp(res$plot.dat$mid[1]), 3)),
                                " as threshold \n")))
  dat = res$qhat
  colnames(dat)[ncol(dat)] = " 95% ci"
  tbl = dplyr::as_tibble(cbind(Model = rownames(dat), dat))
  old = options(pillar.bold=TRUE)
  tbl.format = format(tbl)[-c(1,3)]
  tbl.format[-1] = lapply(tbl.format[-1], 
                          function(x) stringr::str_sub(x, 19))
  tbl.format[1] = stringr::str_sub(tbl.format[1], 3)
  tbl.format = lapply(tbl.format, function(x){
    stringr::str_replace_all(x, "`", " ")})
  tbl.format = lapply(tbl.format, function(x){
    stringr::str_replace_all(x, "-", "\U2013")})
  cat(do.call(c, tbl.format), sep="\n")
  options(old)
  
  if (min(as.numeric(res$qhat[,3]))<15 | max(as.numeric(res$qhat[,3]))>85 |
      min(as.numeric(res$qhat[,5]))<15 | max(as.numeric(res$qhat[,5]))>85){
    warning("Some estimated proportions are p<15 or p>85.",
            " In this case, the closed-form asymptotic CIs returned by this function",
            " should be interpreted with caution.")
  }
}



#' Plot method for objects of class 'proportionMID'
#'
#' Plot S3 method for objects of class \code{proportionMID}.
#'
#' @param x An object of class \code{proportionMID}.
#' @param ... Additional arguments.
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @importFrom stats dnorm
#' @importFrom grDevices adjustcolor
#' @importFrom graphics polygon title
#'
#' @export
#' @method plot proportionMID

plot.proportionMID = function(x, ...){
  
  plot.dat = x$plot.dat
  xlab = x$xlab
  test = x$test
  
  par(mfrow=c(1,1))
  if (nrow(plot.dat) > 1){
    par(mfrow=c(1,2))
  }
  if (nrow(plot.dat) > 2){
    par(mfrow=c(2,2))
  } 
  if (nrow(plot.dat) > 4){
    par(mfrow=c(3,2))
  }
  
  for (i in 1:nrow(plot.dat)){
    X = seq(min(plot.dat), max(plot.dat),length=2000)
    y = dnorm(X, mean=plot.dat$es[i], sd=plot.dat$tau[i])
    plot(X,y,type="l", ylab="", yaxt="n",
         xlab = xlab)
    
    if (identical(test[1], "smaller")){
      x.mid = seq(abs(plot.dat$mid[i]), max(plot.dat), length=2000)
      y.mid = dnorm(x.mid, mean=plot.dat$es[i], sd=plot.dat$tau[i])
      polygon(c(abs(plot.dat$mid[i]),x.mid,max(plot.dat)),
              c(0,y.mid,0), col = adjustcolor("lightblue", alpha.f=0.5))
      
      x.mid.neg = seq(min(plot.dat), -abs(plot.dat$mid.neg[i]), length=2000)
      y.mid.neg = dnorm(x.mid.neg, mean=plot.dat$es[i], sd=plot.dat$tau[i])
      polygon(c(min(plot.dat),x.mid.neg,-abs(plot.dat$mid.neg[i])),
              c(0,y.mid.neg,0), col = adjustcolor("lightyellow", alpha.f=0.5))
      abline(v = plot.dat$mid[i], col="orange3", lty=3, lwd=1.5)
      abline(v = plot.dat$mid.neg[i], col="lightblue", lty=3, lwd=1.5)
      abline(v = plot.dat$es[i], col="black", lwd=2)
    } else {
      x.mid = seq(plot.dat$mid[i], max(plot.dat), length=2000)
      y.mid = dnorm(x.mid, mean=plot.dat$es[i], sd=plot.dat$tau[i])
      polygon(c(plot.dat$mid[i],x.mid,max(plot.dat)),
              c(0,y.mid,0), col = adjustcolor("lightblue", alpha.f=0.5))
      
      x.mid.neg = seq(min(plot.dat), plot.dat$mid.neg[i], length=2000)
      y.mid.neg = dnorm(x.mid.neg, mean=plot.dat$es[i], sd=plot.dat$tau[i])
      polygon(c(min(plot.dat),x.mid.neg,plot.dat$mid.neg[i]),
              c(0,y.mid.neg,0), col = adjustcolor("lightyellow", alpha.f=0.5))
      abline(v = plot.dat$mid[i], col="lightblue", lty=3, lwd=1.5)
      abline(v = plot.dat$mid.neg[i], col="orange3", lty=3, lwd=1.5)
      abline(v = plot.dat$es[i], col="black", lwd=2)
    }
    title(rownames(plot.dat)[i]) 
  }
}







