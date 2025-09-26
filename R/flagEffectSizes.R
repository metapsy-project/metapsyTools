#' Flag extreme and/or implausible effect sizes
#'
#' @description This function flags effect sizes that are too extreme and/or implausible to merit inclusion in
#' meta-analytic models. Effect sizes are flagged based the size of the effect, achieved power, and study quality.
#' \loadmathjax
#'
#' @usage flagEffectSizes(data, 
#'                 reference = c("all", "dep", "psy", "ptsd"),
#'                 es.var = ".g", 
#'                 se.var = ".g_se",
#'                 high.rob.filter = "rob <= 2", 
#'                 power = "reference",
#'                 alpha = 0.05)
#'
#' @param data A `data.frame`, including effect size data following the 
#' [Metapsy data standard](https://docs.metapsy.org/data-preparation/format/), as 
#' created by [calculateEffectSizes()].
#' @param reference A `character` value, indicating the specific disorder for which reference values should be used.
#' Defaults to `"all"`, meaning that the reference values across all included disorders are used. Other supported 
#' values are `"dep"` (depression), `"psy"` (psychosis), and `"ptsd"` (PTSD). To improve calibration, 
#' the `"psy"` reference values may also be used for fields in which low effect sizes are more common,
#' or `"ptsd"` in high-effect size fields.
#' @param es.var A `character`, specifying the name of the variable containing the (pre-calculated) 
#' effect sizes in `data`. Defaults to `".g"`.
#' @param se.var A `character`, specifying the name of the variable containing the 
#' standard errors (square root of the variance) in `data`. Defaults to `".g_se"`.
#' @param high.rob.filter A `character`, including a filtering statement by which to 
#' include studies with high risk of bias/low study quality. Please note that the name of 
#' the variable must be included as a column in `data`.
#' @param power Reference value for defining a "strongly" underpowered study (\mjeqn{1-\beta}{1-\beta}). 
#' Defaults to `"reference"`, which corresponds to the first quartile of a 
#' reference power distribution derived from a large-scale meta-analytic 
#' database of psychological treatments (see Details). Also accepts user-defined 
#' power threshold in the (0,1) range (e.g., `power = 0.8` for 80% power).  
#' @param alpha The \mjeqn{\alpha}{\alpha} level used to determine "incompatibility" of an effect size (0.05 by default). See Details.
#'
#' @return 
#' Returns an object of class `flagEffectSizes`, with the following elements:
#' 
#' \itemize{
#'  \item `smd`. The provided effect sizes.
#'  \item `se`. The provided standard errors of the effect size.
#'  \item `flag.effect`. A logical, indicating if the effect was flagged for the first criterion (incompatible effect size).
#'  \item `flag.power`. A logical, indicating if the effect was flagged for the second criterion (strongly underpowered).
#'  \item `flag.rob`. A logical, indicating if the effect was flagged for the third criterion (risk of bias/low quality).
#'  \item `flags`. The numeric, indicating the number of flags for each effect size (0-3).
#'  \item `lookup`. A data frame, showing the empirical reference values that were used for the computations. Includes the 
#'  name of the reference effect size distribution (`distname`) and its parameters, the reference effect size based on
#'  which power is determined (`mu.fe`), as well as the first quartile of power in the reference sample (`q1.pwr`)
#'  \item `data`. The original data.frame provided in `data`.
#' }
#' 
#' Returned objects can be directly analyzed using [runMetaAnalysis()], which will then
#' remove flagged studies (3 flags) during the computation. 
#' 
#' @details
#' This flagging tool is based on a meta-epidemiological study of >2,000 effect 
#' sizes from the Metapsy databases, covering psychological treatment trials across 
#' 12 mental health problems (depression, generalized anxiety, social anxiety, 
#' specific phobia, panic disorder, obsessive-compulsive disorder, PTSD, complicated grief,
#' psychosis, suicidality). Based on this database, three indicators are 
#' used to flag trials:
#'
#' \itemize{
#'   \item \strong{Extreme effect sizes}: A set of candidate parametric distributions 
#'   (Pareto, exponential, Cauchy, Gamma, logistic, log-normal, and Weibull) was fitted
#'   to effect sizes from low risk-of-bias trials with inactive control groups in the Metapsy databases. 
#'   The best-fitting distribution was then selected. Effect sizes are flagged if they exceeded a critical threshold 
#'   corresponding to a Shannon information value (S-value) of 4.32, 
#'   roughly equating to the improbability of four heads in a row, or 
#'   \eqn{p} = 0.05.
#'
#'   \item \strong{Strong lack of power}: Each study's achieved statistical 
#'   power is computed relative to a bias-corrected reference effect size 
#'   estimated from low risk-of-bias trials (see `mu.fe` in the `lookup` element returned by the function). Instead of 
#'   using the conventional 80% threshold, trials are flagged by default if their power 
#'   falls in the lowest empirical quartile of this empirical distribution (`power = "reference"`). 
#'   Thus, by default, only the most severely underpowered trials are flagged.
#'
#'   \item \strong{Risk of bias / low study quality}: Trials not judged as "low risk" 
#'   (e.g. according to Cochrane's Risk of Bias tool) are flagged.
#' }
#'
#' For each effect size, the three indicators are summed (0–3 flags). Effects showing 
#' all three risk indicators are considered highly implausible and may be 
#' excluded by default from meta-analytic models. This is the default behavior if results
#' of `flagEffectSizes` are forwarded directly to [runMetaAnalysis()].
#' Reference values can be applied across all mental health problems (`"all"`), or calibrated to specific 
#' fields (`"dep"` = depression, `"psy"` = psychosis, `"ptsd"` = PTSD), depending on 
#' the user's selection for the `reference` argument. To improve calibration, 
#' the `"psy"` reference values may also be used for fields in which low effect sizes are more common,
#' or `"ptsd"` in high-effect size fields.
#' 
#' The development and performance of the flagging tool is described in greated detail in Harrer et al. (2025). 
#' Please note that reference values may change over time as the tool is updated using new 
#' versions of the living databases. Exact reference values and distributions 
#' used for the calculations will always be provided in the `lookup` data frame returned by the function. 
#' 
#' @examples
#' \dontrun{
#' data("depressionPsyCtr")
#' 
#' # Use as part of a meta-analysis pipeline
#' # Flagged studies are excluded
#' data %>% 
#'   checkDataFormat() %>% 
#'   calculateEffectSizes() %>% 
#'   filterPoolingData(condition_arm1 == "cbt") %>% 
#'   flagEffectSizes() %>% 
#'   runMetaAnalysis()
#' 
#' # Use depression-specific reference values
#' fe = flagEffectSizes(data, reference = "dep")
#' 
#' # Show results and used reference values
#' fe
#' fe$lookup
#' 
#' # Generate a plot of the empirical effects & 
#' # reference distribution
#' plot(fe, breaks = 20)
#' 
#' # Append number of flags to original data
#' # Then use in meta-regression and subgroup analysis
#' data$flags = fe$flags
#' M = runMetaAnalysis(data)
#' metaRegression(M$model.overall, ~ flags)
#' subgroupAnalysis(M, flags)
#' }
#' 
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @details For more details on the metapsyTools package, see the [Get Started](https://tools.metapsy.org/articles/metapsytools) vignette.
#'
#' @export flagEffectSizes
flagEffectSizes = function(data, reference = c("all", "dep", "psy", "ptsd"),
                           es.var = ".g", se.var = ".g_se",
                           high.rob.filter = "rob <= 2", power = "reference",
                           alpha = 0.05) {
  
  if (inherits(data, "metapsyDatabase")) { 
    x = data[["data"]] } else { x = data }
  
  if (!reference[1] %in% c("all", "dep", "psy", "ptsd")) {
    stop("'reference' muste be one of 'all', 'dep', 'psy', or 'ptsd'.")
  }
  
  robVar = strsplit(high.rob.filter, " ")[[1]][1]
  if (!robVar %in% colnames(x)){
    stop("'", robVar, "' variable not found in 'data'.", call. = FALSE)
  }
  
  if (identical(power[1], "reference")) {
    power = "q1"
  }
  
  robFilter = paste0("x$", high.rob.filter, " & !is.na(x$", robVar, ")")
  robMask = eval(parse(text = robFilter))
  
  if (sum(robMask)==0) {
    warning("No studies with risk of bias detected.",
            " Check if 'high.rob.filter' has been correctly applied.")
  }
  
  flagStudy(smd = abs(x[[es.var]]), se = x[[se.var]], rob = robMask, 
            power = power, pval = alpha, 
            reference = reference[1]) -> res
  res[["data"]] = x
  class(res) = c("flagEffectSizes", "list")
  return(res)
}


#' Print method for objects of class 'flagEffectSizes'
#'
#' Print S3 method for objects of class \code{flagEffectSizes}, as created by [flagEffectSizes()].
#'
#' @param x An object of class \code{flagEffectSizes}.
#' @param ... Additional arguments.
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers 
#' \email{p.cuijpers@@vu.nl}
#'
#' @export
#' @method print flagEffectSizes
print.flagEffectSizes = function(x, ...) {
  print.flagStudy(x, ...)
}

#' Plot method for objects of class 'flagEffectSizes'
#'
#' Plot S3 method for objects of class \code{flagEffectSizes}, as created by [flagEffectSizes()].
#'
#' @param x An object of class \code{flagEffectSizes}.
#' @param ... Additional arguments.
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @export
#' @method plot flagEffectSizes
plot.flagEffectSizes = function(x, ...) {
  plot.flagStudy(x, ...)
}








