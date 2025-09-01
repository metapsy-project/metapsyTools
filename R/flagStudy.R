#' Flagging engine for internal use
#' @keywords internal
#' @importFrom stats qnorm qweibull qgamma qlnorm qlogis qcauchy qexp 
#' @importFrom crayon green yellow cyan bold
flagStudy = function(smd, se, rob, power = "q1", pval = 0.05,
                     reference = c("all", "dep", "psy", "ptsd"), ...){
  
  # Argument checking
  if (!reference[1] %in% c("all", "dep", "psy", "ptsd"))
    stop("\"reference\" must be either \"all\", \"dep\", \"psy\" or \"ptsd\".")
  if (!is.numeric(pval[1]) || identical(pval[1],0) || pval[1]>=1)
    stop("\"pval\" must be a numeric between 0 and 1.")
  if (!is.logical(rob))
    stop("\"rob\" must be a logical vector.")
  if (!all(length(smd) == length(se), length(se) == length(rob)))
    stop("lengths of \"smd\", \"se\", and \"rob\" do not match.")
  if (!(identical(power[1], "q1") ||
        all(is.numeric(power[1]), power[1]<1, power[1]>0)))
    stop("\"power\" must by a numeric [0,1], or set to \"q1\".")
  
  # Lookup table of pre-computed thresholds by dataset type
  lookup = structure(list(type = c("all", "dep", "gad", "psy", "ptsd"), 
                          shape = c(1.22780106290059, 1.49556995998092, 2.08900202717736, 
                                    1.2491159032556, 2.69699379305318), rate = c(1.82800375781207, 
                                                                                 2.343560247288, 2.78725113029231, 2.78807085068211, 2.02655752359899
                                    ), mu.fe = c(0.455964509009351, 0.484125449169705, 0.678273877252083, 
                                                 0.30882193290888, 1.06951362680109), q1.pwr = c(0.3240959446545, 
                                                                                                 0.410179676205664, 0.57613016526012, 0.218204329492164, 0.748225230452593
                                                 ), distname = c("gamma", "gamma", "gamma", "gamma", "gamma"
                                                 )), row.names = c(NA, -5L), class = "data.frame")
  rownames(lookup) = lookup$type
  lookup$type = NULL
  
  message("- ", crayon::green("[OK] "), "Using reference values for database \"", 
          reference[1], "\".", appendLF = TRUE)
  
  # Threshold for extreme SMD (based on Weibull 1-pval quantile)
  list(x = list(
    distname = lookup[reference[1],"distname"],
    estimate = lookup[reference[1],!colnames(lookup) %in%
                        c("mu.fe", "q1.pwr", "distname")]),
    p = (1-pval)) -> args
  do.call(get.critical, args) -> smd.crit
  message("- ", crayon::green("[OK] "), "Flagging all SMDs greater than ", 
          round(smd.crit,2), " as extreme.", appendLF = TRUE)
  extreme = as.numeric(abs(smd) > smd.crit)
  
  # Power thresholding
  if (identical(power, "q1"))
    power = lookup[reference[1],"q1.pwr"]
  message("- ", crayon::green("[OK] "), "Assuming minimum power (1-beta) of ", 
          round(power*100,1), "%.", appendLF = TRUE)
  lookup[reference[1],"mu.fe"]/(qnorm(1 - pval / 2) + qnorm(power)) -> se.crit
  underpowered = as.numeric(!(se.crit >= se))
  biased = as.numeric(rob) 
  flags = extreme + underpowered + biased
  flagged = flags==3
  
  # Output with custom class for printing
  ret = list(smd = smd, se = se,
             flag.effect = extreme==1,
             flag.power = underpowered==1,
             flag.rob = biased==1,
             flags = flags, 
             lookup = lookup[reference[1],])
  class(ret) = c("flagStudy", "list")
  return(ret)
}

#' Custom print method for flagStudy objects
#' @keywords internal
#' @method print flagStudy
print.flagStudy = function(x, ...) {
  data.frame(smd = x$smd, se = x$se, flag.effect = x$flag.effect, 
             flag.power = x$flag.power, flag.rob = x$flag.rob, flags = x$flags) -> x
  cat("---------------------------------------------------------\n")
  cat("N effect sizes with 1+ flag:  ", sum(x$flags>=1, na.rm=TRUE),
      " (", round((sum(x$flags>=1, na.rm=TRUE)/nrow(x))*100,1), "%) \n",sep = "")
  cat("N effect sizes with 2+ flags: ", sum(x$flags>=2, na.rm=TRUE),
      " (", round((sum(x$flags>=2, na.rm=TRUE)/nrow(x))*100,1), "%) \n",sep = "")
  cat("N effect sizes with 3 flags:  ", sum(x$flags>=3, na.rm=TRUE),
      " (", round((sum(x$flags>=3, na.rm=TRUE)/nrow(x))*100,1), "%) \n",sep = "")
  cat("---------------------------------------------------------\n")
  cat("---------------------------------------------------------\n")
  cat("N effect sizes with flag 1 (extreme effect): ", sum(x$flag.effect, na.rm=TRUE),
      " (", round((sum(x$flag.effect, na.rm=TRUE)/nrow(x))*100,1), "%) \n",sep = "")
  cat("N effect sizes with flag 2 (underpowered):   ", sum(x$flag.power, na.rm=TRUE),
      " (", round((sum(x$flag.power, na.rm=TRUE)/nrow(x))*100,1), "%) \n",sep = "")
  cat("N effect sizes with flag 3 (risk of bias):   ", sum(x$flag.rob, na.rm=TRUE),
      " (", round((sum(x$flag.rob, na.rm=TRUE)/nrow(x))*100,1), "%) \n",sep = "")
  cat("---------------------------------------------------------\n\n")
  within(x, {
    flag.effect = ifelse(flag.effect, "yes", "no")
    flag.power = ifelse(flag.power, "yes", "no")
    flag.rob = ifelse(flag.rob, "yes", "no")
  }) -> x
  print.data.frame(head(x, ...))
  invisible(x)
}

#' Custom plot method for flagStudy objects
#' @keywords internal
#' @method plot flagStudy
plot.flagStudy = function(x, ...) {
  message("Generating reference distribution plot...")
  plotDensityHist(x, ...)
}
