#' Calculate effect sizes
#'
#' This is a function to calculate effect sizes of meta-analysis data 
#' prepared in the [Metapsy data format](https://docs.metapsy.org/data-preparation/format/).
#'
#' @usage calculateEffectSizes(data,
#'                      funcs.g = list(g.m.sd = g.m.sd,
#'                                     g.change.m.sd = g.change.m.sd,
#'                                     g.binary = g.binary,
#'                                     g.precalc = g.precalc),
#'                      funcs.rr = list(rr.binary = rr.binary,
#'                                      rr.precalc = rr.precalc),
#'                      include.switched.arms = FALSE,
#'                      change.sign = NULL,
#'                      impute.response = FALSE,
#'                      vars.for.id = c("study", "outcome_type",
#'                                      "instrument", "time",
#'                                      "time_weeks",
#'                                      "rating"),
#'                      .condition = "condition",
#'                      .condition.specification = "multi",
#'                      .groups.column.indicator = c("_arm1", "_arm2"),
#'                      .trt.indicator = "arm",
#'                      .impute.response.vars = c(m.trt.pre = "baseline_m_arm1", 
#'                                                m.trt.post = "mean_arm1", 
#'                                                sd.trt.post = "sd_arm1", 
#'                                                n.trt = "n_arm1",
#'                                                m.ctr.pre = "baseline_m_arm2", 
#'                                                m.ctr.post = "mean_arm2", 
#'                                                sd.ctr.post = "sd_arm2", 
#'                                                n.ctr = "n_arm2"))
#'
#'
#' @param data Meta-analysis data set formatted using the 
#' [Metapsy guidelines](https://tools.metapsy.org/articles/metapsyTools.html#required-data-structure).
#' @param funcs.g \code{list} of functions. These functions will be used to calculate the 
#' effect sizes (Hedges' _g_) based on the raw data (see Details).
#' @param funcs.rr \code{list} of functions. These functions will be used to calculate risk ratios 
#' based on the raw event data (see Details).
#' @param include.switched.arms \code{logical}. Should all unique arm \emph{comparisons} (in lieu of unique arm \emph{combinations}) be
#' calculated? Default is \code{FALSE}. 
#' @param change.sign \code{character}. Name of a \code{logical} column in \code{data}, encoding if the
#' sign of a calculated effect size should be reversed (\code{TRUE}) or not (\code{FALSE}). Set to \code{NULL} (default) if
#' no changes should be made.
#' @param impute.response `logical`. When calculating the (log)-risk ratios, should response rates be computed
#' using the [imputeResponse()] function? `FALSE` by default. If defined, the column specified in `change.sign`
#' will also be considered when calculating the response.
#' @param vars.for.id \code{character} vector, containing column names of all variables 
#' used to construct unique comparison IDs.
#' @param .condition \code{character}. The prefix of the two variables in \code{data} in 
#' which the conditions (e.g. "guided iCBT", "waitlist") of the trial arm comparison are stored.
#' @param .condition.specification \code{character}. The prefix of the two variables in the dataset
#' which provide a "specification" of the trial arm condition in multiarm trials.
#' @param .groups.column.indicator \code{character}. A character vector with two elements, 
#' representing the suffix used to differentiate between the first and second arm in a comparison.
#' @param .trt.indicator \code{character}. A character specifying the name used to indicate the treatment arm.
#' @param .impute.response.vars `list`. Named list with the names of columns in `data` and the specific argument in
#' [imputeResponse()] they should be used for.
#' 
#' @return \code{calculateEffectSizes} returns the meta-analysis data set as 
#' class \code{data.frame} in wide format (if results are saved to a variable). 
#' It also generates the following columns, wich are added to the data:
#' 
#' - `.id`: Unique identifier for a trial arm comparison/row.
#' - `.g`: Calculated effect size (Hedges' _g_).
#' - `.g_se`: Standard error of Hedges' _g_.
#' - `.log_rr`: Calculated effect size (_logRR_).
#' - `.log_rr_se`: Standard error of _logRR_.
#' - `.event_arm1`: Number of events (responders, remission, deterioration cases) in the first trial arm.
#' - `.event_arm2`: Number of events (responders, remission, deterioration cases) in the second trial arm.
#' - `.totaln_arm1`: Total sample size in the first trial arm.
#' - `.totaln_arm2`: Total sample size in the second trial arm.
#' 
#' @examples
#' \dontrun{
#' data("depressionPsyCtr")
#' depressionPsyCtr %>%
#'     checkDataFormat() %>%
#'     checkConflicts() %>% 
#'     calculateEffectSizes()
#' }
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com}, 
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, 
#' Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{checkDataFormat}}, \code{\link{imputeResponse}}
#' 
#' @details By default, `calculateEffectSizes` calculates the 
#' small-sample bias corrected standardized mean difference  (Hedges' _g_)
#' and log-risk ratios, as well their respective standard errors, if 
#' adequate raw effect size data is available for a comparison.
#' 
#' It is essential that the data set in `data` contains a few required columns
#' for this to work. An overview of the required data format is provided
#' on the ["Get Started"](https://tools.metapsy.org/articles/metapsytools) 
#' page of the `metapsyTools` documentation. 
#' 
#' **Standardized mean differences (Hedges' _g_)** can be calculated from the following
#' column types:
#' 
#' - (1) Continuous Outcome Data
#'   - **`mean_arm1`**: Mean of the outcome in the first arm at the measured time point.
#'   - **`mean_arm2`**: Mean of the outcome in the second arm at the measured time point.
#'   - **`sd_arm1`**: Standard deviation of the outcome in the first arm at the measured time point.
#'   - **`sd_arm2`**: Standard deviation of the outcome in the second arm at the measured time point.
#'   - **`n_arm1`**: Sample size in the first trial arm.
#'   - **`n_arm2`**: Sample size in the second trial arm.
#' - (2) Change Score Data
#'   - **`mean_change_arm1`**: Mean score change between baseline and the measured time point in the first arm.
#'   - **`mean_change_arm2`**: Mean score change between baseline and the measured time point in the second arm.
#'   - **`sd_change_arm1`**: Standard deviation of the mean change in the first arm.
#'   - **`sd_change_arm2`**: Standard deviation of the mean change in the second arm.
#'   - **`n_change_arm1`**: Sample size in the first trial arm.
#'   - **`n_change_arm2`**: Sample size in the second trial arm.
#' - (3) Dichotomous Outcome Data
#'   - **`event_arm1`**: Number of events (responders, remission, deterioration cases) in the first trial arm.
#'   - **`event_arm2`**: Number of events (responders, remission, deterioration cases) in the second trial arm.
#'   - **`totaln_arm1`**: Sample size in the first trial arm.
#'   - **`totaln_arm2`**: Sample size in the second trial arm.
#' - (4) Pre-calculated Hedges' _g_
#'   - **`precalc_g`**: The pre-calculated value of Hedges' _g_ (small-sample bias corrected standardized mean difference; [Hedges, 1981](https://journals.sagepub.com/doi/10.3102/10769986006002107)).
#'   - **`precalc_g_se`**: Standard error of _g_.
#' 
#' The **log-risk ratio** and its standard error can be calculated from the followin column types:
#' 
#' - (1) Dichotomous Outcome Data
#'   - **`event_arm1`**: Number of events (responders, remission, deterioration cases) in the first trial arm.
#'   - **`event_arm2`**: Number of events (responders, remission, deterioration cases) in the second trial arm.
#'   - **`totaln_arm1`**: Sample size in the first trial arm.
#'   - **`totaln_arm2`**: Sample size in the second trial arm.
#' - (2) Pre-calculated log-risk ratio
#'   - **`precalc_log_rr`**: The pre-calculated value of the log-risk ratio logRR, comparing events in the first arm to events in the second arm.
#'   - **`precalc_log_rr_se`**: The standard error of the log-risk ratio logRR, comparing events in the first arm to events in the second arm.
#' 
#' 
#' Other functions can be added to the list provided to `funcs.g` and `funcs.rr`. 
#' However, results of the function must result in a \code{data.frame} that contains
#' the following columns:
#' 
#' - `.id`: Unique identifier for a trial arm comparison/row.
#' - `.g`: Calculated effect size (Hedges' _g_).
#' - `.g_se`: Standard error of Hedges' _g_.
#' - `.log_rr`: Calculated effect size (logRR).
#' - `.log_rr_se`: Standard error of logRR.
#' - `.event_arm1`: Number of events (responders, remission, deterioration cases) in the first trial arm.
#' - `.event_arm2`: Number of events (responders, remission, deterioration cases) in the second trial arm.
#' - `.totaln_arm1`: Total sample size in the first trial arm.
#' - `.totaln_arm2`: Total sample size in the second trial arm.
#' 
#' It is possible to set one or several of these column entries to `NA`; but the columns
#' themselves must be included.
#' 
#' For more details see the [Get Started](https://tools.metapsy.org/articles/metapsytools) vignette.
#'
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom dplyr all_of select filter mutate arrange
#' @importFrom purrr pmap_dfr map map2
#' @importFrom esc esc_mean_sd
#' @importFrom stringr str_remove str_replace str_replace_all
#' @import dplyr
#' @importFrom crayon yellow green
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export calculateEffectSizes

calculateEffectSizes = function(data,
                                funcs.g = list(g.m.sd = g.m.sd,
                                               g.change.m.sd = g.change.m.sd,
                                               g.binary = g.binary,
                                               g.precalc = g.precalc),
                                funcs.rr = list(rr.binary = rr.binary,
                                                rr.precalc = rr.precalc),
                                include.switched.arms = FALSE,
                                change.sign = NULL,
                                impute.response = FALSE,
                                vars.for.id = c("study", "outcome_type",
                                                "instrument", "time",
                                                "time_weeks",
                                                "rating"),
                                .condition = "condition",
                                .condition.specification = "multi",
                                .groups.column.indicator = c("_arm1", "_arm2"),
                                .trt.indicator = "arm",
                                .impute.response.vars = c(m.trt.pre = "baseline_m_arm1", 
                                                          m.trt.post = "mean_arm1", 
                                                          sd.trt.post = "sd_arm1", 
                                                          n.trt = "n_arm1",
                                                          m.ctr.pre = "baseline_m_arm2", 
                                                          m.ctr.post = "mean_arm2", 
                                                          sd.ctr.post = "sd_arm2", 
                                                          n.ctr = "n_arm2")){

  
  # Check for data conflicts flagged in checkConflicts
  if ("checkConflicts" %in% class(data)){
    stop("Data format conflicts were detected. Run 'checkConflicts' to diagnose potential problems.")
  }
  
  # Remove metapsyTools variables, if they exist
  within(data, {
    .g = NULL
    .g_se = NULL
    .log_rr = NULL
    .log_rr_se = NULL
    .event_arm1 = NULL
    .event_arm2 = NULL
    .totaln_arm1 = NULL
    .totaln_arm2 = NULL
  }) -> data
  
  # Convert to data.frame (conventional)
  data = data.frame(data)

  # Get id variables
  data %>% 
    dplyr::select(
      all_of(c(vars.for.id, 
               paste0(.condition,
                      .groups.column.indicator),
               paste0(.condition.specification,
                      .groups.column.indicator)))) %>% 
    colnames() -> id.vars
  
  # Create comparison IDs
  apply(data, 1,
        function(x){paste(as.character(x[id.vars]), 
                          collapse = "_")}) %>%
    stringr::str_replace_all(",", "_") %>% 
    stringr::str_remove_all(" ") -> data$.id

  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                 #
  #   1. Hedges' g                                                  #
  #                                                                 #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  # Loop through ES functions
  data.wide = data
  es.res = list()
  for (i in 1:length(funcs.g)){
    es.res[[i]] = try({funcs.g[[i]](data.wide)}, silent = TRUE)
  }
  
  # Check if all functions could be applied
  es.res %>% purrr::map(~class(.)) %>%
    unlist() %>% {. == "try-error"} -> error.mask
  
  es.res = do.call(cbind, es.res[!error.mask])
  
  if (sum(error.mask) > 0){
    message("- ", crayon::yellow("[!] "), 
            "Function(s) ", paste(names(funcs.g)[error.mask], collapse = ", "),
            " not applied. Check for potential data/function problems.")
    message("- ", crayon::yellow("[!] "), 
            "All other Hedges' g calculation functions were applied successfully.")
  } else {
    message("- ", crayon::green("[OK] "), 
            "Hedges' g calculated successfully.")
  }
  
  # Now, bind all calculated ES together,
  # then bind together with wide dataset
  es.res %>%
    extractG() %>% 
    cbind(data.wide, .) -> dat.final
  
  colnames(dat.final)[c(ncol(dat.final)-1, 
                        ncol(dat.final))] = c(".g", ".g_se")
  
  # Change sign
  if (!is.null(change.sign)){
    change.mask = ifelse(as.logical(dat.final[[change.sign]]), -1, 1)
    dat.final$.g = dat.final$.g * change.mask
  }
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                 #
  #   2. Risk Ratio                                                 #
  #                                                                 #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  # Impute response
  if (impute.response[1]) {
    if (is.null(change.sign)) {data$change.sign = FALSE; change.sign = "change.sign"}
    purrr::map2_dfr(split(data, seq(nrow(data))), 
                    as.list(as.logical(data[[change.sign]])),
                    function(x,y){
                      imputeResponse(m.trt.pre = x[[.impute.response.vars["m.trt.pre"]]], 
                                     m.trt.post = x[[.impute.response.vars["m.trt.post"]]],
                                     sd.trt.post = x[[.impute.response.vars["sd.trt.post"]]], 
                                     n.trt = x[[.impute.response.vars["n.trt"]]],
                                     lower.is.better = ifelse(y,FALSE,TRUE)) %>% 
                        suppressMessages()}) -> imp.response.trt
    purrr::map2_dfr(split(data, seq(nrow(data))), 
                    as.list(as.logical(data[[change.sign]])),
                    function(x,y){
                      imputeResponse(m.trt.pre = x[[.impute.response.vars["m.ctr.pre"]]], 
                                     m.trt.post = x[[.impute.response.vars["m.ctr.post"]]],
                                     sd.trt.post = x[[.impute.response.vars["sd.ctr.post"]]], 
                                     n.trt = x[[.impute.response.vars["n.ctr"]]],
                                     lower.is.better = ifelse(y,FALSE,TRUE)) %>% 
                        suppressMessages()}) -> imp.response.ctr
    within(data, {
      is.na(event_arm1) -> mask_arm1; is.na(event_arm2) -> mask_arm2
      event_arm1[mask_arm1] = imp.response.trt$trtResponder[mask_arm1]
      totaln_arm1[mask_arm1] = imp.response.trt$nTrt[mask_arm1]
      event_arm2[mask_arm2] = imp.response.ctr$trtResponder[mask_arm2]
      totaln_arm2[mask_arm2] = imp.response.ctr$nTrt[mask_arm2]
    }) -> data 
  }
  
  # Loop through ES functions
  data.wide = data
  es.res = list()
  for (i in 1:length(funcs.rr)){
    es.res[[i]] = try({funcs.rr[[i]](data.wide)}, silent = TRUE)
  }
  
  # Check if all functions could be applied
  es.res %>% purrr::map(~class(.)) %>%
    unlist() %>% {. == "try-error"} -> error.mask

  es.res = do.call(cbind, es.res[!error.mask])
  
  if (sum(error.mask) > 0){
    message("- ", crayon::yellow("[!] "), 
            "Function(s) ", paste(names(funcs.rr)[error.mask], collapse = ", "),
            " not applied. Check for potential data/function problems.")
    message("- ", crayon::yellow("[!] "), 
            "All other log-risk ratio calculation functions were applied successfully.")
  } else {
    message("- ", crayon::green("[OK] "), 
            "Log-risk ratios calculated successfully.")
  }
  
  # Now, bind all calculated ES together,
  # then bind together with wide dataset
  es.res[colnames(es.res) %in% c("es", "se")] %>%
    apply(., 1, function(x){
      if (length(x[!isNAorNaN(x)]) <= 1){return(c(NA, NA))} else {
        return(x[!isNAorNaN(x)])}}) %>% t() %>%
    cbind(dat.final, .) -> dat.final
  
  colnames(dat.final)[c(ncol(dat.final)-1, 
                        ncol(dat.final))] = c(".log_rr", ".log_rr_se")
  
  # Add cell counts for Mantel-Haenszel Method (if possible)
  if (".event_arm1" %in% colnames(es.res)){
    cbind(dat.final, 
          es.res[c(".event_arm1", ".event_arm2",
                   ".totaln_arm1", ".totaln_arm2")]) -> dat.final
    dat.final[is.na(dat.final[[".log_rr"]]),
              c(".event_arm1", ".event_arm2",
                ".totaln_arm1", ".totaln_arm2")] = NA
  }
  
  
  # Check if required Metapsy variables are missing;
  # Add as NA if necessary
  mp.standard.vars = c(".id", ".g", ".g_se", 
                       ".log_rr", ".log_rr_se", ".event_arm1", 
                       ".event_arm2", ".totaln_arm1", 
                       ".totaln_arm2") 

  mp.standard.vars[!mp.standard.vars %in% 
                     colnames(dat.final)] %>% 
    {matrix(NA, nrow(dat.final), length(.),
      dimnames = list(NULL, .))} %>% 
    cbind(dat.final, .) -> dat.final
  
  # Set Inf values to NA
  dat.final[mp.standard.vars][dat.final[mp.standard.vars] == Inf] = NA
  dat.final[mp.standard.vars][dat.final[mp.standard.vars] == -Inf] = NA

  with(dat.final, {
    .event_arm1 > .totaln_arm1 |
      .event_arm2 > .totaln_arm2
  }) -> totaln.mask
  
  dat.final[!is.na(totaln.mask) & totaln.mask, 
            mp.standard.vars[-1]] = NA
  
  # Set effect sizes with variance 0 to NA
  dat.final[dat.final[[".g_se"]] == 0 & 
              !is.na(dat.final[[".g_se"]]), ".g"] = NA
  dat.final[dat.final[[".g_se"]] == 0 & 
              !is.na(dat.final[[".g_se"]]), ".g_se"] = NA
  dat.final[dat.final[[".log_rr_se"]] == 0 & 
              !is.na(dat.final[[".log_rr_se"]]), ".log_rr"] = NA
  dat.final[dat.final[[".log_rr_se"]] == 0 &
              !is.na(dat.final[[".log_rr_se"]]), ".log_rr_se"] = NA
  
  # When one variables (es, se is NA, set all to NA)
  dat.final[is.na(dat.final[[".g_se"]]), ".g"] = NA
  dat.final[is.na(dat.final[[".g"]]), ".g_se"] = NA
  dat.final[is.na(dat.final[[".log_rr_se"]]), ".log_rr"] = NA
  dat.final[is.na(dat.final[[".rr_se"]]), ".log_rr_se"] = NA
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                 #
  #   3. Add switched arms                                          #
  #                                                                 #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  # Add switched reference arm rows
  if (include.switched.arms[1] == TRUE){
    if (.trt.indicator != "arm"){
      stop("'.trt.indicator' must be set to 'arm'",
           " when including switched trial arms.")
    } else {
      # Remove suffix for _arm variables that only appear once
      dat.final %>%
        dplyr::select(dplyr::ends_with(c("_arm1", "_arm2"))) %>%
        colnames() %>%
        stringr::str_replace_all("_arm(1|2)", "") %>%
        table() %>% {names(.[.==1])} -> one.vars
      
      dat.final %>%
        dplyr::select(dplyr::starts_with(one.vars)) %>%
        colnames() -> select.vars
      
      colnames(dat.final)[colnames(dat.final) ==
                            select.vars] = one.vars
      
      dat.final = includeSwitchedArms(dat.final)
    }}
  
  # Return
  return(dat.final)

}









