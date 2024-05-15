#' Calculate Hedges' g using means and standard deviations
#'
#' Calculate Hedges' g using Mean and Standard Deviation. Only meant to be used as
#' part of \code{\link{calculateEffectSizes}}.
#'
#' @param x data
#' @param ... Effect size data. Data frame must include columns `mean_arm1`, `mean_arm2`, 
#' `sd_arm1`, `sd_arm2`, `n_arm1`, `n_arm2`. See the [Metapsy data standard](https://docs.metapsy.org/data-preparation/format/).
#' @usage g.m.sd(x, ...)
#' @importFrom dplyr select
#' @importFrom purrr pmap_dfr
#' @importFrom esc esc_mean_sd
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export g.m.sd
#' @keywords internal

g.m.sd = function(x, ...){
  x %>%
    purrr::pmap_dfr(function(mean_arm1, mean_arm2, sd_arm1,
                             sd_arm2, n_arm1, n_arm2, ...)
    {esc::esc_mean_sd(mean_arm1, sd_arm1, n_arm1,
                      mean_arm2, sd_arm2, n_arm2, es.type = "g") %>%
        as.data.frame() %>% dplyr::select(es, se) %>%
        suppressWarnings()})
}

#' Calculate Hedges' g using binary outcome data
#'
#' Calculates Hedges' g from binary outcome data. Only meant to be used as
#' part of \code{\link{calculateEffectSizes}}. 
#'
#' @param x data
#' @param cc Should a continuity correction for zero cells be applied? 
#' Either `FALSE` or the increment to be added. Default is 0.5.
#' @param ... Binary effect size data. Data frame must include columns `event_arm1`, `event_arm2`, 
#' `totaln_arm1`, `totaln_arm2`. See the [Metapsy data standard](https://docs.metapsy.org/data-preparation/format/).
#' @usage g.binary(x, cc = 0.5, ...)
#' @importFrom dplyr select mutate
#' @importFrom purrr pmap_dfr
#' @importFrom esc esc_2x2
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export g.binary
#' @keywords internal

g.binary = function(x, cc = 0.5, ...){
  
  if (identical(cc, FALSE)){
    x %>%
      purrr::pmap_dfr(function(event_arm1, event_arm2, 
                               totaln_arm1, totaln_arm2, ...)
      {esc::esc_2x2(event_arm1,
                    totaln_arm1 - event_arm1,
                    event_arm2,
                    totaln_arm2 - event_arm2,
                    es.type = "g") %>%
          as.data.frame() %>% dplyr::select(es, se) %>%
          suppressWarnings() %>%
          dplyr::mutate(es = es*-1)})
  } else {
    x %>%
      purrr::pmap_dfr(function(event_arm1, event_arm2, 
                               totaln_arm1, totaln_arm2, ...)
      {
        if (identical(event_arm1, 0) ||
            identical(event_arm2, 0)) {
          esc::esc_2x2(event_arm1 + cc,
                       totaln_arm1 - event_arm1 + cc,
                       event_arm2 + cc,
                       totaln_arm2 - event_arm2 + cc,
                       es.type = "g") %>%
            as.data.frame() %>% dplyr::select(es, se) %>%
            suppressWarnings() %>%
            dplyr::mutate(es = es*-1)
        } else {
          esc::esc_2x2(event_arm1,
                       totaln_arm1 - event_arm1,
                       event_arm2,
                       totaln_arm2 - event_arm2,
                       es.type = "g") %>%
            as.data.frame() %>% dplyr::select(es, se) %>%
            suppressWarnings() %>%
            dplyr::mutate(es = es*-1)
        }
      })
  }
}


#' Forward pre-calculated values of Hedges' g
#'
#' Forwards pre-calculated values of Hedges' g. Only meant to be used as
#' part of \code{\link{calculateEffectSizes}}. 
#'
#' @param x data
#' @param ... Pre-calculated effect size data. Data frame must include columns \code{precalc_g} and 
#' \code{precalc_g_se}. See the [Metapsy data standard](https://docs.metapsy.org/data-preparation/format/).
#' @usage g.precalc(x, ...)
#' @importFrom purrr pmap_dfr
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export g.precalc
#' @keywords internal

g.precalc = function(x, ...){
  if (is.null(x$precalc_g) | is.null(x$precalc_g_se)){
    return(data.frame(es = NA, se = NA))
  }
  x %>%
    purrr::pmap_dfr(function(precalc_g, precalc_g_se, ...)
    {data.frame(es = precalc_g, se = precalc_g_se) %>%
        suppressWarnings()})
}


#' Forward pre-calculated log-risk ratios
#'
#' Forwards pre-calculated log-risk ratios. Only meant to be used as
#' part of \code{\link{calculateEffectSizes}}.
#'
#' @param x data
#' @param ... Pre-calculated effect size data. Data frame must include columns  
#' \code{precalc_log_rr} and \code{precalc_log_rr_se}. See the [Metapsy data standard](https://docs.metapsy.org/data-preparation/format/).
#' @usage rr.precalc(x, ...)
#' @importFrom purrr pmap_dfr
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export rr.precalc
#' @keywords internal

rr.precalc = function(x, ...){
  if (is.null(x$precalc_log_rr) | is.null(x$precalc_log_rr_se)){
    return(data.frame(es = NA, se = NA))
  }
  x %>%
    purrr::pmap_dfr(function(precalc_log_rr, precalc_log_rr_se, ...)
    {data.frame(es = precalc_log_rr, se = precalc_log_rr_se) %>%
        suppressWarnings()})
}



#' Calculate the log-risk ratio using binary outcome data
#'
#' Calculate the log risk ratio using binary outcome data. Only meant to be used as
#' part of \code{\link{calculateEffectSizes}}.
#'
#' @param x data
#' @param cc Should a continuity correction for zero cells be applied? 
#' Either `FALSE` or the increment to be added. Default is 0.5.
#' @param ... Binary effect size data. Data frame must include columns `event_arm1`, `event_arm2`, 
#' `totaln_arm1`, `totaln_arm2`. See the [Metapsy data standard](https://docs.metapsy.org/data-preparation/format/).
#' @usage rr.binary(x, cc = 0.5, ...)
#' @importFrom dplyr select mutate
#' @importFrom purrr pmap_dfr
#' @importFrom esc esc_2x2
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export rr.binary
#' @keywords internal

rr.binary = function(x, cc = 0.5, ...){
  
  if (identical(cc, FALSE)){
    x %>%
      purrr::pmap_dfr(function(event_arm1, event_arm2, 
                               totaln_arm1, totaln_arm2, ...)
      {data.frame(es = log((event_arm1/totaln_arm1)/
                             (event_arm2/totaln_arm2)),
                  se = sqrt((1/event_arm1) + (1/event_arm2) - 
                              (1/totaln_arm1) - (1/totaln_arm2)),
                  .event_arm1 = event_arm1,
                  .event_arm2 = event_arm2,
                  .totaln_arm1 = totaln_arm1,
                  .totaln_arm2 = totaln_arm2) %>% 
          suppressWarnings()})
  } else {
    x %>%
      purrr::pmap_dfr(function(event_arm1, event_arm2, 
                               totaln_arm1, totaln_arm2, ...)
      {
        if (identical(event_arm1, 0) ||
            identical(event_arm2, 0)) {
          data.frame(es = log(((event_arm1 + cc)/(totaln_arm1 + cc))/
                                ((event_arm2 + cc)/(totaln_arm2 + cc))),
                     se = sqrt((1/(event_arm1 + cc)) + 
                                 (1/(event_arm2 + cc)) - 
                                 (1/(totaln_arm1 + cc)) - 
                                 (1/(totaln_arm2 + cc))),
                     .event_arm1 = event_arm1,
                     .event_arm2 = event_arm2,
                     .totaln_arm1 = totaln_arm1,
                     .totaln_arm2 = totaln_arm2) %>% 
            suppressWarnings()
        } else {
          data.frame(es = log((event_arm1/totaln_arm1)/
                                (event_arm2/totaln_arm2)),
                     se = sqrt((1/event_arm1) + (1/event_arm2) - 
                                 (1/totaln_arm1) - (1/totaln_arm2)),
                     .event_arm1 = event_arm1,
                     .event_arm2 = event_arm2,
                     .totaln_arm1 = totaln_arm1,
                     .totaln_arm2 = totaln_arm2) %>% 
            suppressWarnings()
        }
      })
  }
}


#' Calculate Hedges' g using within-group change data
#'
#' Calculate Hedges' g based on change data. Only meant to be used as
#' part of \code{\link{calculateEffectSizes}}. 
#'
#' @param x data
#' @param ... Change score effect size data. Data frame must include columns `mean_change_arm1`, 
#' `mean_change_arm2`, `sd_change_arm1`, `sd_change_arm2`, 
#' `n_change_arm1`, `n_change_arm2`. See the [Metapsy data standard](https://docs.metapsy.org/data-preparation/format/).
#' @usage g.change.m.sd(x, ...)
#' @importFrom dplyr select
#' @importFrom purrr pmap_dfr
#' @importFrom esc esc_mean_sd
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export g.change.m.sd
#' @keywords internal

g.change.m.sd = function(x, ...){
  x %>%
    purrr::pmap_dfr(function(mean_change_arm1, mean_change_arm2, sd_change_arm1,
                             sd_change_arm2, n_change_arm1, n_change_arm2, ...)
    {esc::esc_mean_sd(mean_change_arm1, sd_change_arm1, n_change_arm1,
                      mean_change_arm2, sd_change_arm2, n_change_arm2, 
                      es.type = "g") %>%
        as.data.frame() %>% dplyr::select(es, se) %>%
        suppressWarnings()})
}


#' Calculate NNTs
#'
#' Calculate NNTs (extracted from \code{dmetar})
#'
#' @param d A single numeric or concatenated vector of numerics representing the effect size expressed as
#' Cohen's \eqn{d} or Hedges' \eqn{g}. If this is the only parameter specified in the function, the method by
#' Kraemer and Kupfer is used automatically to calculate \eqn{NNT}s.
#' @param CER The control group event ratio. Furukawa's method (Furukawa & Leucht, 2011) to calculate \code{NNT}s
#' from \code{d} requires that the assumed response ("event") ratio in the control group (\eqn{\frac{n_{responders}}{N_{total}}})
#' is specified. The CER can assume values from 0 to 1. If a value is specified for \code{CER}, Furukawa's method is
#' used automatically. Argument \code{method} has to be set to \code{"KraemerKupfer"} to override this.
#' @param event.e Single number or numeric vector. The number of (favourable) events in the experimental group.
#' @param n.e Single number or numeric vector. The number participants in the experimental group.
#' @param event.c Single number or numeric vector. The number of (favourable) events in the control group.
#' @param n.c Single number or numeric vector. The number of participants in the control group.
#' @param names Optional. Character vector of equal length as the vector supplied to \code{d} or \code{event.e} containing
#' study/effect size labels.
#' @param method The method to be used to calculate the NNT from \code{d}. Either \code{"KraemerKupfer"} for the
#' method proposed by Kraemer and Kupfer (2006) or \code{"Furukawa"} for the Furukawa method (Furukawa & Leucht, 2011).
#' Please note that the Furukawa's method can only be used when \code{CER} is specified.
#' @usage metapsyNNT(d, CER, event.e, n.e, event.c, n.c, names, method)
#' @export metapsyNNT
#' @keywords internal

metapsyNNT = function(d, CER = NULL, event.e = NULL, n.e = NULL, event.c = NULL,
          n.c = NULL, names = NULL, method = NULL)
{
  if (!missing(CER) & is.numeric(CER)) {
    if (sum(CER < 0 | CER > 1) > 0) {
      stop("'CER' range must be between 0 and 1.")
    }
  }
  if (missing(event.e)) {
    if (missing(CER)) {
      NNT = 1/((2 * pnorm(d/sqrt(2)) - 1))
      class(NNT) = c("NNT", "kk", "numeric")
    }
    else {
      if (missing(method)) {
        NNT = 1/(pnorm(d + qnorm(CER)) - CER)
        class(NNT) = c("NNT", "fl", "numeric")
      }
      else {
        if (method == "KraemerKupfer") {
          NNT = 1/((2 * pnorm(d/sqrt(2)) - 1))
          class(NNT) = c("NNT", "kk", "numeric")
        }
        else {
        }
        if (method == "Furukawa") {
          if (missing(CER) | class(CER) != "numeric") {
            stop("To use Furukawa's method, provide a numeric value for CER. \n")
          }
          else {
            NNT = 1/(pnorm(d + qnorm(CER)) - CER)
            class(NNT) = c("NNT", "fl", "numeric")
          }
        }
        else {
        }
      }
    }
  }
  else {
    if (class(event.e) == "numeric" & missing(d)) {
      EER = event.e/n.e
      CER = event.c/n.c
      NNT = abs(1/(EER - CER))
      class(NNT) = c("NNT", "raw", "numeric")
      if (missing(method)) {
      }
      else {
        if (method %in% c("KraemerKupfer", "Furukawa")) {
          class(NNT) = c("NNT", "unspecified", "numeric")
        }
        else {
        }
      }
    }
  }
  if (missing(names)) {
    return(NNT)
  }
  else {
    data = data.frame(Name = names, NNT = NNT)
    class(data) = c(class(NNT)[1:2], "data.frame")
    class(data$NNT) = "numeric"
    NNT = data
    return(NNT)
  }
}


#' Expander function for multiarm trials
#'
#' Expands multiarm trial.
#'
#' @param study data of one study, for one unique assessment point.
#' @param condition.specification trial condition specification.
#' @param group.indicator group indicator (IG or CG).
#' @param group.names group names (list).
#' @export multiarmExpander
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @keywords internal

multiarmExpander = function(study,
                            condition.specification,
                            group.indicator,
                            group.names){

  # Replace missings with empty char
  study[[condition.specification]][is.na(study[[condition.specification]])] = "N/A"


  # Get combination of study arms
  varType = study[[group.indicator]]
  names(varType) = study[[condition.specification]]
  utils::combn(study[[condition.specification]], 2,
        simp = FALSE) -> combinations

  # Expand by looping through all combinations;
  # IGs come first if compared to CGs
  trialExpand = list()

  for (i in 1:length(combinations)){

    if (varType[[combinations[[i]][1]]] == group.names[["cg"]] &
        varType[[combinations[[i]][2]]] == group.names[["ig"]]){

      rbind(study[study[[condition.specification]] == combinations[[i]][2],],
            study[study[[condition.specification]] == combinations[[i]][1],]) -> trial

    } else {

      rbind(study[study[[condition.specification]] == combinations[[i]][1],],
            study[study[[condition.specification]] == combinations[[i]][2],]) -> trial

    }

    trial$trt = c("trt1", "trt2")
    trial$id = paste0(trial$id, "_", paste(combinations[[i]], collapse = "_"))
    trialExpand[[i]] = trial
  }

  do.call(rbind, trialExpand) -> res

  return(res)
}


#' String detection wrapper
#'
#' Wrapper for \code{str_detect} in the \code{stringr} package.
#'
#' @param ... Arguments passed to \code{str_detect}.
#' @usage Detect(...)
#' @export Detect
#' @importFrom stringr str_detect
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @keywords internal

Detect = function(...){
  stringr::str_detect(...)
}



#' Include information of rows with switched reference group
#'
#' Adds effect size data and study information of rows with switched reference arms.
#'
#' @param dat Data set created by \code{calculateEffectSizes} in the wider format, which includes
#' calculated effect sizes and standard errors in columns \code{es} and \code{se}, respectively. Only
#' data sets created by \code{calculateEffectSizes} with \code{trt.indicator} set to \code{"trt"}
#' can be used.
#' @param ... Further arguments (not used).
#' @usage includeSwitchedArms(dat, ...)
#' @export includeSwitchedArms
#' @importFrom stringr str_replace_all
#' @importFrom dplyr mutate select
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @keywords internal

includeSwitchedArms = function(dat, ...){
  dat.orig = dat
  stringr::str_replace_all(colnames(dat),
                           c("_arm1" = "_armX", "_arm2" = "_arm1",
                             "_armX" = "_arm2")) -> colnames(dat)

  dat %>% dplyr::select(colnames(dat.orig)) %>%
    dplyr::mutate(.g = .g*-1,
                  .log_rr = .log_rr*-1,
                  .id = paste0(.id, "_arm_switched")) %>%
    rbind(., dat.orig) -> dat.expand

  dat.expand = dat.expand[order(dat.expand$.id),]
  rownames(dat.expand) = NULL

  return(dat.expand)
}


#' `tryCatch` alternative that saves the error message
#'
#' @param expr R expression.
#' @usage tryCatch2(expr)
#' @export 
#' @keywords internal
tryCatch2 = function(expr) {
  warn = err = NULL
  value = withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, 
       error=err, has.error = is.null(value))
}


#' Domain 1 RoB Algorithm
#' @keywords internal 
domain1Algorithm = function(d1_1, d1_2, d1_3, d1_4) {
  if (identical(d1_2, "No")) {"high risk"}
  else if (identical(d1_2, "NI") & (identical(d1_3, "No")|identical(d1_4, "No"))) {"high risk"}
  else if (identical(d1_2, "NI") & ((identical(d1_3, "Yes")|identical(d1_3, "NI"))|(identical(d1_4, "Yes")|identical(d1_4, "NI")))) {"some concerns"}
  else if (identical(d1_2, "Yes") & identical(d1_1, "No")) {"some concerns"}
  else if (identical(d1_2, "Yes") & (identical(d1_1, "Yes")|identical(d1_1, "NI")) & (identical(d1_3, "No")|identical(d1_4, "No"))) {"some concerns"}
  else if (identical(d1_2, "Yes") & (identical(d1_1, "Yes")|identical(d1_1, "NI")) & ((identical(d1_3, "Yes")|identical(d1_3, "NI"))|(identical(d1_4, "Yes")|identical(d1_4, "NI")))) {"low risk"}
  else NA
}

#' Domain 2 RoB Algorithm
#' @keywords internal 
domain2Algorithm = function(d2_5, d2_6, d2_7, d2_8, d2_9) {
  if (identical(d2_5, "Yes") & identical(d2_6, "Yes")) {"low risk"}
  else if ((identical(d2_7, "NI")|identical(d2_7, "Yes")) & identical(d2_8, "Yes")) {"low risk"}
  else if ((identical(d2_7, "NI")|identical(d2_7, "Yes")) & (identical(d2_8, "No")|identical(d2_8, "NI")) & (identical(d2_9, "Yes"))) {"some concerns"}
  else if ((identical(d2_7, "NI")|identical(d2_7, "Yes")) & (identical(d2_8, "No")|identical(d2_8, "NI")) & (identical(d2_9, "No")|identical(d2_9, "NI"))) {"high risk"}
  else if (identical(d2_7, "No") & (d2_8 %in% c("Yes", "No", "NI")) & (d2_9 %in% c("Yes", "No", "NI"))) {"high risk"}
  else NA
}

#' Domain 3 RoB Algorithm
#' @keywords internal 
domain3Algorithm = function(d3_10, d3_11, d3_12, d3_13, d3_14) {
  if (identical(d3_10, "Yes")) {"low risk"}
  else if ((identical(d3_11, "No")|identical(d3_11, "NI")) & identical(d3_12, "Yes") & identical(d3_13, "Yes") & identical(d3_14, "Yes")) {"some concerns"}
  else if ((identical(d3_11, "No")|identical(d3_11, "NI")) & identical(d3_12, "Yes") & (identical(d3_13, "No")|identical(d3_13, "NI")) & identical(d3_14, "Yes")) {"some concerns"}
  else if ((identical(d3_11, "No")|identical(d3_11, "NI")) & identical(d3_12, "Yes") & identical(d3_13, "Yes") & (identical(d3_14, "No")|identical(d3_14, "NI"))) {"some concerns"}
  else if ((identical(d3_11, "No")|identical(d3_11, "NI")) & identical(d3_12, "Yes") & (identical(d3_13, "No")|identical(d3_13, "NI")) & (identical(d3_14, "No")|identical(d3_14, "NI"))) {"high risk"}
  else if ((identical(d3_11, "No")|identical(d3_11, "NI")) & (identical(d3_12, "No")|identical(d3_12, "NI")) & identical(d3_13, "Yes") & identical(d3_14, "Yes")) {"some concerns"}
  else if ((identical(d3_11, "No")|identical(d3_11, "NI")) & (identical(d3_12, "No")|identical(d3_12, "NI")) & (identical(d3_13, "No")|identical(d3_13, "NI")) & identical(d3_14, "Yes")) {"high risk"}
  else if ((identical(d3_11, "No")|identical(d3_11, "NI")) & (identical(d3_12, "No")|identical(d3_12, "NI")) & identical(d3_13, "Yes") & (identical(d3_14, "No")|identical(d3_14, "NI"))) {"high risk"}
  else if ((identical(d3_11, "No")|identical(d3_11, "NI")) & (identical(d3_12, "No")|identical(d3_12, "NI")) & (identical(d3_13, "No")|identical(d3_13, "NI")) & (identical(d3_14, "No")|identical(d3_14, "NI"))) {"high risk"}
  else NA
}

#' Domain 4 RoB Algorithm
#' @keywords internal 
domain4Algorithm = function(d4_15, d4_16, d4_17, d4_18) {
  if (identical(d4_16, "Yes") & identical(d4_17, "No")) {"low risk"}
  else if (identical(d4_16, "No") & identical(d4_17, "Yes") & identical(d4_18, "Yes")) {"low risk"}
  else if (identical(d4_16, "No") & identical(d4_17, "Yes") & (identical(d4_18, "No")|identical(d4_18, "NI"))) {"high risk"}
  else if (identical(d4_16, "Yes") & identical(d4_17, "Yes") & identical(d4_18, "Yes")) {"low risk"}
  else if (identical(d4_16, "Yes") & identical(d4_17, "Yes") & (identical(d4_18, "No")|identical(d4_18, "NI"))) {"high risk"}
  else NA
}

#' Domain 5 RoB Algorithm
#' @keywords internal 
domain5Algorithm = function(d5_19, d5_20, d5_21, d5_22, d5_23, d5_24) {
  if ((identical(d5_19, "NI")|identical(d5_19, "No"))) {"some concerns"}
  else if ((identical(d5_21, "NI")|identical(d5_21, "No"))) {"some concerns"}
  else if (identical(d5_22, "Yes") & identical(d5_24, "No")) {"high risk"}
  else if (identical(d5_22, "Yes")) {"low risk"}
  else if (identical(d5_22, "NI") & identical(d5_24, "No")) {"high risk"}
  else if (identical(d5_22, "NI")) {"some concerns"}
  else if (identical(d5_22, "No")) {"high risk"}
  else NA
}

#' Overall RoB Algorithm
#' @keywords internal 
overallAlgorithm = function(d1, d2, d3, d4, d5) {
  ratings = c(d1, d2, d3, d4, d5)
  if (sum(is.na(ratings))>0) { NA } else {
    lows = sum(ratings=="low risk"); highs = sum(ratings=="high risk")
    concerns = sum(ratings=="some concerns")
    if (lows==5) {"low risk"}
    else if (highs > 0 | concerns > 2) {"high risk"}
    else if (concerns > 0 | concerns <= 2) {"some concerns"}
    else NA
  }
}


#' Tests if value is `NA` or `NaN`
#'
#' @param expr R expression.
#' @usage isNAorNaN(x)
#' @export 
#' @keywords internal
isNAorNaN = function(x){
  is.na(x) | is.nan(x)
}


#' Extract Hedges _g_ values
#'
#' @param x `data.frame` with effect size results.
#' @usage extractG(x)
#' @export 
#' @keywords internal
extractG = function(x){
  mask = apply(
    x, 1, function(f) !isNAorNaN(f), simplify = FALSE) %>% 
    do.call(rbind, .)
  index = apply(mask, 1, which, simplify = FALSE)
  lapply(as.list(1:length(index)), 
         function(i) {
           y = x[i, index[[i]]]
           if (length(y) <= 1) { 
             return(data.frame(es = NA, se = NA)[,c("es", "se")]) } 
           else {return(y[,c("es", "se")])}}) %>% 
    do.call(rbind, .) -> x.select
  return(x.select)
}


#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom tidyr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL


#' Select available arguments in function
#' @keywords internal 
selectArguments = function(f, dots){
  func.args = names(as.list(args(f)))
  return(dots[names(dots) %in% func.args])
}

#' Try function and catch errors
#' @keywords internal 
tryCatch2 = function(expr) {
  warn = err = NULL
  value = withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, 
       error=err, has.error = is.null(value))
}


#' Send message after fitting model
#' @keywords internal 
sendMessage = function(obj, model.name = NULL, 
                       which.run = NULL, .type.es = NULL, 
                       es.type = NULL){
  if (identical(obj, "start")){
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
  } else {
    if (model.name %in% which.run){
      if (!obj$has.error) {
        message(crayon::green("DONE"))
      } else {
        message(crayon::yellow("ERROR"))
        message("- ", crayon::yellow("[!] "), "'", model.name, "' ", 
                "model could not be calculated. Error: ",
                obj$message)
      }
    }
}}


#' Fit 'overall' model
#' @keywords internal 
fitOverallModel = function(data, es.var, se.var, arm.var.1, arm.var.2,
                           measure.var, study.var, .raw.bin.es, .type.es, hakn,
                           method.tau.meta, method.tau.ci, method.tau,
                           dots, es.binary.raw.vars, round.digits,
                           nnt.cer, which.run, rob.data){
  
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
    mGeneralArgs = list(
      TE = data[[es.var[1]]],
      seTE = data[[se.var[1]]],
      studlab = data[[study.var]],
      data = data,
      sm = .type.es,
      hakn = hakn,
      method.tau = method.tau.meta,
      method.tau.ci = method.tau.ci,
      prediction = TRUE,
      common = ifelse(method.tau == "FE", TRUE, FALSE),
      fixed = ifelse(method.tau == "FE", TRUE, FALSE),
      random = ifelse(method.tau == "FE", FALSE, TRUE)) %>% 
      append(selectArguments(meta::metagen, dots))
    
    mGeneral = 
      tryCatch2(do.call(meta::metagen, mGeneralArgs))
  } else {
    mGeneralArgs = list(
      event.e = data[[es.binary.raw.vars[1]]],
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
      common = ifelse(method.tau == "FE", TRUE, FALSE),
      fixed = ifelse(method.tau == "FE", TRUE, FALSE),
      random = ifelse(method.tau == "FE", FALSE, TRUE)) %>% 
      append(selectArguments(meta::metabin, dots))
    
    mGeneral = 
      tryCatch2(do.call(meta::metabin, mGeneralArgs))
    
    if (!mGeneral$has.error){
      with(updateMeta(mGeneral$value, sm="RD"), {
        ifelse(isTRUE(common) & !isTRUE(random),
               abs(TE.common)^-1, abs(TE.random)^-1)
      }) -> nnt.raw.bin.es
    } else {
      nnt.raw.bin.es = NULL
    }
  }
  
  # Collect results
  if (mGeneral$has.error){
    mGeneralRes = 
      data.frame(
        k = NA, g = NA, g.ci = NA, p = NA, i2 = NA, i2.ci = NA, 
        prediction.ci = NA, nnt = NA, excluded = "none")
  } else {
    mGeneralRes = with(mGeneral$value, {
      data.frame(
        k = k,
        g = ifelse(
          isTRUE(common) & !isTRUE(random), 
          TE.common, TE.random) %>% 
          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
          round(round.digits),
        g.ci = paste0("[", 
                      ifelse(isTRUE(common) & !isTRUE(random),
                             lower.common, lower.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "; ",
                      ifelse(isTRUE(common) & !isTRUE(random),
                             upper.common, upper.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "]"),
        p = ifelse(isTRUE(common) & !isTRUE(random), 
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
          ifelse(isTRUE(common) & !isTRUE(random), 
                 abs(TE.common), abs(TE.random)), nnt.cer) %>%
          ifelse(identical(.type.es, "RR"), NA, .) %>% 
          ifelse(isTRUE(.raw.bin.es), 
                 nnt.raw.bin.es, .) %>% 
          round(round.digits) %>% abs(),
        excluded = "none"
      )
    })
  }
  
  # Add rob data if specified
  if (!is.null(rob.data[1])) {
    mGeneral$value = addRobData(mGeneral$value, rob.data)
  }
  
  return(list(m = mGeneral$value, 
              res = mGeneralRes,
              has.error = mGeneral$has.error,
              message = mGeneral$error$message))
}



#' Fit 'lowest' model
#' @keywords internal 
fitLowestModel = function(data, study.var, multi.study,
                          mGeneral, .type.es, round.digits,
                          .raw.bin.es, nnt.cer, rob.data){
  
  message("- ", crayon::green("[OK] "), 
          "Calculating effect size using only lowest effect... ",
          appendLF = FALSE)
  
  multi.study = names(table(data[[study.var]])
                      [table(data[[study.var]]) > 1])
  
  if (length(multi.study) > 0){
    if (length(multi.study) > 0){
      data$.TE = mGeneral$m$TE
      data %>%
        split(.[[study.var]]) %>% 
        {.[unique(data[[study.var]])]} %>% 
        purrr::map(function(x){
          tiebreaker = rnorm(nrow(x), sd=1e-10)
          x$.TE + tiebreaker == suppressWarnings(min(x$.TE + tiebreaker))}) %>%
        do.call(c,.) -> lowest
      data$.TE = NULL
    } else {
      lowest = NULL
    }
    mLowest = {
      mLowest = tryCatch2(
        updateMeta(mGeneral$m, exclude = !lowest, id = NULL))
      if (!mLowest$has.error){
        with(updateMeta(mLowest$value, sm="RD"), {
          ifelse(isTRUE(common) & !isTRUE(random),
                 abs(TE.common)^-1, abs(TE.random)^-1)
        }) -> nnt.raw.bin.es
      }
      # Collect results
      if (mLowest$has.error){
        mLowestRes = 
          data.frame(
            k = NA, g = NA, g.ci = NA, p = NA, i2 = NA, i2.ci = NA, 
            prediction.ci = NA, nnt = NA, excluded = "none")
      } else {
        mLowestRes = with(mLowest$value,{
          data.frame(
            k = k,
            g = ifelse(
              isTRUE(common) & !isTRUE(random), 
              TE.common, TE.random) %>% 
              ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
              round(round.digits),
            g.ci = paste0("[", 
                          ifelse(isTRUE(common) & !isTRUE(random),
                                 lower.common, lower.random) %>% 
                            ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                            round(round.digits), "; ",
                          ifelse(isTRUE(common) & !isTRUE(random),
                                 upper.common, upper.random) %>% 
                            ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                            round(round.digits), "]"),
            p = ifelse(isTRUE(common) & !isTRUE(random), 
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
              ifelse(isTRUE(common) & !isTRUE(random), 
                     abs(TE.common), abs(TE.random)), nnt.cer) %>%
              ifelse(identical(.type.es, "RR"), NA, .) %>% 
              ifelse(isTRUE(.raw.bin.es), 
                     nnt.raw.bin.es, .) %>% 
              round(round.digits) %>% abs(),
            excluded = "none")
        })
      }
      mLowest = list(m = mLowest$value, 
                     res = mLowestRes,
                     has.error = mLowest$has.error,
                     message = mLowest$error$message)
    }
  } else {
    mLowest = mGeneral
  }
  
  if (is.null(mLowest$m$exclude[1])){
    mLowest$res$excluded = "none"
  } else {
    mLowest$res$excluded = paste(
      mLowest$m$data$comparison[mLowest$m$exclude],
      collapse = "; ")
  }
  
  # Add rob data if specified
  if (!is.null(rob.data[1])) {
    mLowest$value = addRobData(mLowest$m, rob.data)
  }
  
  return(mLowest)
}


#' Fit 'highest' model
#' @keywords internal 
fitHighestModel = function(data, study.var, multi.study,
                           mGeneral, .type.es, round.digits,
                           .raw.bin.es, nnt.cer, rob.data){
  
  message("- ", crayon::green("[OK] "), 
          "Calculating effect size using only highest effect... ",
          appendLF = FALSE)
  
  multi.study = names(
    table(data[[study.var]])
    [table(data[[study.var]]) > 1])
  
  if (length(multi.study) > 0){
    if (length(multi.study) > 0){
      data$.TE = mGeneral$m$TE
      data %>%
        split(.[[study.var]]) %>%
        {.[unique(data[[study.var]])]} %>% 
        purrr::map(function(x){
          tiebreaker = rnorm(nrow(x), sd=1e-10)
          x$.TE + tiebreaker == 
            suppressWarnings(max(x$.TE + tiebreaker))}) %>%
        do.call(c,.) -> highest
      data$.TE = NULL
    } else {
      highest = NULL
    }
    # Fit model
    mHighest = {
      
      mHighest = tryCatch2(
        updateMeta(mGeneral$m, exclude = !highest, id = NULL))
      
      if (!mHighest$has.error){
        with(updateMeta(mHighest$value, sm="RD"), {
          ifelse(isTRUE(common) & !isTRUE(random),
                 abs(TE.common)^-1, abs(TE.random)^-1)
        }) -> nnt.raw.bin.es
      }
      
      # Collect results
      if (mHighest$has.error){
        mHighestRes = 
          data.frame(
            k = NA, g = NA, g.ci = NA, p = NA, i2 = NA, i2.ci = NA, 
            prediction.ci = NA, nnt = NA, excluded = "none")
      } else {
        mHighestRes = with(mHighest$value,{
          data.frame(
            k = k,
            g = ifelse(
              isTRUE(common) & !isTRUE(random), 
              TE.common, TE.random) %>% 
              ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
              round(round.digits),
            g.ci = paste0("[", 
                          ifelse(isTRUE(common) & !isTRUE(random),
                                 lower.common, lower.random) %>% 
                            ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                            round(round.digits), "; ",
                          ifelse(isTRUE(common) & !isTRUE(random),
                                 upper.common, upper.random) %>% 
                            ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                            round(round.digits), "]"),
            p = ifelse(isTRUE(common) & !isTRUE(random), 
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
              ifelse(isTRUE(common) & !isTRUE(random), 
                     abs(TE.common), abs(TE.random)), nnt.cer) %>%
              ifelse(identical(.type.es, "RR"), NA, .) %>% 
              ifelse(isTRUE(.raw.bin.es), 
                     nnt.raw.bin.es, .) %>% 
              round(round.digits) %>% abs(),
            excluded = "none")
        })
      }
      list(m = mHighest$value, 
           res = mHighestRes,
           has.error = mHighest$has.error,
           message = mHighest$error$message)
    }
  } else {
    mHighest = mGeneral
  }
  
  if (is.null(mHighest$m$exclude[1])){
    mHighest$res$excluded = "none"
  } else {
    mHighest$res$excluded = 
      paste(mHighest$m$data$comparison[mHighest$m$exclude],
            collapse = "; ")
  }
  
  # Add rob data if specified
  if (!is.null(rob.data[1])) {
    mHighest$value = addRobData(mHighest$m, rob.data)
  }
  
  return(mHighest)
}


#' Fit 'combined' model
#' @keywords internal 
fitCombinedModel = function(which.combine, which.combine.var,
                            data, study.var, multi.study, es.var, se.var,
                            mGeneral, .type.es, round.digits, hakn,
                            .raw.bin.es, nnt.cer, rho.within.study,
                            method.tau, method.tau.ci, dots,
                            es.binary.raw.vars, phi.within.study,
                            rob.data){
  
  message("- ", crayon::green("[OK] "), 
          "Calculating effect size using combined effects (rho=", 
          rho.within.study, 
          appendLF = FALSE)
  
  # Select methods
  method.tau.ci = ifelse(identical(method.tau.ci,"Q-Profile"), 
                         "QP", method.tau.ci)
  method.tau.meta = ifelse(identical(method.tau[1], "FE"),
                           "REML", method.tau[1])
  
  # Aggregate arm-wise or trial-wise
  if (identical(which.combine[1], "arms")){
    if (!is.null(which.combine.var)){
      data$study.var.comb = 
        paste(data[[study.var]],
              ifelse(is.na(data[[which.combine.var]]), "", 
                     data[[which.combine.var]]))
      study.var.comb = "study.var.comb"
    } else {
      data$study.var.comb = data[[study.var]]
      study.var.comb = study.var
    }
  }
  if (identical(which.combine[1], "studies")){
    data$study.var.comb = data[[study.var]]
    study.var.comb = study.var
  }
  
  # Second part of message
  message("; ", ifelse(
    identical(study.var.comb[1], "study.var.comb"),
    "arm-wise", "study-wise"), ")... ",
    appendLF = FALSE)
  
  
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
  if (phi.within.study[1] >.91){
    message("- ", crayon::yellow("[!] "), 
            "'phi.within.study' is very close to 1.",
            " This can lead to non-positive-definite variance-covariance matrices.",
            " If the V matrix is not positive-definite, assume a lower value...")
    warn.end = TRUE}
  if (phi.within.study[1] <.1){
    message("- ", crayon::yellow("[!] "), 
            "'phi.within.study' is very close to 0.",
            " This can lead to non-positive-definite variance-covariance matrices.",
            " If the V matrix is not positive-definite, assume a higher value...")
    warn.end = TRUE
  }
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
  
  # Fit model
  mComb = {
    mCombArgs = list(
      TE = data.comb$yi,
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
      common = ifelse(method.tau == "FE", TRUE, FALSE),
      fixed = ifelse(method.tau == "FE", TRUE, FALSE),
      random = ifelse(method.tau == "FE", FALSE, TRUE)) %>% 
      append(selectArguments(meta::metagen, dots))
    
    mComb = tryCatch2(
      do.call(meta::metagen, mCombArgs))
    
    if (isTRUE(.raw.bin.es)){
      mCombArgsBin = list(
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
        common = ifelse(method.tau == "FE", TRUE, FALSE),
        fixed = ifelse(method.tau == "FE", TRUE, FALSE),
        random = ifelse(method.tau == "FE", FALSE, TRUE))
      
      do.call(meta::metabin, 
              mCombArgsBin) %>% 
        with(., {
          ifelse(isTRUE(common) & !isTRUE(random),
                 abs(TE.common)^-1, abs(TE.random)^-1)
        }) -> nnt.raw.bin.es
    }
    
    if (mComb$has.error){
      mCombRes = 
        data.frame(
          k = NA, g = NA, g.ci = NA, p = NA,
          i2 = NA, i2.ci = NA, prediction.ci = NA,
          nnt = NA, excluded = "none")
    } else {
      within(mComb$value$data, 
             {.event.e = data.comb[[es.binary.raw.vars[1]]]
             .n.e = data.comb[[es.binary.raw.vars[3]]]
             .event.c = data.comb[[es.binary.raw.vars[2]]]
             .n.c = data.comb[[es.binary.raw.vars[4]]]}) -> mComb$value$data
      
      mCombRes = with(mComb$value,{
        data.frame(
          k = k,
          g = ifelse(
            isTRUE(common) & !isTRUE(random), 
            TE.common, TE.random) %>% 
            ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
            round(round.digits),
          g.ci = paste0("[", 
                        ifelse(isTRUE(common) & !isTRUE(random),
                               lower.common, lower.random) %>% 
                          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                          round(round.digits), "; ",
                        ifelse(isTRUE(common) & !isTRUE(random),
                               upper.common, upper.random) %>% 
                          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                          round(round.digits), "]"),
          p = ifelse(isTRUE(common) & !isTRUE(random), 
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
            ifelse(isTRUE(common) & !isTRUE(random), 
                   abs(TE.common), abs(TE.random)), nnt.cer) %>%
            ifelse(identical(.type.es, "RR"), NA, .) %>% 
            ifelse(isTRUE(.raw.bin.es), 
                   nnt.raw.bin.es, .) %>% 
            round(round.digits) %>% abs(),
          excluded = "none")
      })
    }
    
    rownames(mCombRes) = "Combined"
    mCombRes$excluded = 
      paste(
        "combined",
        ifelse(identical(study.var.comb[1], "study.var.comb"),
               "(arm-level):", "(study-level):"),
        ifelse(length(multi.study.comb) > 0,
               paste(multi.study.comb, collapse = "; "), 
               "none"))
    
    # Add rho, which.combine, combine var to mComb model
    mComb$rho = rho.within.study
    mComb$which.combine = 
      ifelse(identical(study.var.comb[1], "study.var.comb"), "arms", "studies")
    mComb$which.combine.var = which.combine.var
    
    # Add rob data if specified
    if (!is.null(rob.data[1])) {
      mComb$value = addRobData(mComb$value, rob.data)
    }
    
    # return
    return(list(m = mComb$value, 
                res = mCombRes,
                has.error = mComb$has.error,
                message = mComb$error$message))
  }
}

#' Fit 'combined' complex model
#' @keywords internal 
fitCombinedHACEModel = function(which.combine, which.combine.var, measure.var,
                                data, study.var, multi.study, es.var, se.var,
                                mGeneral, .type.es, round.digits, hakn,
                                .raw.bin.es, nnt.cer, rho.within.study,
                                method.tau, method.tau.ci, dots,
                                es.binary.raw.vars, arm.var.1, arm.var.2,
                                phi.within.study, n.var.arm1, 
                                n.var.arm2, w1.var, w2.var, time.var,
                                near.pd, rob.data){
  
  message("- ", crayon::green("[OK] "), 
          "Calculating effect size using combined effects (rho=", 
          rho.within.study, "; phi=", phi.within.study, "; multiarm corr.",
          appendLF = FALSE)
  
  # Select methods
  method.tau.ci = ifelse(identical(method.tau.ci,"Q-Profile"), 
                         "QP", method.tau.ci)
  method.tau.meta = ifelse(identical(method.tau[1], "FE"),
                           "REML", method.tau[1])
  
  # Aggregate arm-wise or trial-wise
  if (identical(which.combine[1], "arms")){
    if (!is.null(which.combine.var)){
      data$study.var.comb = 
        paste(data[[study.var]],
              ifelse(is.na(data[[which.combine.var]]), "", 
                     data[[which.combine.var]]))
      study.var.comb = "study.var.comb"
    } else {
      data$study.var.comb = data[[study.var]]
      study.var.comb = study.var
    }
  }
  if (identical(which.combine[1], "studies")){
    data$study.var.comb = data[[study.var]]
    study.var.comb = study.var
  }
  
  # Second part of message
  message("; ", ifelse(
    identical(study.var.comb[1], "study.var.comb"),
    "arm-wise", "study-wise"), ")... ",
    appendLF = FALSE)
  
  
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
  if (phi.within.study[1] >.90){
    message("- ", crayon::yellow("[!] "), 
            "'phi.within.study' is very close to 1.",
            " This can lead to non-positive-definite variance-covariance matrices.",
            " If the V matrix is not positive-definite, assume a lower value...")
    warn.end = TRUE}
  if (phi.within.study[1] <.1){
    message("- ", crayon::yellow("[!] "), 
            "'phi.within.study' is very close to 0.",
            " This can lead to non-positive-definite variance-covariance matrices.",
            " If the V matrix is not positive-definite, assume a higher value...")
    warn.end = TRUE
  }
  if (rho.within.study[1] <.01){
    message("- ", crayon::yellow("[!] "), 
            "'rho.within.study' is very close to 0.",
            " This can lead to non-positive-definite variance-covariance matrices.",
            " If the V matrix is not positive-definite, assume a higher value...")
    warn.end = TRUE
  }
  
  data.che = data
  data.che$TE = data$yi
  data.che$seTE = sqrt(data$vi)
  
  # Define required variables
  data.che$study = data.che[[study.var]]
  data.che$instrument = data.che[[measure.var]]
  data.che$condition_arm1 = data.che[[arm.var.1]]
  data.che$condition_arm2 = data.che[[arm.var.2]]
  data.che$multi_arm1 = data.che[[which.combine.var]]
  data.che$n_arm1 = data.che[[w1.var]]
  data.che$n_arm2 = data.che[[w2.var]]
  data.che$time_weeks = data.che[[time.var]]
  within(data.che, {
    instrument[is.na(instrument)] = "missing"
    vcalc_arm1 = paste(
      study, condition_arm1, multi_arm1, "arm1")
    vcalc_arm2 = paste(
      study, condition_arm2, "arm2") 
    time_weeks[is.na(time_weeks)] = median(time_weeks, na.rm=T)
    n_arm1[is.na(n_arm1)] = 20
    n_arm2[is.na(n_arm2)] = 20
    es = TE; V = seTE^2
  }) -> dat.che
  
  # Approximate Vcovs
  metafor::vcalc(vi = V, cluster = svc, 
                 obs = instrument, time1 = time_weeks,
                 grp1 = vcalc_arm1, grp2 = vcalc_arm2,
                 w1 = n_arm1, w2 = n_arm2, 
                 rho = rho.within.study, 
                 phi = phi.within.study, data = dat.che, 
                 nearpd = near.pd) -> Vcov
  tryCatch2(
    metafor::aggregate.escalc(data, cluster = svc, V = Vcov)) -> data.comb
  
  if (data.comb$has.error){
    return(list(m = NA, 
                res = NA,
                has.error = TRUE,
                message = data.comb))
  } else {
    data.comb = data.comb$value
  }
  
  # Get studlabs of studies with multiple arms
  table(data.comb[[study.var]], 
        data.comb[[study.var.comb]]) %>%
    {names(rowSums(.)[rowSums(.) > 1])} -> multi.entries
  
  # Fit model
  mComb = {
    mCombArgs = list(
      TE = data.comb$yi,
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
      common = ifelse(method.tau == "FE", TRUE, FALSE),
      fixed = ifelse(method.tau == "FE", TRUE, FALSE),
      random = ifelse(method.tau == "FE", FALSE, TRUE)) %>% 
      append(selectArguments(meta::metagen, dots))
    
    mComb = tryCatch2(
      do.call(meta::metagen, mCombArgs))
    
    if (isTRUE(.raw.bin.es)){
      mCombArgsBin = list(
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
        common = ifelse(method.tau == "FE", TRUE, FALSE),
        fixed = ifelse(method.tau == "FE", TRUE, FALSE),
        random = ifelse(method.tau == "FE", FALSE, TRUE))
      
      do.call(meta::metabin, 
              mCombArgsBin) %>% 
        with(., {
          ifelse(isTRUE(common) & !isTRUE(random),
                 abs(TE.common)^-1, abs(TE.random)^-1)
        }) -> nnt.raw.bin.es
    }
    
    if (mComb$has.error){
      mCombRes = 
        data.frame(
          k = NA, g = NA, g.ci = NA, p = NA,
          i2 = NA, i2.ci = NA, prediction.ci = NA,
          nnt = NA, excluded = "none")
    } else {
      within(mComb$value$data, 
             {.event.e = data.comb[[es.binary.raw.vars[1]]]
             .n.e = data.comb[[es.binary.raw.vars[3]]]
             .event.c = data.comb[[es.binary.raw.vars[2]]]
             .n.c = data.comb[[es.binary.raw.vars[4]]]}) -> mComb$value$data
      
      mCombRes = with(mComb$value,{
        data.frame(
          k = k,
          g = ifelse(
            isTRUE(common) & !isTRUE(random), 
            TE.common, TE.random) %>% 
            ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
            round(round.digits),
          g.ci = paste0("[", 
                        ifelse(isTRUE(common) & !isTRUE(random),
                               lower.common, lower.random) %>% 
                          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                          round(round.digits), "; ",
                        ifelse(isTRUE(common) & !isTRUE(random),
                               upper.common, upper.random) %>% 
                          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                          round(round.digits), "]"),
          p = ifelse(isTRUE(common) & !isTRUE(random), 
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
            ifelse(isTRUE(common) & !isTRUE(random), 
                   abs(TE.common), abs(TE.random)), nnt.cer) %>%
            ifelse(identical(.type.es, "RR"), NA, .) %>% 
            ifelse(isTRUE(.raw.bin.es), 
                   nnt.raw.bin.es, .) %>% 
            round(round.digits) %>% abs(),
          excluded = "none")
      })
    }
    
    rownames(mCombRes) = "Combined"
    mCombRes$excluded = 
      paste(
        "combined",
        ifelse(identical(study.var.comb[1], "study.var.comb"),
               "(arm-level):", "(study-level):"),
        ifelse(length(multi.study.comb) > 0,
               paste(multi.study.comb, collapse = "; "), 
               "none"))
    
    # Add rho, which.combine, combine var to mComb model
    mComb$rho = rho.within.study
    mComb$which.combine = 
      ifelse(identical(study.var.comb[1], "study.var.comb"), "arms", "studies")
    mComb$which.combine.var = which.combine.var
    
    # Add rob data if specified
    if (!is.null(rob.data[1])) {
      mComb$value = addRobData(mComb$value, rob.data)
    }
    
    # return
    return(list(m = mComb$value, 
                res = mCombRes,
                has.error = mComb$has.error,
                message = mComb$error$message))
  }
}


#' Fit 'outliers' model
#' @keywords internal 
fitOutliersModel = function(data, study.var, multi.study,
                            mGeneral, .type.es, round.digits,
                            .raw.bin.es, nnt.cer, which.run,
                            which.outliers, method.tau,
                            m.for.outliers, rob.data){
  if (method.tau == "FE"){
    mOutliers = 
      tryCatch2(
        metapsyFindOutliers(m.for.outliers)) %>% 
      {.$value = .$value$m.fixed;.}
  } else {
    mOutliers = 
      tryCatch2(
        metapsyFindOutliers(m.for.outliers)) %>% 
      {.$value = .$value$m.random;.}
  }
  if (isTRUE(.raw.bin.es) &&
      !mOutliers$has.error){
    with(updateMeta(mOutliers$value, 
                           sm = "RD", id = NULL), {
                             ifelse(isTRUE(common) & !isTRUE(random),
                                    abs(TE.common)^-1, abs(TE.random)^-1)
                           }) -> nnt.raw.bin.es
  } else {
    nnt.raw.bin.es = NA
  }
  if (mOutliers$has.error){
    mOutliersRes = with(mOutliers$value,{
      data.frame(
        k = NA, g = NA, g.ci = NA,
        p = NA, i2 = NA, i2.ci = NA,
        prediction.ci = NA, nnt = NA,
        excluded = "none")
    })
  } else {
    mOutliersRes = with(mOutliers$value,{
      data.frame(
        k = k,
        g = ifelse(
          isTRUE(common) & !isTRUE(random), 
          TE.common, TE.random) %>% 
          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
          round(round.digits),
        g.ci = paste0("[", 
                      ifelse(isTRUE(common) & !isTRUE(random),
                             lower.common, lower.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "; ",
                      ifelse(isTRUE(common) & !isTRUE(random),
                             upper.common, upper.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "]"),
        p = ifelse(isTRUE(common) & !isTRUE(random), 
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
          ifelse(isTRUE(common) & !isTRUE(random), 
                 abs(TE.common), abs(TE.random)), nnt.cer) %>%
          ifelse(identical(.type.es, "RR"), NA, .) %>% 
          ifelse(isTRUE(.raw.bin.es), 
                 nnt.raw.bin.es, .) %>% 
          round(round.digits) %>% abs(),
        excluded = "none")
    })
  }
  rownames(mOutliersRes) = "Outliers removed"
  if (length(mOutliers$value$data$comparison[mOutliers$value$exclude]) != 0 ||
      sum(mOutliers$value$exclude) > 0){
    if (identical(which.outliers[1], "combined")){
      mOutliersRes$excluded = paste(mOutliers$value$studlab
                                    [mOutliers$value$exclude], collapse = "; ")
    } else {
      mOutliersRes$excluded = paste(mOutliers$value$data$comparison
                                    [mOutliers$value$exclude], collapse = "; ")
    }
  } else {
    mOutliersRes$excluded = "no outliers detected"
  }
  # Add rob data if specified
  if (!is.null(rob.data[1])) {
    mOutliers$value = addRobData(mOutliers$value, rob.data)
  }
  return(list(m = mOutliers$value, 
              res = mOutliersRes,
              has.error = mOutliers$has.error,
              message = mOutliers$error$message))
}



#' Select which model should be used for outliers analysis
#' @keywords internal 
selectOutlierModel = function(
    which.run, mGeneral, mComb,
    which.outliers){
  
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
    m.for.outliers = mGeneral$m
  }
  
  if (which.outliers[1] == "combined"){
    if (is.null(mComb)) {
      stop("If which.outliers = 'combined', the 'combined' model must be",
           " included in which.run.", call. = FALSE)
    }
    m.for.outliers = mComb$m
  }
  
  if (!(which.outliers[1] %in% c("general", "combined", "overall"))){
    stop("'which.outliers' must either be 'overall' or 'combined'.")
  }
  return(list(m.for.outliers  = m.for.outliers,
              which.outliers = which.outliers))
}


#' Fit 'influence' model
#' @keywords internal 
fitInfluenceModel = function(which.influence, mComb,
                             mGeneral, which.run, method.tau,
                             .raw.bin.es, .type.es, round.digits,
                             nnt.cer, rob.data){
  
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
    m.for.influence = mGeneral$m
  }
  
  if (which.influence[1] == "combined"){
    if (is.null(mComb)) {
      stop("If which.influence = 'combined', the 'combined' model must be",
           " included in which.run.", call. = FALSE)
    }
    m.for.influence = mComb$m
  }
  
  if (!(which.influence[1] %in% c("general", "combined", "overall"))){
    stop("'which.influence' must either be 'overall' or 'combined'.")
  }
  
  influenceRes =
    tryCatch2(
      metapsyInfluenceAnalysis(
        m.for.influence,
        random = ifelse(method.tau == "FE",
                        FALSE, TRUE)))
  influenceRes = influenceRes$value
  
  mInfluence = 
    tryCatch2(
      updateMeta(
        m.for.influence,
        exclude = influenceRes$Data$is.infl == "yes",
        id = NULL))
  
  if (isTRUE(.raw.bin.es) &&
      !mInfluence$has.error){
    with(updateMeta(mInfluence$value, sm="RD", id = NULL), {
      ifelse(isTRUE(common) & !isTRUE(random),
             abs(TE.common)^-1, abs(TE.random)^-1)
    }) -> nnt.raw.bin.es
  } else {
    nnt.raw.bin.es = NA
  }
  
  if (mInfluence$has.error){
    mInfluenceRes = 
      data.frame(
        k = NA, g = NA, g.ci = NA,
        p = NA, i2 = NA, i2.ci = NA,
        prediction.ci = NA, nnt = NA,
        excluded = "none")
  } else {
    mInfluenceRes = with(mInfluence$value, {
      data.frame(
        k = k,
        g = ifelse(
          isTRUE(common) & !isTRUE(random), 
          TE.common, TE.random) %>% 
          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
          round(round.digits),
        g.ci = paste0("[", 
                      ifelse(isTRUE(common) & !isTRUE(random),
                             lower.common, lower.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "; ",
                      ifelse(isTRUE(common) & !isTRUE(random),
                             upper.common, upper.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "]"),
        p = ifelse(isTRUE(common) & !isTRUE(random), 
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
          ifelse(isTRUE(common) & !isTRUE(random), 
                 abs(TE.common), abs(TE.random)), nnt.cer) %>%
          ifelse(identical(.type.es, "RR"), NA, .) %>% 
          ifelse(isTRUE(.raw.bin.es), 
                 nnt.raw.bin.es, .) %>% 
          round(round.digits) %>% abs(),
        excluded = "none")
    })
  }
  rownames(mInfluenceRes) = "Influence Analysis"
  
  if (identical(which.influence[1], "combined")){
    mInfluenceRes$excluded = 
      paste("removed as influential cases:",
            ifelse(sum(influenceRes$Data$is.infl == "yes") > 0,
                   paste(mInfluence$value$studlab[influenceRes$Data$is.infl == "yes"],
                         collapse = "; "), "none"))
  } else {
    mInfluenceRes$excluded = 
      paste("removed as influential cases:",
            ifelse(sum(influenceRes$Data$is.infl == "yes") > 0,
                   paste(mInfluence$value$data$comparison
                         [influenceRes$Data$is.infl == "yes"],
                         collapse = "; "), "none"))
  }
  # Add rob data if specified
  if (!is.null(rob.data[1])) {
    mInfluence$value = addRobData(mInfluence$value, rob.data)
  }
  
  return(list(m = mInfluence$value, 
              res = mInfluenceRes,
              has.error = mInfluence$has.error,
              message = mInfluence$error$message,
              influenceRes = influenceRes))
}


#' Fit 'rob' model
#' @keywords internal 
fitRobModel = function(which.run, which.rob, which.outliers,
                       mGeneral, mComb, low.rob.filter,
                       method.tau, .raw.bin.es, .type.es, round.digits,
                       nnt.cer, rob.data){
  
  has.bs = FALSE
  
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
    m.for.rob = mGeneral$m
    data.for.rob = mGeneral$m$data
  }
  
  if (which.rob[1] == "combined"){
    if (is.null(mComb)) {
      stop("If which.influence = 'overall', the 'overall' model must be",
           " included in which.run.", call. = FALSE)
    }
    m.for.rob = mComb$m
    data.for.rob = mComb$m$data
  }
  
  if (!(which.rob[1] %in% c("general", "combined", "overall"))){
    stop("'which.rob' must either be 'general' or 'combined'.")
  }
  
  if ("rob" %in% which.run){
    robVar = strsplit(low.rob.filter, " ")[[1]][1]
    if (!robVar %in% colnames(data.for.rob)){
      stop("'", robVar, 
           "' variable not found in data to conduct analyses with low RoB only.",
           call. = FALSE)
    }
    
    m.for.rob$data[[robVar]] = 
      as.numeric(m.for.rob$data[[robVar]])
    robFilter = paste0("data.for.rob$", low.rob.filter, 
                       " & !is.na(data.for.rob$", robVar, ")")
    robMask = eval(parse(text = robFilter))
    if (sum(robMask) == 0){
      message("\n- ", crayon::yellow("[!] "), 
              "No low risk of bias studies detected! Switching to 'general'... ",
              appendLF = FALSE)
      if (mGeneral$has.error){
        robMask = NA
      } else {
        robMask = rep(TRUE, nrow(mGeneral$m$data))
      }
      mRob = tryCatch2(
        updateMeta(
          mGeneral$m, exclude = !robMask,id = NULL))
      which.run[!which.run == "rob"] -> which.run
      warn.end = TRUE
    } else {
      mRob = tryCatch2(
        updateMeta(
          m.for.rob, exclude = !robMask, id = NULL))
    }
  } else {
    robMask = rep(TRUE, nrow(mGeneral$m$data))
    mRob = tryCatch2(
      updateMeta(
        mGeneral$m, exclude = !robMask, id = NULL))
  }
  if (isTRUE(.raw.bin.es) &&
      !mRob$has.error) {
    with(updateMeta(
      mRob$value, sm="RD",
      id = NULL), {
        ifelse(isTRUE(common) & !isTRUE(random),
               abs(TE.common)^-1, abs(TE.random)^-1)
      }) -> nnt.raw.bin.es
  } else {
    nnt.raw.bin.es = NA
  }
  if (mRob$has.error){
    mRobRes = 
      data.frame(
        k = NA, g = NA, g.ci = NA,
        p = NA, i2 = NA, i2.ci = NA,
        prediction.ci = NA, nnt = NA,
        excluded = "none")
  } else {
    mRobRes = with(mRob$value,{
      data.frame(
        k = k,
        g = ifelse(
          isTRUE(common) & isTRUE(common), 
          TE.common, TE.random) %>% 
          ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
          round(round.digits),
        g.ci = paste0("[", 
                      ifelse(isTRUE(common) & !isTRUE(random),
                             lower.common, lower.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "; ",
                      ifelse(isTRUE(common) & !isTRUE(random),
                             upper.common, upper.random) %>% 
                        ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                        round(round.digits), "]"),
        p = ifelse(isTRUE(common) & !isTRUE(random), 
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
          ifelse(isTRUE(common) & !isTRUE(random), 
                 abs(TE.common), abs(TE.random)), nnt.cer) %>%
          ifelse(identical(.type.es, "RR"), NA, .) %>% 
          ifelse(isTRUE(.raw.bin.es), 
                 nnt.raw.bin.es, .) %>% 
          round(round.digits) %>% abs(),
        excluded = "none")
    })
  }
  rownames(mRobRes) = paste("Only", low.rob.filter)
  
  if (sum(!robMask) != 0){
    if (identical(which.rob[1], "combined")){
      mRobRes$excluded = 
        paste(mRob$value$studlab[!robMask], 
              collapse = "; ")
    } else {
      mRobRes$excluded = 
        paste(mRob$value$data$comparison[!robMask], 
              collapse = "; ")
    }
  } else {
    mRobRes$excluded = paste0(
      "no studies removed; ",
      low.rob.filter, " applies to all studies")
  }
  # Add rob data if specified
  if (!is.null(rob.data[1])) {
    mRob$value = addRobData(mRob$value, rob.data)
  }
  mRob$value$title2 = paste("Only", low.rob.filter)
  return(list(m = mRob$value, 
              res = mRobRes,
              has.error = mRob$has.error,
              message = mRob$error$message))
}


#' Fit 'threelevel' model
#' @keywords internal 
fitThreeLevelModel = function(data, es.var, se.var, arm.var.1, arm.var.2,
                              measure.var, study.var, .raw.bin.es, .type.es, hakn,
                              method.tau.meta, method.tau.ci, method.tau,
                              dots, es.binary.raw.vars, round.digits,
                              nnt.cer, which.run, mGeneral, mCombined,
                              use.rve, i2.ci.boot, nsim.boot){
  
  has.bs = FALSE
  
  # Define multi.study
  multi.study = names(table(data[[study.var]])
                      [table(data[[study.var]]) > 1])
  
  data$es.id = 1:nrow(data)
  formula = paste0("~ 1 | ", colnames(data[study.var]), " / es.id")
  
  mThreeLevelArgs = list(
    yi = mGeneral$m[["TE"]],
    V = mGeneral$m[["seTE"]]^2,
    slab = mGeneral$m[["studlab"]],
    data = data,
    random = as.formula(formula),
    test = ifelse(hakn == TRUE, "t", "z"),
    method = "REML") %>% 
    append(selectArguments(metafor::rma.mv, dots))
  
  mThreeLevel = 
    tryCatch2(do.call(metafor::rma.mv, mThreeLevelArgs))
  
  model.threelevel.legacy = list(
    slab = data[[study.var]],
    data = data,
    es.var = es.var[1],
    se.var = se.var[1],
    yi = mGeneral$m[["TE"]],
    V = mGeneral$m[["seTE"]]*mGeneral$m[["seTE"]],
    formula.rnd = as.formula(formula))
  mThreeLevel$value$legacy = model.threelevel.legacy
  
  if (isTRUE(.raw.bin.es)){
    # For RR analyses: re-run analyses using g
    # This is needed for NNTs
    event.data = 
      data.frame(event = data[[es.binary.raw.vars[2]]],
                 n = data[[es.binary.raw.vars[4]]])
    tryCatch2(
      meta::metaprop(
        event = event, n = n,
        data = event.data[complete.cases(event.data),])) -> m.metaprop
    
    if (m.metaprop$has.error) {
      cer = NA
    } else {
      m.metaprop$value %>% 
        {.$TE.random} %>% 
        {exp(.)/(1+exp(.))} -> cer    
    }
    
    apply(data, 1, function(x){
      esc::esc_2x2(grp1yes = as.numeric(x[[es.binary.raw.vars[1]]]) + 0.5, 
                   grp1no = as.numeric(x[[es.binary.raw.vars[3]]]) - 
                     as.numeric(x[[es.binary.raw.vars[1]]]) + 0.5,
                   grp2yes = as.numeric(x[[es.binary.raw.vars[2]]]) + 0.5,
                   grp2no = as.numeric(x[[es.binary.raw.vars[4]]]) -
                     as.numeric(x[[es.binary.raw.vars[2]]]) + 0.5,
                   es.type = "g") %>% 
        suppressWarnings() %>% 
        {data.frame(es = .$es, se = .$se)}
    }) %>% 
      do.call(rbind, .) -> data.g
    
    tryCatch2(
      metafor::rma.mv(
        yi = data.g$es, V = data.g$se^2, 
        data = data, random = as.formula(formula),
        test = ifelse(hakn == TRUE, "t", "z"),
        method = "REML")) -> mThreeLevel.g
    
    if (mThreeLevel.g$has.error){
      nnt.g = NA
    } else {
      mThreeLevel.g$value %>% 
        {.[["b"]][1]} %>% 
        {ifelse(.==0, Inf, 
                metapsyNNT(abs(.), cer))} -> nnt.g
    }
  } else {
    nnt.g = NA
  }
  
  if (mThreeLevel$has.error){
    mThreeLevel$value$I2 = NA
    mThreeLevel$value$I2.between.studies = NA
    mThreeLevel$value$I2.within.studies = NA
    mThreeLevel$value$variance.components = NA
  } else {
    # Calculate total I2
    W = diag(1/(mGeneral$m[["seTE"]]^2))
    X = model.matrix(mThreeLevel$value)
    P = W - W %*% X %*% 
      solve(t(X) %*% W %*% X) %*% 
      t(X) %*% W
    with(mThreeLevel$value, {
      100 * sum(sigma2) / 
        (sum(sigma2) + (k-p)/sum(diag(P)))
    }) -> mThreeLevel$value$I2 
    
    # Calculate I2 per level
    with(mThreeLevel$value, {
      (100 * sigma2 / 
         (sum(sigma2) + (k-p)/sum(diag(P))))
    }) -> I2.bw
    mThreeLevel$value$I2.between.studies = I2.bw[1]
    mThreeLevel$value$I2.within.studies = I2.bw[2]
    
    # Get tau and I2
    with(mThreeLevel$value, 
         {data.frame(
           tau2 = c(sigma2, sum(sigma2)),
           i2 = c(I2.between.studies, 
                  I2.within.studies,
                  I2))}) -> mThreeLevel$value$variance.components
    mThreeLevel$value$variance.components$tau2 = 
      round(mThreeLevel$value$variance.components$tau2, 4)
    mThreeLevel$value$variance.components$i2 = 
      round(mThreeLevel$value$variance.components$i2, 1)
    rownames(mThreeLevel$value$variance.components) = 
      c("Between Studies", "Within Studies", "Total")
    
    # Calculate bootstrapped CIs
    if (i2.ci.boot){
      message(
        crayon::cyan(
          crayon::bold(
            "- Parametric bootstrap for i2 confidence intervals (three-level MA model) ...")))
      
      # Run bootstrapping
      sim = metafor::simulate.rma(mThreeLevel$value, nsim=nsim.boot)
      counter = 0
      mThreeLevelArgs.bs = mThreeLevelArgs
      
      sav = lapply(sim, function(x) {
        tmp = try({
          mThreeLevelArgs.bs$yi = x
          do.call(metafor::rma.mv, mThreeLevelArgs.bs)}, silent=TRUE) 
        if (inherits(tmp, "try-error")) { 
          counter <<- counter + 1
          NA
        } else {
          counter <<- counter + 1
          if (counter %in% seq(0, nsim.boot, nsim.boot/100)){
            cat(crayon::green(
              paste0((counter/nsim.boot)*100, "% completed | ")))
          }
          if (identical(counter,nsim.boot)){
            cat(crayon::green("DONE \n"))
          }
          tmp
        }})
      
      sav = sav[!is.na(sav)]
      
      # Extract bootstrap cis (sigma)
      rbind(
        sapply(sav, function(x) x$sigma2[1]) %>% 
          quantile(c(.025, .975)),
        sapply(sav, function(x) x$sigma2[2]) %>% 
          quantile(c(.025, .975)),
        sapply(sav, function(x) x$sigma2[1] + x$sigma2[2]) %>% 
          quantile(c(.025, .975))) -> sigma2.ci
      
      # Extract bootstrap SD (sigma)
      rbind(
        sapply(sav, function(x) x$sigma2[1]) %>% sd(),
        sapply(sav, function(x) x$sigma2[2]) %>% sd(),
        sapply(sav, function(x) x$sigma2[1] + x$sigma2[2]) %>% 
          sd()) -> se.sigma2
      
      # Extract bootstrap cis (i2)
      rbind(
        sapply(sav, function(x) 100 * x$sigma2[1] / 
                 (sum(x$sigma2) + (mThreeLevel$value$k-mThreeLevel$value$p)/
                    sum(diag(P)))) %>% 
          quantile(c(.025, .975)),
        sapply(sav, function(x) 100 * x$sigma2[2] / 
                 (sum(x$sigma2) + (mThreeLevel$value$k-mThreeLevel$value$p)/
                    sum(diag(P)))) %>% 
          quantile(c(.025, .975)),
        sapply(sav, function(x) 100 * sum(x$sigma2) / 
                 (sum(x$sigma2) + (mThreeLevel$value$k-mThreeLevel$value$p)/
                    sum(diag(P)))) %>% 
          quantile(c(.025, .975))) -> i2.ci
      
      mThreeLevel$value$variance.components$tau2.ci = 
        c(paste0("[", round(sigma2.ci[1,], round.digits+1) %>% 
                   paste(collapse = "; "), "]"),
          paste0("[", round(sigma2.ci[2,], round.digits+1) %>% 
                   paste(collapse = "; "), "]"), 
          paste0("[", round(sigma2.ci[3,], round.digits+1) %>% 
                   paste(collapse = "; "), "]"))
      
      mThreeLevel$value$variance.components$i2.ci = 
        c(paste0("[", round(i2.ci[1,], round.digits) %>% 
                   paste(collapse = "; "), "]"),
          paste0("[", round(i2.ci[2,], round.digits) %>% 
                   paste(collapse = "; "), "]"), 
          paste0("[", round(i2.ci[3,], round.digits) %>% 
                   paste(collapse = "; "), "]"))
      
      mThreeLevel$value$variance.components =
        mThreeLevel$value$variance.components[,c("tau2", "tau2.ci", "i2", "i2.ci")] 
      
      mThreeLevel$value$se.sigma2 = se.sigma2
      
      has.bs = TRUE
    } 
  }
  
  if (use.rve[1] == FALSE){
    
    if ("threelevel" %in% which.run){
      message("- ", crayon::green("[OK] "), 
              "Calculating effect size using three-level MA model... ",
              appendLF = FALSE)
    }
    
    if (mThreeLevel$has.error){
      mThreeLevelRes = 
        data.frame(k = NA, g = NA, g.ci = NA,
                   p = NA, i2 = NA, i2.ci = NA,
                   prediction.ci = NA, nnt = NA,
                   excluded = "none")
    } else {
      mThreeLevelRes = with(mThreeLevel$value, {
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
                   i2.ci = ifelse(has.bs, variance.components$i2.ci[3], "-"),
                   prediction.ci = paste0("[", 
                                          round(predict(mThreeLevel$value)$pi.lb %>% 
                                                  ifelse(identical(.type.es, "RR"), exp(.), .), 
                                                round.digits), "; ",
                                          round(predict(mThreeLevel$value)$pi.ub %>% 
                                                  ifelse(identical(.type.es, "RR"), exp(.), .), 
                                                round.digits), "]"),
                   nnt = ifelse(identical(.type.es, "RR"), 
                                NA, metapsyNNT(abs(as.numeric(b[,1])), nnt.cer)) %>% 
                     ifelse(isTRUE(.raw.bin.es), nnt.g, .) %>% 
                     round(round.digits) %>% abs(),
                   excluded = "none")
      })
    }
    rownames(mThreeLevelRes) = "Three-Level Model"
    
    if (!mThreeLevel$has.error){
      mThreeLevelRes$excluded = 
        paste("Number of clusters/studies:", 
              mThreeLevel$value$s.nlevels[1]) 
    }
    
    sendMessage(mThreeLevel, "threelevel")
    
    if (length(multi.study) == 0 & 
        "threelevel" %in% which.run){
      message("\n- ", crayon::yellow("[!] "), 
              "All included ES seem to be independent.",
              " A three-level model is not adequate and ",
              " tau/I2 estimates are not trustworthy! ",
              appendLF = FALSE)
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
    
    if (mThreeLevel$has.error){
      mThreeLevelRes.RVE = 
        with(mThreeLevel$value, {
          data.frame(k = NA, g = NA, g.ci = NA,
                     p = NA, i2 = NA, i2.ci = NA,
                     prediction.ci = NA, nnt = NA,
                     excluded = "none")})
    } else {
      # Get results: RVE used
      crit = qt(0.025, mThreeLevel$value$ddf[[1]], lower.tail = FALSE)
      tau2 = sum(mThreeLevel$value$sigma2)
      SE = clubSandwich::conf_int(mThreeLevel$value, vcov = "CR2")[["SE"]]
      pi.lb.rve = 
        clubSandwich::conf_int(mThreeLevel$value, "CR2")[["beta"]] -
        crit * sqrt(tau2 + (SE^2))
      pi.ub.rve = 
        clubSandwich::conf_int(mThreeLevel$value, "CR2")[["beta"]] +
        crit * sqrt(tau2 + (SE^2))
      if (isTRUE(.raw.bin.es)){
        as.numeric(clubSandwich::conf_int(
          mThreeLevel.g$value, "CR2")[["beta"]]) %>% 
          {ifelse(.==0, Inf, metapsyNNT(abs(.), cer))} -> nnt.g
      }
      mThreeLevelRes.RVE = with(mThreeLevel$value, {
        data.frame(k = k.all,
                   g = as.numeric(
                     clubSandwich::conf_int(
                       mThreeLevel$value, "CR2")[["beta"]]) %>%
                     ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                     round(round.digits),
                   g.ci = paste0(
                     "[",
                     clubSandwich::conf_int(
                       mThreeLevel$value, "CR2")[["CI_L"]] %>%
                       ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                       round(round.digits), "; ",
                     clubSandwich::conf_int(
                       mThreeLevel$value, "CR2")[["CI_U"]] %>%
                       ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                       round(round.digits), "]"),
                   p = clubSandwich::coef_test(
                     mThreeLevel$value, "CR2")[["p_Satt"]] %>%
                     scales::pvalue(),
                   i2 = round(I2, 1),
                   i2.ci = ifelse(has.bs, variance.components$i2.ci[3], "-"),
                   prediction.ci = paste0(
                     "[", 
                     round(pi.lb.rve %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "; ",
                     round(pi.ub.rve %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "]"),
                   nnt = ifelse(identical(.type.es, "RR"), 
                                NA, metapsyNNT(abs(as.numeric(
                                  clubSandwich::conf_int(
                                    mThreeLevel$value, "CR2")[["beta"]])), 
                                  nnt.cer)) %>% 
                     ifelse(isTRUE(.raw.bin.es), nnt.g, .) %>% 
                     round(round.digits) %>% abs(),
                   excluded = "none")
      })
      mThreeLevelRes.RVE$excluded = 
        paste0("Number of clusters/studies: ", 
               mThreeLevel$value$s.nlevels[1], 
               "; robust variance estimation (RVE) used.")
    }
    rownames(mThreeLevelRes.RVE) = "Three-Level Model"
    mThreeLevelRes = mThreeLevelRes.RVE
    
    sendMessage(mThreeLevel, "threelevel")
    
    if (length(multi.study) == 0 & "threelevel" %in% which.run){
      message("\n- ", crayon::yellow("[!] "), 
              "All included ES seem to be independent.",
              " A three-level model is not adequate and ",
              " tau/I2 estimates are not trustworthy! ",
              appendLF = FALSE)
      warn.end = TRUE
      which.run[!which.run == "threelevel"] -> which.run
    }
  }
  
  if (has.bs){
    mThreeLevel$value$i2.ci.boot = TRUE
  } else {
    mThreeLevel$value$i2.ci.boot = FALSE
  }
  
  return(list(m = mThreeLevel$value, 
              res = mThreeLevelRes,
              has.error = mThreeLevel$has.error,
              message = mThreeLevel$error$message))
}


#' Fit 'threelevel.che' model
#' @keywords internal 
fitThreeLevelCHEModel = function(data, es.var, se.var, arm.var.1, arm.var.2,
                                 measure.var, study.var, .raw.bin.es, .type.es, hakn,
                                 method.tau.meta, method.tau.ci, method.tau,
                                 dots, es.binary.raw.vars, round.digits,
                                 nnt.cer, which.run, mGeneral, mCombined,
                                 use.rve, rho.within.study, phi.within.study,
                                 i2.ci.boot, nsim.boot){
  
  has.bs = FALSE
  
  if (rho.within.study[1] >.99){
    message("- ", crayon::yellow("[!] "), 
            "'rho.within.study' is very close to 1.",
            " This can lead to non-positive-definite variance-covariance matrices.",
            " If the V matrix is not positive-definite, assume a lower value...")
    warn.end = TRUE}
  if (phi.within.study[1] >.91){
    message("- ", crayon::yellow("[!] "), 
            "'phi.within.study' is very close to 1.",
            " This can lead to non-positive-definite variance-covariance matrices.",
            " If the V matrix is not positive-definite, assume a lower value...")
    warn.end = TRUE}
  if (phi.within.study[1] <.1){
    message("- ", crayon::yellow("[!] "), 
            "'phi.within.study' is very close to 0.",
            " This can lead to non-positive-definite variance-covariance matrices.",
            " If the V matrix is not positive-definite, assume a higher value...")
    warn.end = TRUE
  }
  if (rho.within.study[1] <.01){
    message("- ", crayon::yellow("[!] "), 
            "'rho.within.study' is very close to 0.",
            " This can lead to non-positive-definite variance-covariance matrices.",
            " If the V matrix is not positive-definite, assume a higher value...")
    warn.end = TRUE
  }
  
  # Define multi.study
  multi.study = names(table(data[[study.var]])
                      [table(data[[study.var]]) > 1])
  
  data$es.id = 1:nrow(data)
  formula.fixed = as.formula("TE ~ 1")
  formula.rnd = as.formula(paste0("~ 1 | ", 
                                  colnames(data[study.var]), "/ es.id"))
  data.che = data
  data.che$TE = mGeneral$m[["TE"]]
  
  Vmat = clubSandwich::impute_covariance_matrix(
    mGeneral$m[["seTE"]]^2,
    cluster = data[[study.var]],
    r = rho.within.study,
    smooth_vi = TRUE)
  
  mCHEArgs = list(
    formula.fixed,
    V = Vmat,
    slab = data[[study.var]],
    data = data.che,
    random = formula.rnd,
    test = ifelse(hakn == TRUE, "t", "z"),
    method = "REML", 
    sparse = TRUE) %>% 
    append(selectArguments(metafor::rma.mv, dots))
  
  mCHE = 
    tryCatch2(
      do.call(metafor::rma.mv, mCHEArgs))
  
  model.threelevel.che.legacy = list(
    slab = data.che[[study.var]],
    data = data.che,
    formula.rnd = formula.rnd,
    formula.fixed = formula.fixed,
    Vmat = Vmat)
  mCHE$value$legacy = 
    model.threelevel.che.legacy
  
  if (isTRUE(.raw.bin.es)){
    # For RR analyses: re-run analyses using g
    # This is needed for NNTs
    event.data = 
      data.frame(event = data.che[[es.binary.raw.vars[2]]],
                 n = data.che[[es.binary.raw.vars[4]]])
    tryCatch2(
      meta::metaprop(
        event = event, n = n,
        data = event.data[complete.cases(event.data),])) -> m.metaprop
    
    if (m.metaprop$has.error) {
      cer = NA
    } else {
      m.metaprop$value %>% 
        {.$TE.random} %>% 
        {exp(.)/(1+exp(.))} -> cer    
    }
    
    apply(data.che, 1, function(x){
      esc::esc_2x2(grp1yes = as.numeric(x[[es.binary.raw.vars[1]]]) + 0.5, 
                   grp1no = as.numeric(x[[es.binary.raw.vars[3]]]) - 
                     as.numeric(x[[es.binary.raw.vars[1]]]) + 0.5,
                   grp2yes = as.numeric(x[[es.binary.raw.vars[2]]]) + 0.5,
                   grp2no = as.numeric(x[[es.binary.raw.vars[4]]]) -
                     as.numeric(x[[es.binary.raw.vars[2]]]) + 0.5,
                   es.type = "g") %>% 
        suppressWarnings() %>% 
        {data.frame(es = .$es, se = .$se)}
    }) %>% 
      do.call(rbind, .) -> data.g
    
    Vmat.g = clubSandwich::impute_covariance_matrix(
      data.g$se^2,
      cluster = data.che[[study.var]],
      r = rho.within.study,
      smooth_vi = TRUE)
    
    data.che$TE = data.g$es
    tryCatch2(
      metafor::rma.mv(
        formula.fixed,
        V = Vmat.g,
        slab = data.che[[study.var]],
        data = data.che,
        random = formula.rnd,
        test = ifelse(hakn == TRUE, "t", "z"),
        method = "REML", 
        sparse = TRUE)) -> mCHE.g
    
    if (mCHE.g$has.error){
      nnt.g = NA
    } else {
      mCHE.g$value %>% 
        {.[["b"]][1]} %>% 
        {ifelse(.==0, Inf, 
                metapsyNNT(abs(.), cer))} -> nnt.g
    }
  } else {
    nnt.g = NA
  }
  
  if (mCHE$has.error){
    mCHE$value$I2 = NA
    mCHE$value$I2.between.studies = NA
    mCHE$value$I2.within.studies = NA
    mCHE$value$variance.components = NA
  } else {
    # Calculate total I2
    W = diag(1/(mGeneral$m[["seTE"]]^2))
    X = model.matrix(mCHE$value)
    P = W - W %*% X %*% 
      solve(t(X) %*% W %*% X) %*% 
      t(X) %*% W
    with(mCHE$value, {
      100 * sum(sigma2) / 
        (sum(sigma2) + (k-p)/sum(diag(P)))
    }) -> mCHE$value$I2 
    
    # Calculate I2 per level
    with(mCHE$value, {
      (100 * sigma2 / 
         (sum(sigma2) + (k-p)/sum(diag(P))))
    }) -> I2.bw
    
    mCHE$value$I2.between.studies = I2.bw[1]
    mCHE$value$I2.within.studies = I2.bw[2]
    
    # Get tau and I2
    with(mCHE$value, 
         {data.frame(
           tau2 = c(sigma2, sum(sigma2)),
           i2 = c(I2.between.studies, 
                  I2.within.studies,
                  I2))}) -> mCHE$value$variance.components
    
    mCHE$value$variance.components$tau2 = 
      round(mCHE$value$variance.components$tau2, 4)
    mCHE$value$variance.components$i2 = 
      round(mCHE$value$variance.components$i2, 1)
    rownames(mCHE$value$variance.components) = 
      c("Between Studies", "Within Studies", "Total")
    
    # Calculate bootstrapped CIs
    if (i2.ci.boot){
      message(
        crayon::cyan(
          crayon::bold(
            "- Parametric bootstrap for i2 confidence intervals (three-level CHE model) ...")))
      
      # Run bootstrapping
      sim = metafor::simulate.rma(mCHE$value, nsim=nsim.boot)
      counter = 0
      mCHEArgs.bs = mCHEArgs
      
      sav = lapply(sim, function(x) {
        tmp = try({
          mCHEArgs.bs$data$TE = x
          do.call(metafor::rma.mv, mCHEArgs.bs)}, silent=TRUE) 
        if (inherits(tmp, "try-error")) { 
          counter <<- counter + 1
          NA
        } else {
          counter <<- counter + 1
          if (counter %in% seq(0, nsim.boot, nsim.boot/100)){
            cat(crayon::green(
              paste0((counter/nsim.boot)*100, "% completed | ")))
          }
          if (identical(counter,nsim.boot)){
            cat(crayon::green("DONE \n"))
          }
          tmp
        }})
      
      sav = sav[!is.na(sav)]
      
      # Extract bootstrap cis (sigma)
      rbind(
        sapply(sav, function(x) x$sigma2[1]) %>% 
          quantile(c(.025, .975)),
        sapply(sav, function(x) x$sigma2[2]) %>% 
          quantile(c(.025, .975)),
        sapply(sav, function(x) x$sigma2[1] + x$sigma2[2]) %>% 
          quantile(c(.025, .975))) -> sigma2.ci
      
      # Extract bootstrap SD (sigma)
      rbind(
        sapply(sav, function(x) x$sigma2[1]) %>% sd(),
        sapply(sav, function(x) x$sigma2[2]) %>% sd(),
        sapply(sav, function(x) x$sigma2[1] + x$sigma2[2]) %>% 
          sd()) -> se.sigma2
      
      # Extract bootstrap cis (i2)
      rbind(
        sapply(sav, function(x) 100 * x$sigma2[1] / 
                 (sum(x$sigma2) + (mCHE$value$k-mCHE$value$p)/
                    sum(diag(P)))) %>% 
          quantile(c(.025, .975)),
        sapply(sav, function(x) 100 * x$sigma2[2] / 
                 (sum(x$sigma2) + (mCHE$value$k-mCHE$value$p)/
                    sum(diag(P)))) %>% 
          quantile(c(.025, .975)),
        sapply(sav, function(x) 100 * sum(x$sigma2) / 
                 (sum(x$sigma2) + (mCHE$value$k-mCHE$value$p)/
                    sum(diag(P)))) %>% 
          quantile(c(.025, .975))) -> i2.ci
      
      mCHE$value$variance.components$tau2.ci = 
        c(paste0("[", round(sigma2.ci[1,], round.digits+1) %>% 
                   paste(collapse = "; "), "]"),
          paste0("[", round(sigma2.ci[2,], round.digits+1) %>% 
                   paste(collapse = "; "), "]"), 
          paste0("[", round(sigma2.ci[3,], round.digits+1) %>% 
                   paste(collapse = "; "), "]"))
      
      mCHE$value$variance.components$i2.ci = 
        c(paste0("[", round(i2.ci[1,], round.digits) %>% 
                   paste(collapse = "; "), "]"),
          paste0("[", round(i2.ci[2,], round.digits) %>% 
                   paste(collapse = "; "), "]"), 
          paste0("[", round(i2.ci[3,], round.digits) %>% 
                   paste(collapse = "; "), "]"))
      
      mCHE$value$variance.components =
        mCHE$value$variance.components[,c("tau2", "tau2.ci", "i2", "i2.ci")] 
      
      mCHE$value$se.sigma2 = se.sigma2
      
      has.bs = TRUE
    } 
    else {
      has.bs = FALSE
    }
  }
  
  
  # Get results: no RVE
  if (use.rve[1] == FALSE){
    if ("threelevel.che" %in% which.run){
      message("- ", crayon::green("[OK] "), 
              "Calculating effect size using three-level CHE model... ",
              appendLF = FALSE)
    }
    
    if (mCHE$has.error){
      mCHERes = 
        data.frame(k = NA, g = NA, g.ci = NA,
                   p = NA, i2 = NA, i2.ci = NA,
                   prediction.ci = NA, nnt = NA,
                   excluded = "none")
    } else {
      mCHERes = with(mCHE$value, {
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
                   i2.ci = ifelse(has.bs, variance.components$i2.ci[3], "-"),
                   prediction.ci = paste0(
                     "[", 
                     round(predict(mCHE$value)$pi.lb %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "; ",
                     round(predict(mCHE$value)$pi.ub %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "]"),
                   nnt = ifelse(identical(.type.es, "RR"), 
                                NA, metapsyNNT(abs(as.numeric(b[,1])), nnt.cer)) %>% 
                     ifelse(isTRUE(.raw.bin.es), nnt.g, .) %>% 
                     round(round.digits) %>% abs(),
                   excluded = "none")
      })
      mCHERes$excluded = paste("Number of clusters/studies:", 
                               mCHE$value$s.nlevels[1])
    }
    rownames(mCHERes) = "Three-Level Model (CHE)"
    sendMessage(mCHE, "threelevel.che")
    
    if (length(multi.study) == 0 & "threelevel.che" %in% which.run){
      message("\n- ", crayon::yellow("[!] "), 
              "All included ES seem to be independent.",
              " A three-level CHE model is not adequate and ",
              "tau/I2 estimates are not trustworthy! ",
              appendLF = FALSE)
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
    
    if (mCHE$has.error){
      mCHERes.RVE = 
        data.frame(k = NA, g = NA, g.ci = NA,
                   p = NA, i2 = NA, i2.ci = NA,
                   prediction.ci = NA, nnt = NA,
                   excluded = "none")
    } else {
      # Get results: RVE used
      crit = qt(0.025, mCHE$value$ddf[[1]], lower.tail = FALSE)
      tau2 = sum(mCHE$value$sigma2)
      SE = clubSandwich::conf_int(mCHE$value, vcov = "CR2")[["SE"]]
      pi.lb.rve = clubSandwich::conf_int(mCHE$value, "CR2")[["beta"]] -
        crit * sqrt(tau2 + (SE^2))
      pi.ub.rve = clubSandwich::conf_int(mCHE$value, "CR2")[["beta"]] +
        crit * sqrt(tau2 + (SE^2))
      if (isTRUE(.raw.bin.es)){
        as.numeric(clubSandwich::conf_int(
          mCHE.g$value, "CR2")[["beta"]]) %>% 
          {ifelse(.==0, Inf, metapsyNNT(abs(.), cer))} -> nnt.g
      }
      mCHERes.RVE = with(mCHE$value, {
        data.frame(k = k.all,
                   g = as.numeric(
                     clubSandwich::conf_int(
                       mCHE$value, "CR2")[["beta"]]) %>%
                     ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                     round(round.digits),
                   g.ci = paste0(
                     "[",
                     clubSandwich::conf_int(
                       mCHE$value, "CR2")[["CI_L"]] %>%
                       ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                       round(round.digits), "; ",
                     clubSandwich::conf_int(
                       mCHE$value, "CR2")[["CI_U"]] %>%
                       ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                       round(round.digits), "]"),
                   p = clubSandwich::coef_test(
                     mCHE$value, "CR2")[["p_Satt"]] %>%
                     scales::pvalue(),
                   i2 = round(I2, 1),
                   i2.ci = ifelse(has.bs, variance.components$i2.ci[3], "-"),
                   prediction.ci = paste0(
                     "[", 
                     round(pi.lb.rve %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "; ",
                     round(pi.ub.rve %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "]"),
                   nnt = ifelse(identical(.type.es, "RR"), 
                                NA, metapsyNNT(abs(as.numeric(
                                  clubSandwich::conf_int(
                                    mCHE$value, "CR2")[["beta"]])), nnt.cer)) %>% 
                     ifelse(isTRUE(.raw.bin.es), nnt.g, .) %>% 
                     round(round.digits) %>% abs(),
                   excluded = "none")
      })
      mCHERes.RVE$excluded = 
        paste0("Number of clusters/studies: ", mCHE$value$s.nlevels[1],
               "; robust variance estimation (RVE) used.")
    }
    
    rownames(mCHERes.RVE) = "Three-Level Model (CHE)"
    mCHERes = mCHERes.RVE
    sendMessage(mCHE, "threelevel.che")
    
    if (length(multi.study) == 0 & "threelevel.che" %in% which.run){
      message("\n- ", crayon::yellow("[!] "), 
              "All included ES seem to be independent.",
              " A three-level CHE model is not adequate and ",
              "tau/I2 estimates are not trustworthy! ",
              appendLF = FALSE)
      warn.end = TRUE
      which.run[!which.run == "threelevel.che"] -> which.run
    }
  }
  if (has.bs){
    mCHE$value$i2.ci.boot = TRUE
  } else {
    mCHE$value$i2.ci.boot = FALSE
  }
  
  return(list(m = mCHE$value, 
              res = mCHERes,
              has.error = mCHE$has.error,
              message = mCHE$error$message))
}



#' Fit 'threelevel' complex model
#' @keywords internal 
fitThreeLevelHACEModel = function(data, es.var, se.var, arm.var.1, arm.var.2,
                                  measure.var, study.var, .raw.bin.es, .type.es, hakn,
                                  method.tau.meta, method.tau.ci, method.tau,
                                  dots, es.binary.raw.vars, round.digits,
                                  nnt.cer, which.run, mGeneral, mCombined,
                                  use.rve, rho.within.study, which.combine.var,
                                  phi.within.study, n.var.arm1, 
                                  n.var.arm2, w1.var, w2.var, time.var,
                                  near.pd, i2.ci.boot, nsim.boot){
  
  has.bs = FALSE
  
  if (rho.within.study[1] >.99){
    message("- ", crayon::yellow("[!] "), 
            "'rho.within.study' is very close to 1.",
            " This can lead to non-positive-definite variance-covariance matrices.",
            " If the V matrix is not positive-definite, assume a lower value...")
    warn.end = TRUE}
  if (phi.within.study[1] >.91){
    message("- ", crayon::yellow("[!] "), 
            "'phi.within.study' is very close to 1.",
            " This can lead to non-positive-definite variance-covariance matrices.",
            " If the V matrix is not positive-definite, assume a lower value...")
    warn.end = TRUE}
  if (phi.within.study[1] <.1){
    message("- ", crayon::yellow("[!] "), 
            "'phi.within.study' is very close to 0.",
            " This can lead to non-positive-definite variance-covariance matrices.",
            " If the V matrix is not positive-definite, assume a higher value...")
    warn.end = TRUE
  }
  if (rho.within.study[1] <.01){
    message("- ", crayon::yellow("[!] "), 
            "'rho.within.study' is very close to 0.",
            " This can lead to non-positive-definite variance-covariance matrices.",
            " If the V matrix is not positive-definite, assume a higher value...")
    warn.end = TRUE
  }
  
  # Define multi.study
  multi.study = names(table(data[[study.var]])
                      [table(data[[study.var]]) > 1])
  
  data$es.id = 1:nrow(data)
  formula.fixed = as.formula("TE ~ 1")
  formula.rnd = as.formula(paste0("~ 1 | ", 
                                  colnames(data[study.var]), "/ es.id"))
  data.che = data
  data.che$TE = mGeneral$m[["TE"]]
  data.che$seTE = mGeneral$m[["seTE"]]
  
  # Define required variables
  data.che$study = data.che[[study.var]]
  data.che$instrument = data.che[[measure.var]]
  data.che$condition_arm1 = data.che[[arm.var.1]]
  data.che$condition_arm2 = data.che[[arm.var.2]]
  data.che$multi_arm1 = data.che[[which.combine.var]]
  data.che$n_arm1 = data.che[[w1.var]]
  data.che$n_arm2 = data.che[[w2.var]]
  data.che$time_weeks = data.che[[time.var]]
  within(data.che, {
    instrument[is.na(instrument)] = "missing"
    vcalc_arm1 = paste(
      study, condition_arm1, multi_arm1, "arm1")
    vcalc_arm2 = paste(
      study, condition_arm2, "arm2") 
    time_weeks[is.na(time_weeks)] = median(time_weeks, na.rm=T)
    n_arm1[is.na(n_arm1)] = 20
    n_arm2[is.na(n_arm2)] = 20
    es = TE; V = seTE^2
  }) -> dat.che
  
  # Approximate Vcovs
  tryCatch2(metafor::vcalc(vi = V, cluster = study, 
                           obs = instrument, time1 = time_weeks,
                           grp1 = vcalc_arm1, grp2 = vcalc_arm2,
                           w1 = n_arm1, w2 = n_arm2, 
                           rho = rho.within.study, 
                           phi = phi.within.study, data = dat.che, 
                           nearpd = near.pd)) -> vcalc.res
  if (!is.null(vcalc.res$warning[1])){
    append.error = TRUE
  } else {
    append.error = FALSE
  }
  metafor::blsplit(vcalc.res$value, dat.che$study) -> Vcov
  
  mCHEArgs = list(
    formula.fixed,
    V = Vcov,
    slab = data[[study.var]],
    data = data.che,
    random = formula.rnd,
    test = ifelse(hakn == TRUE, "t", "z"),
    method = "REML", 
    sparse = TRUE) %>% 
    append(selectArguments(metafor::rma.mv, dots))
  
  mCHE = 
    tryCatch2(
      do.call(metafor::rma.mv, mCHEArgs))
  
  model.threelevel.che.legacy = list(
    slab = data.che[[study.var]],
    data = data.che,
    formula.rnd = formula.rnd,
    formula.fixed = formula.fixed,
    Vmat = Vcov)
  mCHE$value$legacy = 
    model.threelevel.che.legacy
  
  if (isTRUE(.raw.bin.es)){
    # For RR analyses: re-run analyses using g
    # This is needed for NNTs
    event.data = 
      data.frame(event = data.che[[es.binary.raw.vars[2]]],
                 n = data.che[[es.binary.raw.vars[4]]])
    tryCatch2(
      meta::metaprop(
        event = event, n = n,
        data = event.data[complete.cases(event.data),])) -> m.metaprop
    
    if (m.metaprop$has.error) {
      cer = NA
    } else {
      m.metaprop$value %>% 
        {.$TE.random} %>% 
        {exp(.)/(1+exp(.))} -> cer    
    }
    
    apply(data.che, 1, function(x){
      esc::esc_2x2(grp1yes = as.numeric(x[[es.binary.raw.vars[1]]]) + 0.5, 
                   grp1no = as.numeric(x[[es.binary.raw.vars[3]]]) - 
                     as.numeric(x[[es.binary.raw.vars[1]]]) + 0.5,
                   grp2yes = as.numeric(x[[es.binary.raw.vars[2]]]) + 0.5,
                   grp2no = as.numeric(x[[es.binary.raw.vars[4]]]) -
                     as.numeric(x[[es.binary.raw.vars[2]]]) + 0.5,
                   es.type = "g") %>% 
        suppressWarnings() %>% 
        {data.frame(es = .$es, se = .$se)}
    }) %>% 
      do.call(rbind, .) -> data.g
    
    # Approximate Vcovs
    tryCatch2(
      metafor::vcalc(vi = data.g$se^2, cluster = study, 
                     obs = instrument, time1 = time_weeks,
                     grp1 = vcalc_arm1, grp2 = vcalc_arm2,
                     w1 = n_arm1, w2 = n_arm2, 
                     rho = rho.within.study, 
                     phi = phi.within.study, data = dat.che, 
                     nearpd = near.pd)) -> vcalc.res
    if (!is.null(vcalc.res$warning[1])){
      append.error = TRUE
    } else {
      append.error = FALSE
    }
    metafor::blsplit(vcalc.res$value, dat.che$study) -> Vmat.g
    
    data.che$TE = data.g$es
    tryCatch2(
      metafor::rma.mv(
        formula.fixed,
        V = Vmat.g,
        slab = data.che[[study.var]],
        data = data.che,
        random = formula.rnd,
        test = ifelse(hakn == TRUE, "t", "z"),
        method = "REML", 
        sparse = TRUE)) -> mCHE.g
    
    if (mCHE.g$has.error){
      nnt.g = NA
    } else {
      mCHE.g$value %>% 
        {.[["b"]][1]} %>% 
        {ifelse(.==0, Inf, 
                metapsyNNT(abs(.), cer))} -> nnt.g
    }
  } else {
    nnt.g = NA
  }
  
  if (mCHE$has.error){
    mCHE$value$I2 = NA
    mCHE$value$I2.between.studies = NA
    mCHE$value$I2.within.studies = NA
    mCHE$value$variance.components = NA
  } else {
    # Calculate total I2
    W = diag(1/(mGeneral$m[["seTE"]]^2))
    X = model.matrix(mCHE$value)
    P = W - W %*% X %*% 
      solve(t(X) %*% W %*% X) %*% 
      t(X) %*% W
    with(mCHE$value, {
      100 * sum(sigma2) / 
        (sum(sigma2) + (k-p)/sum(diag(P)))
    }) -> mCHE$value$I2 
    
    # Calculate I2 per level
    with(mCHE$value, {
      (100 * sigma2 / 
         (sum(sigma2) + (k-p)/sum(diag(P))))
    }) -> I2.bw
    
    mCHE$value$I2.between.studies = I2.bw[1]
    mCHE$value$I2.within.studies = I2.bw[2]
    
    # Get tau and I2
    with(mCHE$value, 
         {data.frame(
           tau2 = c(sigma2, sum(sigma2)),
           i2 = c(I2.between.studies, 
                  I2.within.studies,
                  I2))}) -> mCHE$value$variance.components
    
    mCHE$value$variance.components$tau2 = 
      round(mCHE$value$variance.components$tau2, 4)
    mCHE$value$variance.components$i2 = 
      round(mCHE$value$variance.components$i2, 1)
    rownames(mCHE$value$variance.components) = 
      c("Between Studies", "Within Studies", "Total")
    
    # Calculate bootstrapped CIs
    if (i2.ci.boot){
      message(
        crayon::cyan(
          crayon::bold(
            "- Parametric bootstrap for i2 confidence intervals (three-level CHE model) ...")))
      
      # Run bootstrapping
      sim = metafor::simulate.rma(mCHE$value, nsim=nsim.boot)
      counter = 0
      mCHEArgs.bs = mCHEArgs
      
      sav = lapply(sim, function(x) {
        tmp = try({
          mCHEArgs.bs$data$TE = x
          do.call(metafor::rma.mv, mCHEArgs.bs)}, silent=TRUE) 
        if (inherits(tmp, "try-error")) { 
          counter <<- counter + 1
          NA
        } else {
          counter <<- counter + 1
          if (counter %in% seq(0, nsim.boot, nsim.boot/100)){
            cat(crayon::green(
              paste0((counter/nsim.boot)*100, "% completed | ")))
          }
          if (identical(counter,nsim.boot)){
            cat(crayon::green("DONE \n"))
          }
          tmp
        }})
      
      sav = sav[!is.na(sav)]
      
      # Extract bootstrap cis (sigma)
      rbind(
        sapply(sav, function(x) x$sigma2[1]) %>% 
          quantile(c(.025, .975)),
        sapply(sav, function(x) x$sigma2[2]) %>% 
          quantile(c(.025, .975)),
        sapply(sav, function(x) x$sigma2[1] + x$sigma2[2]) %>% 
          quantile(c(.025, .975))) -> sigma2.ci
      
      # Extract bootstrap SD (sigma)
      rbind(
        sapply(sav, function(x) x$sigma2[1]) %>% sd(),
        sapply(sav, function(x) x$sigma2[2]) %>% sd(),
        sapply(sav, function(x) x$sigma2[1] + x$sigma2[2]) %>% 
          sd()) -> se.sigma2
      
      # Extract bootstrap cis (i2)
      rbind(
        sapply(sav, function(x) 100 * x$sigma2[1] / 
                 (sum(x$sigma2) + (mCHE$value$k-mCHE$value$p)/
                    sum(diag(P)))) %>% 
          quantile(c(.025, .975)),
        sapply(sav, function(x) 100 * x$sigma2[2] / 
                 (sum(x$sigma2) + (mCHE$value$k-mCHE$value$p)/
                    sum(diag(P)))) %>% 
          quantile(c(.025, .975)),
        sapply(sav, function(x) 100 * sum(x$sigma2) / 
                 (sum(x$sigma2) + (mCHE$value$k-mCHE$value$p)/
                    sum(diag(P)))) %>% 
          quantile(c(.025, .975))) -> i2.ci
      
      mCHE$value$variance.components$tau2.ci = 
        c(paste0("[", round(sigma2.ci[1,], round.digits+1) %>% 
                   paste(collapse = "; "), "]"),
          paste0("[", round(sigma2.ci[2,], round.digits+1) %>% 
                   paste(collapse = "; "), "]"), 
          paste0("[", round(sigma2.ci[3,], round.digits+1) %>% 
                   paste(collapse = "; "), "]"))
      
      mCHE$value$variance.components$i2.ci = 
        c(paste0("[", round(i2.ci[1,], round.digits) %>% 
                   paste(collapse = "; "), "]"),
          paste0("[", round(i2.ci[2,], round.digits) %>% 
                   paste(collapse = "; "), "]"), 
          paste0("[", round(i2.ci[3,], round.digits) %>% 
                   paste(collapse = "; "), "]"))
      
      mCHE$value$variance.components =
        mCHE$value$variance.components[,c("tau2", "tau2.ci", "i2", "i2.ci")]
      
      mCHE$value$se.sigma2 = se.sigma2
      
      has.bs = TRUE
    } 
    else {
      has.bs = FALSE
    }
  }
  
  
  # Get results: no RVE
  if (use.rve[1] == FALSE){
    if ("threelevel.che" %in% which.run){
      message("- ", crayon::green("[OK] "), 
              "Calculating effect size using 3L-CHE model (rho=", 
              rho.within.study, "; phi=", phi.within.study,
              "; multiarm corr.)... ",
              appendLF = FALSE)
    }
    
    if (mCHE$has.error){
      mCHERes = 
        data.frame(k = NA, g = NA, g.ci = NA,
                   p = NA, i2 = NA, i2.ci = NA,
                   prediction.ci = NA, nnt = NA,
                   excluded = "none")
    } else {
      mCHERes = with(mCHE$value, {
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
                   i2.ci = ifelse(has.bs, variance.components$i2.ci[3], "-"),
                   prediction.ci = paste0(
                     "[", 
                     round(predict(mCHE$value)$pi.lb %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "; ",
                     round(predict(mCHE$value)$pi.ub %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "]"),
                   nnt = ifelse(identical(.type.es, "RR"), 
                                NA, metapsyNNT(abs(as.numeric(b[,1])), nnt.cer)) %>% 
                     ifelse(isTRUE(.raw.bin.es), nnt.g, .) %>% 
                     round(round.digits) %>% abs(),
                   excluded = "none")
      })
      mCHERes$excluded = paste("Number of clusters/studies:", 
                               mCHE$value$s.nlevels[1])
    }
    rownames(mCHERes) = "Three-Level Model (CHE)"
    sendMessage(mCHE, "threelevel.che")
    
    if (length(multi.study) == 0 & "threelevel.che" %in% which.run){
      message("\n- ", crayon::yellow("[!] "), 
              "All included ES seem to be independent.",
              " A three-level CHE model is not adequate and ",
              "tau/I2 estimates are not trustworthy! ",
              appendLF = FALSE)
      warn.end = TRUE
      which.run[!which.run == "threelevel.che"] -> which.run
    }
  } else {
    
    if ("threelevel.che" %in% which.run){
      message("- ", crayon::green("[OK] "), 
              "Calculating effect size using 3L-CHE model (rho=", 
              rho.within.study, "; phi=", phi.within.study,
              "; multiarm corr.)... ",
              appendLF = FALSE)
      message(crayon::green("DONE"))
      message("- ", crayon::green("[OK] "), 
              "Robust variance estimation (RVE) used for three-level CHE model... ",
              appendLF = FALSE)
    }
    
    if (mCHE$has.error){
      mCHERes.RVE = 
        data.frame(k = NA, g = NA, g.ci = NA,
                   p = NA, i2 = NA, i2.ci = NA,
                   prediction.ci = NA, nnt = NA,
                   excluded = "none")
    } else {
      # Get results: RVE used
      crit = qt(0.025, mCHE$value$ddf[[1]], lower.tail = FALSE)
      tau2 = sum(mCHE$value$sigma2)
      SE = clubSandwich::conf_int(mCHE$value, vcov = "CR2")[["SE"]]
      pi.lb.rve = clubSandwich::conf_int(mCHE$value, "CR2")[["beta"]] -
        crit * sqrt(tau2 + (SE^2))
      pi.ub.rve = clubSandwich::conf_int(mCHE$value, "CR2")[["beta"]] +
        crit * sqrt(tau2 + (SE^2))
      if (isTRUE(.raw.bin.es)){
        as.numeric(clubSandwich::conf_int(
          mCHE.g$value, "CR2")[["beta"]]) %>% 
          {ifelse(.==0, Inf, metapsyNNT(abs(.), cer))} -> nnt.g
      }
      mCHERes.RVE = with(mCHE$value, {
        data.frame(k = k.all,
                   g = as.numeric(
                     clubSandwich::conf_int(
                       mCHE$value, "CR2")[["beta"]]) %>%
                     ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                     round(round.digits),
                   g.ci = paste0(
                     "[",
                     clubSandwich::conf_int(
                       mCHE$value, "CR2")[["CI_L"]] %>%
                       ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                       round(round.digits), "; ",
                     clubSandwich::conf_int(
                       mCHE$value, "CR2")[["CI_U"]] %>%
                       ifelse(identical(.type.es, "RR"), exp(.), .) %>% 
                       round(round.digits), "]"),
                   p = clubSandwich::coef_test(
                     mCHE$value, "CR2")[["p_Satt"]] %>%
                     scales::pvalue(),
                   i2 = round(I2, 1),
                   i2.ci = ifelse(has.bs, variance.components$i2.ci[3], "-"),
                   prediction.ci = paste0(
                     "[", 
                     round(pi.lb.rve %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "; ",
                     round(pi.ub.rve %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "]"),
                   nnt = ifelse(identical(.type.es, "RR"), 
                                NA, metapsyNNT(abs(as.numeric(
                                  clubSandwich::conf_int(
                                    mCHE$value, "CR2")[["beta"]])), nnt.cer)) %>% 
                     ifelse(isTRUE(.raw.bin.es), nnt.g, .) %>% 
                     round(round.digits) %>% abs(),
                   excluded = "none")
      })
      mCHERes.RVE$excluded = 
        paste0("Number of clusters/studies: ", mCHE$value$s.nlevels[1],
               "; robust variance estimation (RVE) used.")
    }
    
    rownames(mCHERes.RVE) = "Three-Level Model (CHE)"
    mCHERes = mCHERes.RVE
    sendMessage(mCHE, "threelevel.che")
    
    if (length(multi.study) == 0 & "threelevel.che" %in% which.run){
      message("\n- ", crayon::yellow("[!] "), 
              "All included ES seem to be independent.",
              " A three-level CHE model is not adequate and ",
              "tau/I2 estimates are not trustworthy! ",
              appendLF = FALSE)
      warn.end = TRUE
      which.run[!which.run == "threelevel.che"] -> which.run
    }
  }
  
  # If non-PD was found, append error 
  if (append.error){
    mCHE$has.error = TRUE
    mCHE$error$message = paste(vcalc.res$warning) 
  }
  if (has.bs){
    mCHE$value$i2.ci.boot = TRUE
  } else {
    mCHE$value$i2.ci.boot = FALSE
  }
  
  return(list(m = mCHE$value, 
              res = mCHERes,
              has.error = mCHE$has.error,
              message = mCHE$error$message))
}


#' Add `meta::rob()` element to model
#' @keywords internal
addRobData = function(mdl, rob.data) {
  
  # Get data from model
  data = mdl$data
  
  # Perform data input checks
  if (is.null(rob.data$domains)) 
    stop("Please provide 'domains' for the 'rob.data' specification.")
  if (is.null(rob.data$categories)) 
    stop("Please provide 'categories' for the 'rob.data' specification.")
  domains = try({data[rob.data$domains]}, silent=TRUE)
  if (class(domains)[1] == "try-error")
    stop("None of the provided domains could be found in 'data'.")
  if (!is.null(rob.data$domain.names) && 
      length(rob.data$domains)!=length(rob.data$domain.names)) 
    stop("Number of 'domains' and 'domain.names' do not match.")
  overall.rob = try({data[rob.data$overall.rob]}, silent=TRUE)
  if (class(overall.rob)[1] == "try-error")
    stop("'",rob.data$overall.rob, "' could be found in 'data'.")
  
  # Create args list
  args = lapply(as.list(domains), as.character)
  names(args) = paste0("item", 1:length(args))
  args$domains = if (is.null(rob.data$domain.names)) {
    colnames(domains)} else {rob.data$domain.names}
  args$overall = if (ncol(overall.rob)==0) {
    NULL } else { as.character(unlist(overall.rob)) }
  args$categories = rob.data$categories
  args$symbols = if (is.null(rob.data$symbols)) {
    substring(args$categories, 1, 1)} else {rob.data$symbols}
  args$col = rob.data$colors
  args$data = mdl
  args$overwrite = TRUE
  
  # Run rob function
  res = do.call(meta::rob, args)
  
  return(res)
}


#' Set colnames
#' @keywords internal 
setColnames = function (x, value) 
{
  if (is.data.frame(x)) {
    names(x) <- value
  }
  else {
    dn <- dimnames(x)
    if (is.null(dn)) {
      if (is.null(value)) 
        return(x)
      if ((nd <- length(dim(x))) < 2L) 
        stop("attempt to set 'colnames' on an object with less than two dimensions")
      dn <- vector("list", nd)
    }
    if (length(dn) < 2L) 
      stop("attempt to set 'colnames' on an object with less than two dimensions")
    if (is.null(value)) 
      dn[2L] <- list(NULL)
    else dn[[2L]] <- value
    dimnames(x) <- dn
  }
  x
}


#' Compute I-squared from tau2 using the "generalized" formula
#' (as used, e.g., in metafor). This formula uses an estimate
#' of the "typical" within-study variance (v-tilde).
#' @keywords internal 
i2.gen = function(tau2, k, vi){
  wi = 1/vi
  v.tilde = ((k-1)*sum(wi))/
    (sum(wi)^2-sum(wi^2))
  i2 = 100*(tau2/(tau2+v.tilde))
  return(i2)
}


#' Check for problems in data sets prepared for netmeta.
#' @keywords internal 
checkProblemsNMA = function(dat.netmeta){
  c(unique(dat.netmeta[
    dat.netmeta$condition_arm1 == 
      dat.netmeta$condition_arm2, "study"]),
    table(dat.netmeta$study, 
          with(dat.netmeta, 
               paste0(condition_arm1, condition_arm2))) %>% 
      {rowSums(. > 1) > 0} %>% {names(.)[.]}) %>% unique()
}


#' Add all trial arm combinations for multiarm trials in NMA.
#' @importFrom crayon green yellow
#' @importFrom dplyr select ends_with filter group_map group_by all_of any_of
#' @importFrom stringr str_replace_all str_remove_all
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @keywords internal 
addAllCombinations = 
  function(
    data,
    vars.for.id = c("study", "outcome_type",
                    "instrument", "time",
                    "time_weeks", "rating"),
    vars.for.es = c("mean", "sd", "n", "mean_change", "sd_change", "n_change",   
                    "event", "totaln"),
    condition = "condition",
    condition.specification = "multi",
    groups.column.indicator = c("_arm1", "_arm2")){
    
    
    condition.vars = 
      paste0(condition, 
             groups.column.indicator)
    
    es.vars = c(
      paste0(vars.for.es, groups.column.indicator[1]),
      paste0(vars.for.es, groups.column.indicator[2])) %>% 
      {.[order(.)]}
    
    multiarm.vars = 
      paste0(condition.specification, 
             groups.column.indicator)
    
    apply(data, 1,
          function(x){
            paste(as.character(x[vars.for.id]), 
                  collapse = "_")}) -> data$id.multiarm
    
    table(data$id.multiarm[!apply(data, 1, function(x) 
      x[[condition.vars[1]]] == x[[condition.vars[2]]])]) %>% 
      {names(.[.>1])} -> multiarm.studies
    
    str_remove_all(colnames(data), groups.column.indicator[1]) %>% 
      str_remove_all(groups.column.indicator[2]) %>% table() %>% 
      {.[.==2]} %>% names() -> arm.cols
    
    paste0(rep(arm.cols,2), groups.column.indicator) -> arm.var.cols
    
    
    multiarm.list = list()
    for (i in 1:length(multiarm.studies)){
      
      x = data[data$id.multiarm==multiarm.studies[i],]
      trts = unique(unlist(x[,condition.vars]))
      combinations = t(combn(unique(unlist(x[,condition.vars])), 2))
      
      # Part of df that does not have arm info
      df.equal = x %>% 
        select(!all_of(arm.var.cols))
      
      # Part with info on arm 1
      df.arm1 = x %>% 
        select(ends_with(groups.column.indicator[1]))
      df.arm1.1 = df.arm1 %>% 
        {colnames(.) = gsub(
          groups.column.indicator[1], groups.column.indicator[2], colnames(.));.}
      
      # Part with info on arm 2
      df.arm2 = x %>% 
        select(ends_with(groups.column.indicator[2]))
      df.arm2.1 = df.arm2 %>% 
        {colnames(.) = gsub(
          groups.column.indicator[2], groups.column.indicator[1], colnames(.));.}
      
      # Combine
      df.comb = list(
        arm1 = rbind(
          select(df.arm1, any_of(arm.var.cols)), 
          select(df.arm2.1, any_of(arm.var.cols))),
        arm2 = rbind(
          select(df.arm2, any_of(arm.var.cols)), 
          select(df.arm1.1, any_of(arm.var.cols)))
        )
      
      apply(combinations, 1, function(z){
        cbind(
          df.comb$arm1[df.comb$arm1[[condition.vars[1]]] == z[1],][1,],
          df.comb$arm2[df.comb$arm2[[condition.vars[2]]] == z[2],][1,])
      }) %>% 
        do.call(rbind, .) -> df.arms
      
      cbind(df.equal[rep(1, nrow(df.arms)),], df.arms) %>% 
        dplyr::select(any_of(colnames(x))) -> multiarm.list[[i]]
    }
    
    rbind(data[!data$id.multiarm %in% multiarm.studies,],
          do.call(rbind, multiarm.list)) %>% 
      {.[order(.$id.multiarm),]} -> return.data
    
    rownames(return.data) = NULL
    return(return.data)
}



#' Function to generate forest plots with EB estimates
#' @importFrom crayon green yellow
#' @importFrom dplyr select ends_with filter group_map group_by
#' @importFrom stringr str_replace_all str_remove_all
#' @importFrom metafor rma.uni
#' @importFrom clubSandwich conf_int
#' @importFrom stats coef weights
#' @keywords internal 
forestBlup = function(model, which = NULL, col.line = "#a7a9ac",
                      col.polygon = "#6b58a6", leftlab = "Study",
                      rightlab = "g [95% CI]", summarylab = "Total (95% CI)",
                      sort = TRUE, hetstat = TRUE, eb.labels = TRUE){
  
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
    stop("Forest plots with empirical Bayes estimates or not supported for '",
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
           " Bayes plot for a three-level model,",
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
    forestBlupPlotter(dat, res, sort, col.line, col.polygon,
                      hetstat, leftlab, rightlab, summarylab, M.3l,
                      threeLevel, eb.labels)
    title("Lowest")
    
    # Generate "highest" plot
    M = model[models[[which.run[1]]]][[2]]
    if (!is.null(M$exclude)){
      M$data = M$data[!M$exclude,]
    }
    dat = escalc(yi=.TE, sei=.seTE, data = M$data)
    res = metafor::rma(yi = .TE, sei = .seTE, dat = dat, slab = .studlab,
                       method = M$method.tau, test = ifelse(M$hakn, "knha", "z"))
    forestBlupPlotter(dat, res, sort, col.line, col.polygon,
                      hetstat, leftlab, rightlab, summarylab, M.3l,
                      threeLevel, eb.labels)
    title("Highest")
    
  } else {
    
    M = model[models[[which.run[1]]]][[1]]
    if (!is.null(M$exclude)){
      M$data = M$data[!M$exclude,]
    }
    if (threeLevel) {
      M.3l = model[models[[whichThreeLevel[1]]]][[1]]
    }
    dat = escalc(yi=.TE, sei=.seTE, data = M$data)
    res = metafor::rma(yi = .TE, sei = .seTE, dat = dat, slab = .studlab,
                       method = M$method.tau, test = ifelse(M$hakn, "knha", "z"))
  
    if (threeLevel){
      message("- ", crayon::green("[OK] "), 
              "Generating forest plot ('", whichThreeLevel,"' model).")
    } else {
      message("- ", crayon::green("[OK] "), 
              "Generating forest plot ('", which.run,"' model).")
    }
    forestBlupPlotter(dat, res, sort, col.line, col.polygon,
                      hetstat, leftlab, rightlab, summarylab, M.3l,
                      threeLevel, eb.labels)
  }
}


#' Internal helper function to generate forest plots with EB estimates
#' @importFrom crayon green yellow
#' @importFrom dplyr select ends_with filter group_map group_by
#' @importFrom stringr str_replace_all str_remove_all
#' @importFrom stats dffits model.matrix rnorm rstudent coef weights
#' @importFrom utils combn
#' @importFrom metafor rma.uni blup.rma.uni addpoly
#' @importFrom clubSandwich conf_int
#' @importFrom graphics arrows
#' @keywords internal 
#' @details Parts of the code in this function is based on van Aert et al. (2021; Supplement).
#' @references 
#' van Aert, R. C., Schmid, C. H., Svensson, D., & Jackson, D. (2021). 
#' Study specific prediction intervals for random-effects meta-analysis: A tutorial: 
#' Prediction intervals in meta-analysis. _Research Synthesis Methods, 12_(4), 429-447.
forestBlupPlotter = function(dat, res, sort, col.line, col.polygon,
                             hetstat, leftlab, rightlab, summarylab, M.3l,
                             threeLevel, eb.labels){
  
  
  k = nrow(dat)
  psize = weights(res)
  psize = 1.2 + (psize - min(psize)) / (max(psize) - min(psize))
  
  # Create forest plot
  if (sort[1]) { order = "yi" } else { order = NULL }
  sav = metafor::forest(
    dat$yi, sei = dat$.seTE, ylim=c(-0.5,k), cex=0.88,
    pch=18, psize=psize, efac=0, refline=NA, lty=c(1,0), xlab="",
    rowadj=-.07, slab = dat$study, order = order)
  
  # Create summary data escalc
  sumdat = summary(dat)
  order.te = 1:nrow(dat)
  if (sort[1]){
    sumdat = summary(dat)[order(summary(dat)$.TE),]
    order.te = order(summary(dat)$.TE)
  }
  
  # Add annotations
  segments(0, -1, 0, k+1.6, col=col.line)
  if (threeLevel){
    segments(coef(M.3l), 0, coef(M.3l), k, col=col.polygon, lty="33", lwd=0.8)
  } else {
    segments(coef(res), 0, coef(res), k, col=col.polygon, lty="33", lwd=0.8)
  }
  segments(sumdat$ci.lb, k:1, sumdat$ci.ub, k:1, col=col.polygon, lwd=1.5)
  points(sumdat$yi, k:1, pch=18, cex=psize*1.15, col="white")
  points(sumdat$yi, k:1, pch=18, cex=psize, col=col.polygon)
  axis(side=1, at=seq(-1,1,by=0.5), col=col.line, labels=FALSE)
  par(xpd=NA)
  par(cex=sav$cex, font=2)
  text(sav$xlim[1], k+.75, pos=4, leftlab)
  text(sav$xlim[2], k+.75, pos=2, rightlab)
  par(cex=sav$cex, font=1)
  if (hetstat[1]){
    text(sav$xlim[1], -1, pos=4,
         bquote(paste("Test for heterogeneity: ", tau^2, "=",
                      .(formatC(res$tau2, digits=2, format="f")), "; ", italic(I)^2, "=",
                      .(formatC(res$I2, digits=0, format="f")), "%")))
  }
  
  efac = 1
  rows = sav$rows
  blups = metafor::blup.rma.uni(res)[order.te,]
  
  # Calculate BLUPs
  arrows(x0 = blups$pi.lb, x1 = blups$pi.ub, y0 = rev(rows-0.3), code = 3,
         angle = 90, length = 0.02*efac[1], lty = 1, col="darkgray")
  points(x = blups$pred, y = rev(rows-0.3),
         cex = psize*0.5, bg="lightgray", col="darkgray", pch=21)
  
  ### Add prediction interval for predicted true effect size
  if (threeLevel){
    re_lb = M.3l$b[1] - qt(1-.025, df = M.3l$k-2) * sqrt(M.3l$se^2 + M.3l$sigma2[1])
    re_ub = M.3l$b[1] + qt(1-.025, df = M.3l$k-2) * sqrt(M.3l$se^2 + M.3l$sigma2[1])
  } else {
    re_lb = res$b[1] - qt(1-.025, df = res$k-2) * sqrt(res$se^2 + res$tau2)
    re_ub = res$b[1] + qt(1-.025, df = res$k-2) * sqrt(res$se^2 + res$tau2)
  }
  arrows(x0 = re_lb, x1 = re_ub, y0 = sav$ylim[1]+0.5, code = 3,
         angle = 90, length = 0.02*efac[1], lty = 1, col = "darkgray")
  
  if (eb.labels[1]){
    labels <- paste0(sprintf("%.2f", blups$pred), " [",
                     sprintf("%.2f", blups$pi.lb), ", ",
                     sprintf("%.2f", blups$pi.ub), "]")
    text(x = sav$xlim[2], y = rev(rows-0.35), labels = labels,
         pos = 2, col = "darkgray", cex=0.88)
    
    # Add prediction interval for predicted true effect size in numbers
    label_theta = paste0(sprintf("%.2f", ifelse(threeLevel, 
                                                M.3l$b[1], res$b[1])), " [",
                         sprintf("%.2f", re_lb), ", ",
                         sprintf("%.2f", re_ub), "]")
    text(x = sav$xlim[2], y = sav$ylim[1]+0.2, labels = label_theta,
         pos = 2, col = "gray", cex=0.88)
  }
  
  # Add color segments
  if (threeLevel){
    metafor::addpoly(x = coef(M.3l),
                     ci.lb = clubSandwich::conf_int(M.3l, "CR2")[,5],
                     ci.ub = clubSandwich::conf_int(M.3l, "CR2")[,6],
                     row = 0, efac = 2, col=col.polygon, border=col.polygon,
                     mlab=summarylab)
  } else {
    metafor::addpoly(res, row=0, mlab=summarylab, efac=2,
                     col=col.polygon, border=col.polygon)
  }
  
}


#' Best Linear Unbiased Predictions (BLUPs) for 'runMetaAnalysis' models.
#' 
#' Generates empirical Bayes (EB) estimates, also known as best linear unbiased 
#' predictions (BLUPs), by merging the fitted values obtained from fixed effects 
#' and estimated contributions of random effects. These estimates represent 
#' the study-specific true effect sizes or outcomes and are accompanied by 
#' standard errors and prediction interval bounds.
#'
#'
#' @title eb: Empirical Bayes estimates
#' @param x Model
#' @param ... Other arguments
#' @export
eb <- function (x, ...) {
  UseMethod("eb", x)
}


#' Best Linear Unbiased Predictions (BLUPs) for 'runMetaAnalysis' models.
#' 
#' Generates empirical Bayes (EB) estimates, also known as best linear unbiased 
#' predictions (BLUPs), by merging the fitted values obtained from fixed effects 
#' and estimated contributions of random effects. These estimates represent 
#' the study-specific true effect sizes or outcomes and are accompanied by 
#' standard errors and prediction interval bounds.
#'
#'
#' @title blup: Empirical Bayes estimates
#' @param x Model
#' @param ... Other arguments
#' @export
blup <- function (x, ...) {
  UseMethod("blup", x)
}



#' Catch variables from parent frame
#' @keywords internal 
catch <- function(argname, matchcall, data, encl) {
  eval(matchcall[[match(argname, names(matchcall))]], data, enclos = encl)
}

#' Catch variable names from parent frame
#' @keywords internal 
catchName <- function(argname, matchcall, data, encl) {
  as.character(matchcall[[match(argname, names(matchcall))]])
}


#' Test if the randomized proportion differs from the original allocation ratio
#' 
#' This function returns a one-sample proportion test that allows to 
#' examine if the included sample in a study deviated significantly form
#' the intended allocation ratio.
#' 
#' @param n_arm1 A vector in which the sample size of the first trial arm is stored.
#' @param n_arm2 A vector in which the sample size of the second trial arm is stored.
#' @param study An optional vector with study labels.
#' @param data An optional `data.frame` in which the sample size data is included.
#' @param ratio A single string or character vector of the expected allocation ratio.
#' This ratio should be separated by a colon, e.g. `"1:3`, where the left side
#' indicates the first trial arm (`n_arm1`) and the right side indicates the
#' second trial arm (`n_arm2`). Defaults to `"1:1`.
#' @param ... Additional arguments.
#' 
#' @examples 
#' \dontrun{
#' # Load example data that follows the Metapsy data standard
#' data("depressionPsyCtr")
#' 
#' # This is an unexported function, so we need the metapsyTools::: prefix
#' # Test if all studies follow the 1:1 allocation ratio
#' metapsyTools:::testRandomizedProportion(
#'   n_arm1, n_arm2, data = depressionPsyCtr, study = study) 
#' 
#' # Provide a comparison-specific allocation ratio
#' ratio <- c("1:1", "3:1", "4:5")
#' metapsyTools:::testRandomizedProportion(
#'   n_arm1, n_arm2, data = depressionPsyCtr[1:3,], 
#'   study = study, ratio = ratio) 
#' }
#' 
#' @importFrom dplyr select ends_with filter group_map group_by tibble
#' @importFrom stringr str_replace_all str_remove_all
#' @importFrom stats dffits model.matrix rnorm rstudent p.adjust
#' @importFrom utils combn
#' @keywords internal 
testRandomizedProportion <- 
  function(n_arm1, n_arm2, study, data = NULL, ratio = "1:1", ...){
    
    # get data; extract from parent frame if necessary
    isNullData <- is.null(data)
    parentFrame <- sys.frame(sys.parent())
    mC <- match.call()
    if (isNullData) { data <- parentFrame}
    
    # get data & check missings
    n_arm1 <- try({catch("n_arm1", mC, data, parentFrame)}, silent = TRUE)
    n_arm2 <- try({catch("n_arm2", mC, data, parentFrame)}, silent = TRUE)
    if (!missing(study)) {
      study <- try({catch("study", mC, data, parentFrame)}, silent = TRUE)
      if (identical(class(study), "try-error")) 
        stop(catchName("study", mC, data, parentFrame), 
             " variable not found in data.") 
    } else {
      study = NULL
    }
    if (identical(class(n_arm1), "try-error")) 
      stop(catchName("n_arm1", mC, data, parentFrame), 
           " variable not found in data.") 
    if (identical(class(n_arm2), "try-error"))
      stop(catchName("n_arm2", mC, data, parentFrame),
           " variable not found in data.") 
    
    # check ratio
    if (sum(!grepl(":", ratio)) > 0) 
      stop("'ratio' must be a string including two",
           " numbers separated by a semicolon (e.g. '1:1')") 
    
    # define expected proportion
    e.ratio <- strsplit(ratio, ":")
    e.prop <- unlist(lapply(e.ratio, function(y) 
    {y <- as.numeric(y); (y/sum(y))[1]}))
    
    # calculate observed proportion
    o.prop <- n_arm1/(n_arm1 + n_arm2)
    
    # calculate test of proportion difference
    p.diff <- o.prop - e.prop 
    se <- sqrt(e.prop*(1-e.prop)/(n_arm1 + n_arm2))
    z <- p.diff/se
    p <- pnorm(abs(z), lower.tail = FALSE)*2
    
    # return results
    if (is.null(study)[1])
      study <- 1:length(p.diff)
    
    return(tibble(study = study, n.arm1 = n_arm1, n.arm2 = n_arm2,
                  ratio = ratio, p.expected = e.prop,
                  p.observed = o.prop, p.diff = p.diff,
                  se = se, z = z, p = p))
  }



#' Test for baseline imbalances
#' 
#' This function returns a SMD for baseline imbalance tests, along with
#' 99% confidence intervals. If the `cluster` argument is specified, _p_-values
#' can be adjusted using methods available in [stats::p.adjust()].
#' 
#' @param m_arm1 A vector with the baseline mean(s) in the first arm. 
#' @param m_arm2 A vector with the baseline mean(s) in the second arm. 
#' @param sd_arm1 A vector with the baseline standard deviation(s) in the first arm. 
#' @param sd_arm2 A vector with the baseline standard deviation(s) in the second arm. 
#' @param n_arm1 A vector with the baseline sample size(s) in the first arm. 
#' @param n_arm2 A vector with the baseline sample size(s) in the second arm. 
#' @param cluster An optional vector of the same length as `m_arm1`, `m_arm2`, etc.,
#' encoding to which larger cluster/study a comparison belongs to. If specified,
#' _p_ values of the comparisons will be adjusted according to the method provided in
#' the `p.adj` argument.
#' @param study An optional vector with study labels.
#' @param data An optional `data.frame` that includes the effect size data.
#' @param p.adj A character string, specifying the type of _p_-value adjustment. For available
#' options, see the "Details" section in [stats::p.adjust()]. Defaults to `"holm"`.
#' @param ... Additional arguments.
#' 
#' @examples 
#' \dontrun{
#' # Load example data that follows the Metapsy data standard
#' data("depressionPsyCtr")
#' 
#' # This is an unexported function, so we need the metapsyTools::: prefix
#' # Test for differences without p-value adjustment
#' metapsyTools:::testBaselineImbalance(
#'   bl_mean_arm1, bl_mean_arm2, 
#'   bl_sd_arm1, bl_sd_arm2, 
#'   bl_n_arm1, bl_n_arm2, 
#'   data = depressionPsyCtr,
#'   study = study)
#' 
#' # Provide cluster variable: p-values are adjusted within clusters
#' metapsyTools:::testBaselineImbalance(
#'   bl_mean_arm1, bl_mean_arm2, 
#'   bl_sd_arm1, bl_sd_arm2, 
#'   bl_n_arm1, bl_n_arm2, 
#'   data = depressionPsyCtr,
#'   cluster = study,
#'   study = study)
#' 
#' # Change adjustment method
#' metapsyTools:::testBaselineImbalance(
#'   bl_mean_arm1, bl_mean_arm2, 
#'   bl_sd_arm1, bl_sd_arm2, 
#'   bl_n_arm1, bl_n_arm2, 
#'   data = depressionPsyCtr,
#'   cluster = study,
#'   study = study, p.adj = "bonferroni")
#' }
#' 
#' @importFrom dplyr select ends_with filter group_map group_by tibble
#' @importFrom stringr str_replace_all str_remove_all
#' @importFrom stats dffits model.matrix rnorm rstudent p.adjust
#' @importFrom utils combn
#' @keywords internal 
testBaselineImbalance <- 
  function(m_arm1, m_arm2, sd_arm1, sd_arm2, n_arm1, n_arm2, 
           cluster, study, data = NULL, p.adj = "holm", ...){
    
    # get data; extract from parent frame if necessary
    isNullData <- is.null(data)
    parentFrame <- sys.frame(sys.parent())
    mC <- match.call()
    if (isNullData) { data <- parentFrame}
    
    
    # get data & check missings
    m_arm1 <- try({catch("m_arm1", mC, data, parentFrame)}, silent = TRUE)
    m_arm2 <- try({catch("m_arm2", mC, data, parentFrame)}, silent = TRUE)
    sd_arm1 <- try({catch("sd_arm1", mC, data, parentFrame)}, silent = TRUE)
    sd_arm2 <- try({catch("sd_arm2", mC, data, parentFrame)}, silent = TRUE)
    n_arm1 <- try({catch("n_arm1", mC, data, parentFrame)}, silent = TRUE)
    n_arm2 <- try({catch("n_arm2", mC, data, parentFrame)}, silent = TRUE)
    
    if (!missing(study)) {
      study <- try({catch("study", mC, data, parentFrame)}, silent = TRUE)
      if (identical(class(study), "try-error")) 
        stop(catchName("study", mC, data, parentFrame), 
             " variable not found in data.") 
    } else {
      study = NULL 
    }
    
    if (!missing(cluster)) {
      cluster <- try({catch("cluster", mC, data, parentFrame)}, silent = TRUE)
      if (identical(class(cluster), "try-error")) 
        stop(catchName("cluster", mC, data, parentFrame), 
             " variable not found in data.") 
    } else {
      cluster = NULL 
    }
    
    if (identical(class(m_arm1), "try-error")) 
      stop(catchName("m_arm1", mC, data, parentFrame), 
           " variable not found in data.") 
    if (identical(class(n_arm2), "try-error"))
      stop(catchName("m_arm2", mC, data, parentFrame),
           " variable not found in data.") 
    if (identical(class(sd_arm1), "try-error")) 
      stop(catchName("sd_arm1", mC, data, parentFrame), 
           " variable not found in data.") 
    if (identical(class(sd_arm2), "try-error"))
      stop(catchName("sd_arm1", mC, data, parentFrame),
           " variable not found in data.") 
    if (identical(class(n_arm1), "try-error")) 
      stop(catchName("n_arm1", mC, data, parentFrame), 
           " variable not found in data.") 
    if (identical(class(n_arm2), "try-error"))
      stop(catchName("n_arm1", mC, data, parentFrame),
           " variable not found in data.") 
    
    # Calculate SMDs
    sdp <- sqrt(((n_arm1-1)*sd_arm1^2 + (n_arm2-1)*sd_arm2^2)/
                  ((n_arm1-1) + (n_arm2-1)))
    smd <- (m_arm1 - m_arm2) / sdp
    se <- sqrt(((n_arm1 + n_arm2)/(n_arm1 * n_arm2)) + 
                 (smd^2/(2*(n_arm1+n_arm2))))
    smd.99lower <- smd - (qnorm(1-.005) * se)
    smd.99upper <- smd + (qnorm(1-.005) * se)
    z <- smd/se
    p <- pnorm(abs(z), lower.tail = FALSE)*2
    
    # Encode study if missing
    if (is.null(study)[1])
      study <- 1:length(smd)
    
    
    # Adjust p-values
    if (!is.null(cluster)[1]) {
      message("adjusting p-values based on '", 
              catchName("cluster", mC, data, parentFrame),
              "' variable (", p.adj, " method).")
      adj.p.l <- lapply(split(p, cluster), p.adjust, p.adj)
      adj.p <- setNames(unlist(adj.p.l, use.names=FALSE), 
                        rep(names(adj.p.l), lengths(adj.p.l)))
      
      p.adjusted <- p
      for (i in 1:length(unique(cluster))){
        p.adjusted[which(cluster %in% unique(cluster)[i])] <- 
          adj.p[which(names(adj.p) %in% unique(cluster)[i])]
      }
      
      return(tibble(study = study, smd = smd, 
                    smd.99lower = smd.99lower,
                    smd.99upper = smd.99upper,
                    se = se, z = z, cluster = cluster,
                    p = p, p.adjusted = p.adjusted))
      
    }
    
    return(tibble(study = study, smd = smd, 
                  smd.99lower = smd.99lower,
                  smd.99upper = smd.99upper,
                  se = se, z = z, p = p))
  }


#' Import 'update.meta' for internal use
#' 
#' This function imports the `update.meta` function of the `meta` package,
#' which has ceased being an exported function since version 7.0-0.
#' 
#' @param ... Arguments forwarded to `update.meta`.
#' 
#' @import meta utils
#' @keywords internal 

updateMeta = function(...) {
  func = getFromNamespace("update.meta", "meta")
  dots = list(...)
  do.call(func, dots)
}


#' Create a RevMan-style risk of bias summary chart
#'
#' This function generates summary plots for study quality assessments using the
#' \href{https://bit.ly/2KGQtfG}{Cochrance Risk of Bias Tool}.
#' Summary plots follow the style of \href{https://bit.ly/30eJK29}{RevMan} Risk of Bias (RoB) summary charts.
#'
#' @usage robSummary(data, name.high="High", name.unclear="Unclear",
#'     name.low="Low", studies, name.missing, table = FALSE)
#'
#' @param data A \code{data.frame} containing a column for each risk of bias criterion, where
#' rows represent each individual studies. The risk of bias assessment for each criterion in each
#' study must be coded as a character string. Up to four codes can be used, referring to low risk of bias,
#' unclear risk of bias, high risk of bias, or missing information. The string used to specify the categories
#' must be specified in \code{name.high}, \code{name.unclear}, \code{name.low} and/or \code{name.missing},
#' unless defaults for those parameters are used.
#' @param name.high Character specifying how the "high risk of bias" category was coded in \code{data}
#' (e.g., \code{name.high = "high"}). Default is \code{"High"}.
#' @param name.unclear Character specifying how the "unclear risk of bias" category was coded in \code{data}
#' (e.g., \code{name.unclear = "unclear"}). Default is \code{"Unclear"}.
#' @param name.low Character specifying how the "low risk of bias" category was coded in \code{data}
#' (e.g., \code{name.low = "low"}). Default is \code{"Low"}.
#' @param name.missing Character specifying how missing information was coded in \code{data}
#' (e.g., \code{name.missing} = \code{"missing"}). Default is \code{"Missing"}. All ratings, including missing
#' information, must be coded as strings, so using \code{NA} in \code{data} to signify missing information
#' is not valid.
#' @param studies A vector of the same length as the number of rows in \code{data} specifying the study
#' labels for the risk of bias ratings. Only has to be specified when \code{table = TRUE}.
#' @param table Should an additional RevMan style risk of bias table be produced? If set to \code{TRUE},
#' \code{studies} must be specified. \code{FALSE} by default.
#'
#' @details The function automatically removes separators like "-" or "." from column names/risk of bias criteria. To produce
#' a "clean" plot, you may therefore separate words in the column names of the \code{data} data frame using these
#' symbols (e.g. \code{"Allocation_Concealment"} to return "Allocation Concealment").
#'
#' @references Harrer, M., Cuijpers, P., Furukawa, T.A, & Ebert, D. D. (2019).
#' \emph{Doing Meta-Analysis in R: A Hands-on Guide}. DOI: 10.5281/zenodo.2551803.
#' \href{https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/creating-a-revman-style-risk-of-bias-summary.html}{Chapter 10}
#'
#' @author Mathias Harrer & David Daniel Ebert
#'
#' @import ggplot2
#'
#' @keywords internal

robSummary = function (data, name.high = "High", name.unclear = "Unclear",
                        name.low = "Low", studies, name.missing, table = FALSE) {
  
  
  if (class(data) != "data.frame") {
    stop("'data' must be of class 'data.frame'.")
  }
  if (missing(name.missing)) {
    colnames.rob = character()
    for (i in 1:ncol(data)) {
      vect = as.character(data[, i])
      for (j in 1:length(data[, i])) {
        if (vect[j] %in% c(name.high, name.unclear, name.low)) {
          colnames.rob[i] = TRUE
        }
        else {
          colnames.rob[i] = FALSE
          message(cat("Column '", colnames(data)[i],
                      "' removed from plot because it did not contain the specified RoB ratings (only). \n",
                      sep = ""))
          break
        }
      }
    }
    rob = data[, as.logical(colnames.rob)]
    for (i in 1:ncol(rob)) {
      rob[, i] = as.character(rob[, i])
      rob[rob[, i] == name.high, i] = "High"
      rob[rob[, i] == name.unclear, i] = "Unclear"
      rob[rob[, i] == name.low, i] = "Low"
    }
    if (table == TRUE) {
      if (missing(studies)) {
        stop("'studies' has to be specified when 'table = TRUE'.")
      }
      if (length(as.vector(studies)) != nrow(data)) {
        stop("'studies' vector is not of equal length as the data.")
      }
      if (length(unique(studies)) != length(studies)) {
        stop("'studies' cannot contain duplicate study labels.")
      }
      robby = rob
      robby = data.frame(study = studies, condition = rep(colnames(robby),
                                                          each = length(studies)), measurement = unlist(robby))
      rownames(robby) = NULL
      robby$condition = gsub("_", " ", robby$condition)
      robby$condition = gsub("-", " ", robby$condition)
      robby$condition = gsub("\\.", " ", robby$condition)
      robby[robby$measurement == "Low", "measurement"] = "+"
      robby[robby$measurement == "Unclear", "measurement"] = "?"
      robby[robby$measurement == "High", "measurement"] = "-"
      robby$study = factor(robby$study, levels = unique(studies)[rev(order(unique(robby$study)))])
      rob.table = ggplot(data = robby, aes(y = study, x = condition)) +
        geom_tile(color = "black", fill = "white", size = 0.8) +
        geom_point(aes(color = as.factor(measurement)),
                   size = 20) + geom_text(aes(label = measurement),
                                          size = 8) + scale_x_discrete(position = "top") +
        scale_color_manual(values = c(`?` = "#E2DF07",
                                      `-` = "#BF0000", `+` = "#02C100")) + theme_minimal() +
        coord_equal() + theme(axis.title.x = element_blank(),
                              axis.title.y = element_blank(), axis.ticks.y = element_blank(),
                              axis.text.y = element_text(size = 15, color = "black"),
                              axis.text.x = element_text(size = 13, color = "black",
                                                         angle = 90, hjust = 0), legend.position = "none",
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank())
    }
    rob.long = data.frame(condition = rep(colnames(rob),
                                          each = nrow(rob)), measurement = unlist(rob))
    rownames(rob.long) = NULL
    rob.long$condition = gsub("_", " ", rob.long$condition)
    rob.long$condition = gsub("-", " ", rob.long$condition)
    rob.long$condition = gsub("\\.", " ", rob.long$condition)
    rob.long$measurement = as.factor(rob.long$measurement)
    rob.long$measurement = factor(rob.long$measurement,
                                  levels(rob.long$measurement)[c(1, 3, 2)])
    rob.plot = ggplot(data = rob.long) +
      geom_bar(mapping = aes(x = condition, fill = measurement),
               width = 0.7, position = "fill", color = "black") +
      coord_flip(ylim = c(0, 1)) +
      guides(fill = guide_legend(reverse = TRUE)) +
      scale_fill_manual("Risk of Bias",
                        labels = c(`High` = "    High risk of bias          ",
                                   `Unclear` = "    Unclear risk of bias       ",
                                   `Low` = "    Low risk of bias  "),
                        values = c(High = "#BF0000",
                                   Unclear = "#E2DF07",
                                   Low = "#02C100")) +
      scale_y_continuous(labels = scales::percent) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
            axis.ticks.y = element_blank(), axis.text.y = element_text(size = 18,
                                                                       color = "black"), axis.line.x = element_line(colour = "black",
                                                                                                                    linewidth = 0.5, linetype = "solid"), legend.position = "bottom",
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), legend.background = element_rect(linetype = "solid",
                                                                                 colour = "black"), legend.title = element_blank(),
            legend.key.size = unit(0.75, "cm"), legend.text = element_text(size = 14))
    plot(rob.plot)
    if (table == TRUE) {
      plot(rob.table)
    }
  }
  else {
    data = as.data.frame(data)
    colnames.rob = character()
    for (i in 1:ncol(data)) {
      vect = as.character(data[, i])
      for (j in 1:length(data[, i])) {
        if (vect[j] %in% c(name.high, name.unclear, name.low,
                           name.missing)) {
          colnames.rob[i] = TRUE
        }
        else {
          colnames.rob[i] = FALSE
          message(cat("Column '", colnames(data)[i],
                      "' removed from plot because it did not contain the specified RoB ratings (only). \n",
                      sep = ""))
          break
        }
      }
    }
    rob = data[, as.logical(colnames.rob)]
    for (i in 1:ncol(rob)) {
      rob[, i] = as.character(rob[, i])
      rob[rob[, i] == name.high, i] = "High"
      rob[rob[, i] == name.unclear, i] = "Unclear"
      rob[rob[, i] == name.low, i] = "Low"
      rob[rob[, i] == name.missing, i] = "Missing"
    }
    if (table == TRUE) {
      if (missing(studies)) {
        stop("'studies' has to be specified when 'table = TRUE'.")
      }
      if (length(as.vector(studies)) != nrow(data)) {
        stop("'studies' vector is not of equal length as the data.")
      }
      robby = rob
      robby = data.frame(study = as.factor(studies), condition = rep(colnames(robby),
                                                                     each = length(studies)), measurement = unlist(robby))
      rownames(robby) = NULL
      robby$condition = gsub("_", " ", robby$condition)
      robby$condition = gsub("-", " ", robby$condition)
      robby$condition = gsub("\\.", " ", robby$condition)
      robby[robby$measurement == "Low", "measurement"] = "+"
      robby[robby$measurement == "Unclear", "measurement"] = "?"
      robby[robby$measurement == "High", "measurement"] = "-"
      robby[robby$measurement == "Missing", "measurement"] = " "
      robby$study = factor(robby$study, levels = unique(studies)[rev(order(unique(robby$study)))])
      rob.table = ggplot(data = robby, aes(y = study, x = condition)) +
        geom_tile(color = "black", fill = "white", size = 0.8) +
        geom_point(aes(color = as.factor(measurement)),
                   size = 20) + geom_text(aes(label = measurement),
                                          size = 8) + scale_x_discrete(position = "top") +
        scale_color_manual(values = c(`?` = "#E2DF07",
                                      `-` = "#BF0000", `+` = "#02C100", ` ` = "white")) +
        theme_minimal() + coord_equal() + theme(axis.title.x = element_blank(),
                                                axis.title.y = element_blank(), axis.ticks.y = element_blank(),
                                                axis.text.y = element_text(size = 15, color = "black"),
                                                axis.text.x = element_text(size = 13, color = "black",
                                                                           angle = 90, hjust = 0), legend.position = "none",
                                                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                panel.background = element_blank())
    }
    rob.long = data.frame(condition = rep(colnames(rob),
                                          each = nrow(rob)), measurement = unlist(rob))
    rownames(rob.long) = NULL
    rob.long$condition = gsub("_", " ", rob.long$condition)
    rob.long$condition = gsub("-", " ", rob.long$condition)
    rob.long$condition = gsub("\\.", " ", rob.long$condition)
    rob.long$measurement = as.factor(rob.long$measurement)
    rob.long$measurement = factor(rob.long$measurement, levels(rob.long$measurement)[c(3,
                                                                                       1, 4, 2)])
    rob.plot = ggplot(data = rob.long) + geom_bar(mapping = aes(x = condition,
                                                                fill = measurement), width = 0.7, position = "fill",
                                                  color = "black") + coord_flip(ylim = c(0, 1)) + guides(fill = guide_legend(reverse = TRUE)) +
      scale_fill_manual("Risk of Bias", labels = c(`High` = "    High risk of bias          ",
                                                   `Unclear` = "    Unclear risk of bias       ",
                                                   `Low` = "    Low risk of bias  ",
                                                   `Missing` = "    Missing information      "),
                        values = c(Missing = "white",
                                   High = "#BF0000",
                                   Unclear = "#E2DF07",
                                   Low = "#02C100")) +
      scale_y_continuous(labels = scales::percent) + 
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(), axis.ticks.y = element_blank(),
            axis.text.y = element_text(size = 18, color = "black"),
            axis.line.x = element_line(colour = "black", linewidth = 0.5,
                                       linetype = "solid"), legend.position = "bottom",
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), legend.background = 
              element_rect(linetype = "solid",
              colour = "black"), legend.title = element_blank(),
            legend.key.size = unit(0.75, "cm"), legend.text = element_text(size = 14))
    plot(rob.plot)
    if (table == TRUE) {
      plot(rob.table)
    }
  }
}


