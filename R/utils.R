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
  list(value=value, warning=warn, error=err)
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
           if (length(y) <= 1) { return(c(es = NA, se = NA)) } 
           else {return(y)}}) %>% 
    do.call(rbind, .) -> x.select
  return(x.select)
}


