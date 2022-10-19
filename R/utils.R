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
                           nnt.cer, which.run){
  
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
      fixed = ifelse(method.tau == "FE", TRUE, FALSE),
      random = ifelse(method.tau == "FE", FALSE, TRUE)) %>% 
      append(selectArguments(meta::metabin, dots))
    
    mGeneral = 
      tryCatch2(do.call(meta::metabin, mGeneralArgs))
    
    if (!mGeneral$has.error){
      with(meta::update.meta(mGeneral$value, sm="RD"), {
        ifelse(isTRUE(fixed) & isTRUE(random),
               abs(TE.fixed)^-1, abs(TE.random)^-1)
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
                          .raw.bin.es, nnt.cer){
  
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
        meta::update.meta(mGeneral$m, exclude = !lowest, id = NULL))
      if (!mLowest$has.error){
        with(meta::update.meta(mLowest$value, sm="RD"), {
          ifelse(isTRUE(fixed) & isTRUE(random),
                 abs(TE.fixed)^-1, abs(TE.random)^-1)
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
  return(mLowest)
}


#' Fit 'highest' model
#' @keywords internal 
fitHighestModel = function(data, study.var, multi.study,
                           mGeneral, .type.es, round.digits,
                           .raw.bin.es, nnt.cer){
  
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
        meta::update.meta(mGeneral$m, exclude = !highest, id = NULL))
      
      if (!mHighest$has.error){
        with(meta::update.meta(mHighest$value, sm="RD"), {
          ifelse(isTRUE(fixed) & isTRUE(random),
                 abs(TE.fixed)^-1, abs(TE.random)^-1)
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
  return(mHighest)
}


#' Fit 'combined' model
#' @keywords internal 
fitCombinedModel = function(which.combine, which.combine.var,
                            data, study.var, multi.study, es.var, se.var,
                            mGeneral, .type.es, round.digits, hakn,
                            .raw.bin.es, nnt.cer, rho.within.study,
                            method.tau, method.tau.ci, dots,
                            es.binary.raw.vars, phi.within.study){
  
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
        fixed = ifelse(method.tau == "FE", TRUE, FALSE),
        random = ifelse(method.tau == "FE", FALSE, TRUE))
      
      do.call(meta::metabin, 
              mCombArgsBin) %>% 
        with(., {
          ifelse(isTRUE(fixed) & isTRUE(random),
                 abs(TE.fixed)^-1, abs(TE.random)^-1)
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
                                near.pd){
  
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
        fixed = ifelse(method.tau == "FE", TRUE, FALSE),
        random = ifelse(method.tau == "FE", FALSE, TRUE))
      
      do.call(meta::metabin, 
              mCombArgsBin) %>% 
        with(., {
          ifelse(isTRUE(fixed) & isTRUE(random),
                 abs(TE.fixed)^-1, abs(TE.random)^-1)
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
                            m.for.outliers){
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
    with(meta::update.meta(mOutliers$value, 
                           sm = "RD", id = NULL), {
                             ifelse(isTRUE(fixed) & isTRUE(random),
                                    abs(TE.fixed)^-1, abs(TE.random)^-1)
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
                             nnt.cer){
  
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
      meta::update.meta(
        m.for.influence,
        exclude = influenceRes$Data$is.infl == "yes",
        id = NULL))
  
  if (isTRUE(.raw.bin.es) &&
      !mInfluence$has.error){
    with(meta::update.meta(mInfluence$value, sm="RD", id = NULL), {
      ifelse(isTRUE(fixed) & isTRUE(random),
             abs(TE.fixed)^-1, abs(TE.random)^-1)
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
                       nnt.cer){
  
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
        meta::update.meta(
          mGeneral$m, exclude = !robMask,id = NULL))
      which.run[!which.run == "rob"] -> which.run
      warn.end = TRUE
    } else {
      mRob = tryCatch2(
        meta::update.meta(
          m.for.rob, exclude = !robMask, id = NULL))
    }
  } else {
    robMask = rep(TRUE, nrow(mGeneral$m$data))
    mRob = tryCatch2(
      meta::update.meta(
        mGeneral$m, exclude = !robMask, id = NULL))
  }
  if (isTRUE(.raw.bin.es) &&
      !mRob$has.error) {
    with(meta::update.meta(
      mRob$value, sm="RD",
      id = NULL), {
        ifelse(isTRUE(fixed) & isTRUE(random),
               abs(TE.fixed)^-1, abs(TE.random)^-1)
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
                              use.rve){
  
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
                metapsyNNT(., cer))} -> nnt.g
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
                   i2.ci = "-",
                   prediction.ci = paste0("[", 
                                          round(predict(mThreeLevel$value)$pi.lb %>% 
                                                  ifelse(identical(.type.es, "RR"), exp(.), .), 
                                                round.digits), "; ",
                                          round(predict(mThreeLevel$value)$pi.ub %>% 
                                                  ifelse(identical(.type.es, "RR"), exp(.), .), 
                                                round.digits), "]"),
                   nnt = ifelse(identical(.type.es, "RR"), 
                                NA, metapsyNNT(as.numeric(b[,1]), nnt.cer)) %>% 
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
          {ifelse(.==0, Inf, metapsyNNT(., cer))} -> nnt.g
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
                                    mThreeLevel$value, "CR2")[["beta"]]), 
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
      message("- ", crayon::yellow("[!] "), 
              "All included ES seem to be independent.",
              " A three-level model is not adequate and ",
              " tau/I2 estimates are not trustworthy!")
      warn.end = TRUE
      which.run[!which.run == "threelevel"] -> which.run
    }
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
                                 use.rve, rho.within.study, phi.within.study){
  
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
                metapsyNNT(., cer))} -> nnt.g
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
                   i2.ci = "-",
                   prediction.ci = paste0(
                     "[", 
                     round(predict(mCHE$value)$pi.lb %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "; ",
                     round(predict(mCHE$value)$pi.ub %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "]"),
                   nnt = ifelse(identical(.type.es, "RR"), 
                                NA, metapsyNNT(as.numeric(b[,1]), nnt.cer)) %>% 
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
          {ifelse(.==0, Inf, metapsyNNT(., cer))} -> nnt.g
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
                                    mCHE$value, "CR2")[["beta"]]), nnt.cer)) %>% 
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
      message("- ", crayon::yellow("[!] "), 
              "All included ES seem to be independent.",
              " A three-level CHE model is not adequate and ",
              "tau/I2 estimates are not trustworthy!")
      warn.end = TRUE
      which.run[!which.run == "threelevel.che"] -> which.run
    }
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
                                  near.pd){
  
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
                metapsyNNT(., cer))} -> nnt.g
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
                   i2.ci = "-",
                   prediction.ci = paste0(
                     "[", 
                     round(predict(mCHE$value)$pi.lb %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "; ",
                     round(predict(mCHE$value)$pi.ub %>% 
                             ifelse(identical(.type.es, "RR"), exp(.), .), 
                           round.digits), "]"),
                   nnt = ifelse(identical(.type.es, "RR"), 
                                NA, metapsyNNT(as.numeric(b[,1]), nnt.cer)) %>% 
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
          {ifelse(.==0, Inf, metapsyNNT(., cer))} -> nnt.g
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
                                    mCHE$value, "CR2")[["beta"]]), nnt.cer)) %>% 
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
      message("- ", crayon::yellow("[!] "), 
              "All included ES seem to be independent.",
              " A three-level CHE model is not adequate and ",
              "tau/I2 estimates are not trustworthy!")
      warn.end = TRUE
      which.run[!which.run == "threelevel.che"] -> which.run
    }
  }
  
  # If non-PD was found, append error 
  if (append.error){
    mCHE$has.error = TRUE
    mCHE$error$message = paste(vcalc.res$warning) 
  }
  
  return(list(m = mCHE$value, 
              res = mCHERes,
              has.error = mCHE$has.error,
              message = mCHE$error$message))
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



