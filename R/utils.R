#' Calculate Hedges' g using Mean and Standard Deviation.
#'
#' Calculate Hedges' g using Mean and Standard Deviation. Only meant to be used as
#' part of \code{\link{calculateEffectSizes}}.
#'
#' @param x data
#' @param ... Effect size data. Must be \code{Post_M_trt1}, \code{Post_M_trt2}, \code{Post_SD_trt1},
#' \code{Post_SD_trt2}, \code{Post_N_trt1}, \code{Post_N_trt2} (all \code{numeric}).
#' @usage meanSD(x, ...)
#' @importFrom dplyr select
#' @importFrom purrr pmap_dfr
#' @importFrom esc esc_mean_sd
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export meanSD
#' @keywords internal

meanSD = function(x, ...){
  x %>%
    purrr::pmap_dfr(function(Post_M_trt1, Post_M_trt2, Post_SD_trt1,
                             Post_SD_trt2, Post_N_trt1, Post_N_trt2, ...)
    {esc::esc_mean_sd(Post_M_trt1, Post_SD_trt1, Post_N_trt1,
                      Post_M_trt2, Post_SD_trt2, Post_N_trt2, es.type = "g") %>%
        as.data.frame() %>% dplyr::select(es, se) %>%
        suppressWarnings()})
}

#' Calculate Hedges' g binary outcome data.
#'
#' Calculate Hedges' g binary outcome data. Only meant to be used as
#' part of \code{\link{calculateEffectSizes}}.
#'
#' @param x data
#' @param ... Binary effect size data. Must be \code{Improved_N_trt1}, \code{Improved_N_trt2}, \code{Rand_N_trt1},
#' \code{Rand_N_trt2} (all \code{numeric}).
#' @usage binaryES(x, ...)
#' @importFrom dplyr select mutate
#' @importFrom purrr pmap_dfr
#' @importFrom esc esc_2x2
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export binaryES
#' @keywords internal


binaryES = function(x, ...){
  x %>%
    purrr::pmap_dfr(function(Improved_N_trt1, Improved_N_trt2, Rand_N_trt1, Rand_N_trt2, ...)
    {esc::esc_2x2(Improved_N_trt1,
                  Rand_N_trt1 - Improved_N_trt1,
                  Improved_N_trt2,
                  Rand_N_trt2 - Improved_N_trt2,
                  es.type = "g") %>%
        as.data.frame() %>% dplyr::select(es, se) %>%
        suppressWarnings() %>%
        dplyr::mutate(es = es*-1)})
  }


#' Calculate Hedges' g based on change data.
#'
#' Calculate Hedges' g based on change data. Only meant to be used as
#' part of \code{\link{calculateEffectSizes}}.
#'
#' @param x data
#' @param ... Change score effect size data. Must be \code{Change_m_trt1}, \code{Change_m_trt2}, \code{Change_SD_trt1},
#' \code{Change_SD_trt2}, \code{Change_N_trt1}, \code{Change_N_trt2} (all \code{numeric}).
#' @usage changeES(x, ...)
#' @importFrom dplyr select
#' @importFrom purrr pmap_dfr
#' @importFrom esc esc_mean_sd
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export changeES
#' @keywords internal

changeES = function(x, ...){
  x %>%
    purrr::pmap_dfr(function(Change_m_trt1, Change_m_trt2, Change_SD_trt1,
                             Change_SD_trt2, Change_N_trt1, Change_N_trt2, ...)
    {esc::esc_mean_sd(Change_m_trt1, Change_SD_trt1, Change_N_trt1,
                      Change_m_trt2, Change_SD_trt2, Change_N_trt2, es.type = "g") %>%
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
                           c("_trt1" = "_trtX", "_trt2" = "_trt1",
                             "_trtX" = "_trt2")) -> colnames(dat)

  dat %>% dplyr::select(colnames(dat.orig)) %>%
    dplyr::mutate(es = es*-1,
           id = paste0(id, "_arm_switched")) %>%
    rbind(., dat.orig) -> dat.expand

  dat.expand = dat.expand[order(dat.expand$id),]
  rownames(dat.expand) = NULL

  return(dat.expand)
}


