#' Calculate effect sizes
#'
#' This is a function to calculate effect sizes of meta-analysis data created by \code{expandMultiarmTrials}.
#'
#' @usage calculateEffectSizes(data,
#'        comp.id.indicator = "id",
#'        trt.indicator = "trt",
#'        funcs = list(meanSD = meanSD,
#'                     binaryES = binaryES,
#'                     changeES = changeES),
#'                     change.sign = NULL)
#'
#'
#' @param data Meta-analysis data set of class \code{expandMultiarmTrials}, created using \code{expandMultiarmTrials}.
#' @param comp.id.indicator \code{character}, column name of the variable storing the comparison ID; typically created by \code{expandMultiarmTrials}.
#' @param trt.indicator \code{character}, column name of the variable storing the treatment indicator (treatment 1 or 2); typically created by \code{expandMultiarmTrials}.
#' @param funcs \code{list} of functions. These functions will be used to calculate the effect sizes (Hedges' \emph{g}) based on the raw data (see Details).
#' @param change.sign \code{character}. Name of a \code{logical} column in \code{data}, encoding if the
#' sign of a calculated effect size should be reversed (\code{TRUE}) or not (\code{FALSE}). Set to \code{NULL} (default) if
#' no changes should be made.
#'
#' @return \code{calculateEffectSizes} returns the meta-analysis data set as class \code{data.frame} in wide format (if results are saved to a variable). It also generates the following columns, wich are added to the data:
#' \itemize{
#' \item{\code{es} calculated effect sizes for each comparison.}
#' \item{\code{se} calculated standard error of the effect size.}
#' \item{\code{study.id} a study-specific ID variable.}
#' \item{\code{study} a study-specific variable containing its name. For multiarm studies, this variable specifies the active treatment arm used to calculate the effect.}
#' }
#' @examples
#' \dontrun{
#' #Example 1: use function in default mode; requires data created by expandMultiarmTrials
#' data("inpatients")
#' inpatients %>%
#'     expandMultiarmTrials() %>%
#'     calculateEffectSizes()
#'
#' #Example 2: further use to pool effect sizes
#' library(meta)
#' inpatients %>%
#'   checkDataFormat() %>%
#'   expandMultiarmTrials() %>%
#'   calculateEffectSizes() %>%
#'   filterPoolingData(primary==1) %>%
#'   metagen(es, se, studlab=study, comb.fixed=FALSE, data=.)
#'
#'
#' # Example 3: use for 3-level model
#' library(metafor)
#' library(dplyr)
#' inpatients %>%
#'   checkDataFormat ()%>%
#'   expandMultiarmTrials() %>%
#'   calculateEffectSizes() %>%
#'   dplyr::mutate(es.id = 1:nrow(.)) %>%
#'   metafor::rma.mv(es, se^2, data = ., slab = study.id,
#'                   random = ~ 1 | study.id/es.id, test = "t")
#' }
#'
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com}, Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{expandMultiarmTrials}}
#' @details  By default, the small-sample bias corrected standardized mean difference (Hedges' \emph{g}) is calculated from:
#' \itemize{
#'  \item {the mean, SD and N}
#'  \item {binary outcome data (i.e. response counts) and}
#'  \item {change scores}
#' }
#' Other functions can be added to the list provided to \code{funcs}. However, results of the function must result in a \code{data.frame}
#' with the same number of rows as in \code{data}, and two columns: one for the calculated \emph{g} value (named \code{es}) and its standard error (named \code{se}).
#' In rows for which no fitting raw data was supplied, the resulting \code{data.frame} should contain \code{NA}.
#'
#' For more details see the help vignette: \code{vignette("metapsyTools")}
#'
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom dplyr all_of select filter mutate arrange
#' @importFrom purrr pmap_dfr map
#' @importFrom esc esc_mean_sd
#' @importFrom stringr str_remove str_replace
#' @import magrittr
#' @import dplyr
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export calculateEffectSizes

calculateEffectSizes = function(data,
                                comp.id.indicator = "id",
                                trt.indicator = "trt",
                                funcs = list(meanSD = meanSD,
                                             binaryES = binaryES,
                                             changeES = changeES),
                                change.sign = NULL){

  # check class
  if (class(data)[2] != "expandMultiarmTrials"){
    message(paste("class of 'data' is not 'expandMultiarmTrials'.",
                  "Sure that the data has the right format?"))
  }

  # Convert to data.frame (conventional)
  data = data.frame(data)

  # Define id
  data$id = data[[comp.id.indicator]]


  # Convert to wide
  data %>%
    tidyr::pivot_wider(names_from = dplyr::all_of(trt.indicator),
                       values_from = c(-id, -dplyr::all_of(trt.indicator), -study)) -> data.wide


  # Apply funcs
  message("- Calculating effect sizes...")
  es.res = list()
  for (i in 1:length(funcs)){
    es.res[[i]] = try({funcs[[i]](data.wide)}, silent = TRUE)
  }

  # Check if all functions could be applied
  es.res %>% purrr::map(~class(.)) %>%
    unlist() %>% {. == "try-error"} -> error.mask

  funcs2 = funcs[!error.mask]
  if (length(funcs) == 0){
    stop("None of the supplied functions could be applied correctly.")
  }

  # Recalculate if necessary
  if (sum(error.mask) != 0){
    es.res = list()
    for (i in 1:length(funcs2)){
      es.res[[i]] = try({funcs2[[i]](data.wide)}, silent = TRUE)}
  }

  es.res = do.call(cbind, es.res)

  if (sum(error.mask) > 0){
    message("- [!] Function(s) ", paste(names(funcs)[error.mask], collapse = ", "),
            " not applied. Check for potential data/function problems.")
    message("- [!] All other effect size calculation functions were applied successfully.")
  } else {
    message("- [OK] Effect sizes calculated successfully")
  }

  # Now, bind all calculated ES together,
  # then bind together with wide dataset
  es.res %>%
    apply(., 1, function(x){
      if (length(x[!is.na(x)]) == 0){return(c(NA, NA))} else {
        return(x[!is.na(x)])}}) %>% t() %>%
    cbind(data.wide, .) -> data.wide.es

  # Remove duplicate columns
  data.wide.es[!duplicated(t(data.wide.es))] -> data.wide.es
  for (i in 1:ncol(data.wide.es)){
    if (grepl("_trt1", colnames(data.wide.es)[i]) &
              !(stringr::str_replace(colnames(data.wide.es)[i], "_trt1", "_trt2") %in%
              colnames(data.wide.es))){
      colnames(data.wide.es)[i] = stringr::str_remove(colnames(data.wide.es)[i], "_trt1")
    }
  }
  data.wide.es -> dat.final

  # Change sign
  if (!is.null(change.sign)){
    change.mask = ifelse(dat.final[[change.sign]], -1, 1)
    dat.final$es = dat.final$es * change.mask
  }

  # Return
  return(dat.final)

}









