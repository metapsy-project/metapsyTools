#' Calculate effect sizes
#'
#'
#' This is a function to calculate effect sizes of meta-analysis data created by \code{expandMultiarmTrials}.
#'
#' @param data Meta-analysis dataset of \code{expandMultiarmTrials}, created using \code{expandMultiarmTrials}.
#' @param comp.id.indicator Character, column name of the variable storing the comparison ID.
#' @param study.indicator Character, column name of the variable storing the study name.
#' @param group.indicator Character, column name of the variable storing the condition (intervention or control group).
#' @param pivot.vars Character vector, contains the name of the comparison ID variable, the condition variable, and all raw effect size data variables from which effect sizes should be calculated.
#' @param funcs List of functions. These functions will be used to calculate the effect sizes (Hedges' g) based on the raw data.
#'
#'
#' @return \code{calculateEffectSizes} returns the meta-analysis dataset as class \code{data.frame} (if results are saved to a variable). It generates also the following columns:
#' \itemize{
#' \item{\code{es} calculated effect sizes for each comparison}
#' \item{\code{se} calculated standard error of the effect size}
#' \item{\code{study.id} a study specific ID variable}
#' \item{\code{study} a study specific variable containing its name. For multiarm studies it specifies the active treatment behind the name of the study.}
#' }
#' @examples
#' #Example 1: use function in default mode; requires data created by expandMultiarmTrials
#' inpatients %>%
#' expandMultiarmTrials() %>%
#' calculateEffectSizes()
#'
#' #Example 2: further use: pool effect sizes
#' library(meta)
#'
#' inpatients %>%
#' checkDataFormat ()%>%
#' expandMultiarmTrials() %>%
#' calculateEffectSizes() %>%
#' dplyr::filter(!is.na(es) & primary==1) %>%
#' metagen(es, se, studlab=study, comb.fixed=FALSE, data=.)
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com}, Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{expandMultiarmTrials}}
#' @details  By default, g is calculated from:
#' \itemize{
#'  \item {the mean, SD and N}
#'  \item {binary outcome data (i.e. response counts) and}
#'  \item {change scores}
#' }
#' Other functions can be added to the list. Results of the function must result in a data.frame
#' with the same number of rows as in data, and two columns: one for the calculated g value (named es) and its standard error (named se).
#' In rows for which no fitting raw data was supplied, the resulting data.frame should contain NA.
#'
#' For more details see the help vignette: \code{vignette("help", package = "metapsyTools")}
#'
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom dplyr all_of select filter mutate arrange
#' @importFrom purrr pmap_dfr
#' @importFrom esc esc_mean_sd
#' @import magrittr
#' @import dplyr
#' @export

calculateEffectSizes = function(data,
                                comp.id.indicator = "id",
                                study.indicator = "study",
                                group.indicator = "condition",
                                pivot.vars = c("id", "condition",
                                               "Post_M", "Post_SD", "Post_N",
                                               "Rand_N", "Improved_N",
                                               "Change_m", "Change_SD",
                                               "Change_N"),
                                funcs = list(
                                  mean.sd = function(x, ...){
                                    x %>%
                                      purrr::pmap_dfr(function(Post_M_ig, Post_M_cg, Post_SD_ig,
                                                               Post_SD_cg, Post_N_ig, Post_N_cg, ...)
                                      {esc::esc_mean_sd(Post_M_ig, Post_SD_ig, Post_N_ig,
                                                        Post_M_cg, Post_SD_cg, Post_N_cg, es.type = "g") %>%
                                          as.data.frame() %>% dplyr::select(es, se) %>%
                                          suppressWarnings()})
                                  },
                                  binary = function(x, ...){
                                    x %>%
                                      purrr::pmap_dfr(function(Improved_N_ig, Improved_N_cg, Rand_N_ig, Rand_N_cg, ...)
                                      {esc::esc_2x2(Improved_N_ig,
                                                    Rand_N_ig - Improved_N_ig,
                                                    Improved_N_cg,
                                                    Rand_N_cg - Improved_N_cg,
                                                    es.type = "g") %>%
                                          as.data.frame() %>% dplyr::select(es, se) %>%
                                          suppressWarnings() %>%
                                          mutate(es = es*-1)})
                                  },
                                  change = function(x, ...){
                                    x %>%
                                      purrr::pmap_dfr(function(Change_m_ig, Change_m_cg, Change_SD_ig,
                                                               Change_SD_cg, Change_N_ig, Change_N_cg, ...)
                                      {esc::esc_mean_sd(Change_m_ig, Change_SD_ig, Change_N_ig,
                                                        Change_m_cg, Change_SD_cg, Change_N_cg, es.type = "g") %>%
                                          as.data.frame() %>% dplyr::select(es, se) %>%
                                          suppressWarnings()})})
                                ){

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
    dplyr::select(dplyr::all_of(pivot.vars)) %>%
    tidyr::pivot_wider(names_from = dplyr::all_of(group.indicator),
                values_from = c(-id, -dplyr::all_of(group.indicator))) %>%
    {.[,colSums(is.na(.))<nrow(.)]} -> data.wide

  # Apply funcs
  message("calculating effect sizes...")
  es.res = list()
  for (i in 1:length(funcs)){
    es.res[[i]] = funcs[[i]](data.wide)
  }
  es.res = do.call(cbind, es.res)
  message("SUCCESS")

  # Now, bind all calculated ES together,
  # then bind together with wide dataset
  es.res %>%
    apply(., 1, function(x) x[!is.na(x)]) %>% t() %>%
    cbind(data.wide, .) -> data.wide.es

  # Now, transform the wide format data set with calculated ES
  # back to long and merge back with original version
  data.wide.es %>%
    dplyr::select(id, es, se) %>%
    dplyr::mutate(trt1 = "t.1", trt2 = "t.2") %>%
    tidyr::pivot_longer(-id,
                 names_to = c(".value"),
                 names_pattern = "(..)") %>%
    dplyr::arrange(id) %>%
    dplyr::select(2:3) %>%
    cbind(data %>% arrange(id) %>% select(-id), .) -> dat.final

  # Return
  return(dat.final)
}









