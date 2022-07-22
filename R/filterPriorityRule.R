#' Filter data based on a priority rule
#'
#' This function filters rows of a dataset based on a priority rule for specific variables
#' defined by the user.
#'
#' @usage filterPriorityRule(.data, ..., .study.indicator = "study")
#'
#' @param .data A \code{data.frame} containing the calculated effect sizes, as created by the \code{\link{calculateEffectSizes}} function.
#' @param ... <\link[dplyr]{dplyr_data_masking}>. A number of prioritized filtering rules for variables.
#' Should follow the form \code{variable = c("prio1", "prio2", ...)}. To apply multiple priority filters,
#' simply separate them using a comma. For each study, rows are then selected based on the specified hierarchy for
#' a variable. The priorities are provided as a concatenated vector, representing the variable levels. The level
#' to appear first in this vector has the highest priority, the second one the second-largest priority, and so on.
#' If a study contains none of the variable levels specified in the function call, the study is omitted entirely.
#' @param .study.indicator \code{character}. Name of the variable in which the study IDs are stored.
#'
#' @return \code{filterPriorityRule} returns the filtered data set as class \code{data.frame}.
#' The filtered data set should then be ready for meta-analytic pooling, for example using \link[meta]{metagen}.
#' Further filters can be applied using \code{\link{filterPoolingData}}.
#'
#' @examples
#' \dontrun{
#' # Load data and calculate effect size
#' data("depressionPsyCtr")
#' depressionPsyCtr %>%
#'   checkDataFormat() %>%
#'   checkConflicts() %>%
#'   calculateEffectSizes() -> data
#'
#' # Filter using four priority rules
#' filterPriorityRule(data,
#'                    condition_arm1 = c("cbt", "pst"),
#'                    condition_arm2 = c("cau", "wl", "cbt"),
#'                    instrument = c("cesd", "phq-9", "scl", "hdrs"),
#'                    time = c("post", "fu")) -> res
#' }
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com}, 
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, 
#' Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{filterPoolingData}}
#'
#' For more details see the [Get Started](https://tools.metapsy.org/articles/metapsytools) vignette.
#'
#' @import dplyr
#' @importFrom purrr map_df
#' @importFrom rlang eval_tidy
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export filterPriorityRule


filterPriorityRule = function(.data, ..., 
                              .study.indicator = "study"){
  
  # Check class
  if (!inherits(.data, "data.frame") & 
      !inherits(.data, "metapsyDatabase")){
    stop("'.data' must be a data.frame or 'metapsyDatabase' R6 object.")
  }
  
  if (inherits(.data, "metapsyDatabase")){
    .data = .data[["data"]]
  }

  rules = dplyr::enquos(...)
  vars = names(rules)
  data = .data
  study.indicator = .study.indicator

  for (i in 1:length(rules)){

    data %>%
      split(.[[study.indicator]]) %>%
      purrr::map_df(function(x){

        data.frame(factor = rlang::eval_tidy(rules[[i]]),
                   weight = length(rlang::eval_tidy(rules[[i]])):1) -> weight

        data.frame(factor = x[[vars[i]]],
                   dat = as.numeric(x[[vars[i]]] %in%
                                      rlang::eval_tidy(rules[[i]]))) -> dat

        merge(dat, weight, by = "factor", all.x = TRUE) -> tab
        tab[is.na(tab)] = 0

        if (sum(with(tab, {dat*weight})) > 0){
          with(tab, {dat*weight == max(dat*weight)}) -> mask
          x[x[[vars[i]]] %in% unique(tab[mask, "factor"]),]
        } else {
          x = x[NULL,]
        }

      }) -> data
  }

  return(data)
}



