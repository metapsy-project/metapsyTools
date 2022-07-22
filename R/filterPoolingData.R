#' Filter data to be pooled for meta-analysis
#'
#' This convenience function allows to create a filtered data set, which is then ready
#' to be used for meta-analytic pooling.
#'
#' @usage filterPoolingData(.data, ...,
#'                   .filter.missing.rows = FALSE,
#'                   .es.column = es)
#'
#' @param .data A \code{data.frame} containing the calculated effect sizes, as created by the \code{\link{calculateEffectSizes}} function.
#' @param ... <\link[dplyr]{dplyr_data_masking}>. A number of filtering statements (using variables in \code{.data})
#' that return a logical value. To apply multiple filters, simply separate them using a comma. "OR" statements can be
#' provided using the \code{|} operator. Multiple filter statements separated using commas will be combined using the AND (\code{&}) operator. See "Details".
#' @param .filter.missing.rows \code{logical}. Should rows with no effect sizes be filtered out? Default is \code{FALSE}.
#' @param .es.column Name of the column in \code{.data} to be used for filtering out rows with no effect sizes. Default is \code{es}.
#'
#' @return \code{filterPoolingData} returns the filtered data set as class \code{data.frame}.
#' The filtered data set should then be ready for meta-analytic pooling, for example 
#' using [runMetaAnalysis()].
#'
#' @examples
#' \dontrun{
#'
#' # Example 1: calculate effect sizes and then use multiple AND filters.
#' data("depressionPsyCtr")
#' depressionPsyCtr %>%
#'   calculateEffectSizes() %>%
#'   filterPoolingData(time == "post", instrument == "hdrs")
#'
#' # Example 2: use OR filter
#' data("depressionPsyCtr")
#' depressionPsyCtr %>%
#'   calculateEffectSizes() %>%
#'   filterPoolingData(time == "post" | instrument == "hdrs")
#'
#' # Example 3: use %in% operator
#' data("depressionPsyCtr")
#' depressionPsyCtr %>%
#'   calculateEffectSizes() %>%
#'   filterPoolingData(instrument %in% c("hdrs", "phq-9"))
#'
#' # Example 4: Search for studies using "fuzzy-ish" matching
#' data("depressionPsyCtr")
#' depressionPsyCtr %>%
#'   calculateEffectSizes() %>%
#'   filterPoolingData(Detect(instrument, "bdi"))
#' }
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com}, 
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{filterPriorityRule}}
#'
#' @details The \code{filterPoolingData} function allows to apply several filters to your meta-analysis data set all at once.
#' When used in a pipe (\code{\%>\%}), you only need to supply several filtering statements separated by commas, using the same column names as they appear
#' in the data set (e.g. \code{primary == 1, meanage == 58}). The filtering statements are then connected using "AND" (\code{&}).
#'
#' If you want to apply an "OR" filter, simply use \code{|} instead of a comma (e.g. \code{type == "cbt" | format == 6}). To select all rows
#' that contain one of several values in a variable, use \code{\%in\%}; e.g. \code{study \%in\% c("Bailey, 2017", "Barth 2005")}.
#'
#' The \code{\link{Detect}} function can be used within the function call to search for variable elements which \emph{contain} one or several
#' selected words (separated by \code{|}). To include all rows which contain the word "cbt" or "wl" or "cau" in the "Cond_spec_trt2" variable,
#' we can use \code{Detect(Cond_spect_trt2, "cbt|wl|cau")}. This will also filter out elements like \code{"cbt (online)"}, because "cbt" is
#' included.
#'
#' For more details see the [Get Started](https://tools.metapsy.org/articles/metapsytools) vignette.
#'
#' @import dplyr
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export filterPoolingData

filterPoolingData = function(.data, ...,
                             .filter.missing.rows = FALSE,
                             .es.column = es){
  
  # Check class
  if (!inherits(.data, "data.frame") & 
      !inherits(.data, "metapsyDatabase")){
    stop("'.data' must be a data.frame or 'metapsyDatabase' R6 object.")
  }
  
  if (inherits(.data, "metapsyDatabase")){
    .data = .data[["data"]]
  }

  es = enquo(.es.column)

  if (.filter.missing.rows == TRUE){

    return(dplyr::filter(.data, ..., !is.na(!!es)))

  } else {

    return(dplyr::filter(.data, ...))

  }

}





