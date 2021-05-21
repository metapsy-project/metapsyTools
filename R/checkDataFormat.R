#' Check data format
#'
#'
#' This is a function to check the format of meta-analysis data with class \code{data.frame} for its further applicability to the \code{expandMultiarmTrials} function.
#'
#' @param data Meta-analysis data stored as a data.frame.
#' @param must.contain \code{character} vector, containing all the variable names the data set should contain.
#' @param variable.contains \code{list}, defining which values should only be contained in a variable.
#' @param variable.class \code{list}, defining the required data class for some or all variables. If the class differs, the function will try to convert the variable to the desired class.
#'
#' @usage checkDataFormat(data,
#'            must.contain = c("study", "condition",
#'                             "Cond_spec", "is.multiarm",
#'                             "no.arms", "multiple_arms",
#'                             "Outc_measure", "Time", "primary",
#'                             "Time_weeks"),
#'            variable.contains = list("condition" = c("ig", "cg"),
#'                                     "is.multiarm" = c(0,1)),
#'            variable.class = list("study" = "character",
#'                                  "condition" = "character",
#'                                  "Post_M" = "numeric",
#'                                  "Post_SD" = "numeric",
#'                                  "Post_N" = "numeric",
#'                                  "Rand_N" = "numeric",
#'                                  "Improved_N" = "numeric",
#'                                  "Change_m" = "numeric",
#'                                  "Change_SD" = "numeric",
#'                                  "Change_N" = "numeric"))
#'
#' @return \code{checkDataFormat} returns messages that specify if input variables, values and classes of the variables are as defined.
#' @details The function ckecks if:
#'  \itemize{
#'  \item{the data set contains all relevant variables}
#'  \item{variables contain desired values only and}
#'  \item{variables are of desired type (if not, it tries to convert them).}
#'  }
#' @examples
#' data("inpatients")
#'
#' # Example 1: Check data with default arguments
#' checkDataFormat(inpatients)
#'
#' #Example 2: Check for specified variables and corresponding values and classes
#' checkDataFormat(inpatients,
#'                 must.contain = c("study", "condition",
#'                                  "primary", "year"),
#'                 variable.contains = list("condition" = c("ig","cg")),
#'                 variable.class = list(study = "character",
#'                                       no.arms = "numeric"))
#'
#' #Example 3: Convert variable class as predefined
#' \dontrun{
#' checkDataFormat(inpatients,
#'                 must.contain = c("primary", "study"),
#'                 variable.class = list(primary = "integer"))
#' }
#'
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com}, Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{expandMultiarmTrials}}
#' @export checkDataFormat
#' @importFrom methods as

checkDataFormat = function(data,
                           must.contain = c("study", "condition",
                                            "Cond_spec", "is.multiarm",
                                            "no.arms", "multiple_arms",
                                            "Outc_measure", "Time", "primary",
                                            "Time_weeks"),
                           variable.contains = list("condition" = c("ig", "cg"),
                                                    "is.multiarm" = c(0,1)),
                           variable.class = list("study" = "character",
                                                 "condition" = "character",
                                                 "Post_M" = "numeric",
                                                 "Post_SD" = "numeric",
                                                 "Post_N" = "numeric",
                                                 "Rand_N" = "numeric",
                                                 "Improved_N" = "numeric",
                                                 "Change_m" = "numeric",
                                                 "Change_SD" = "numeric",
                                                 "Change_N" = "numeric")){

  # 1. Check if data set contains all relevant variables.
  if (length(must.contain) > 0){
    if (sum(colnames(data) %in% must.contain) < length(must.contain)){
      must.contain[which(!must.contain %in% colnames(data))] -> miss
      message(paste0("! data set does not contain variable(s) ",
                     paste(miss, collapse = ", "), "."))
    } else {
      message("- data set contains all variables in 'must.contain'.")
    }
  }

  # 2. Check if variables only contain desired values.
  if (length(variable.contains) > 0){
    issue = FALSE
    issueList = list()
    for (i in 1:length(variable.contains)){
      x = unique(data[[names(variable.contains)[i]]])
      if (length(x) == length(variable.contains[[i]]) &
          sum(x %in% variable.contains[[i]]) == length(x)){
        issueList[[i]] = 0
      } else {
        issue = TRUE
        issueList[[i]] = 1
      }
    }

    if (issue == TRUE){
      message(paste0("! ",
                     paste(names(variable.contains)[issueList == 1],
                           collapse = ", "),
                     " not (only) contains the values specified in 'variable.contains'."))
    } else {
      message("- variables contain only the values specified in 'variable.contains'.")
    }
  }

  # 3. Check if variables are of desired type; if not, try to convert
  if (length(variable.class) > 0){
    for (i in 1:length(variable.class)){
      if (class(data[[names(variable.class)[i]]]) == variable.class[[i]]){
        message(paste0("- '", names(variable.class)[i], "' has desired class ",
                       variable.class[[i]], "."))
      } else {
        data[[names(variable.class)[i]]] = as(data[[names(variable.class)[i]]],
                                              variable.class[[i]])
        message(paste0("- '", names(variable.class)[i], "' has been converted to class ",
                       variable.class[[i]], "."))
      }
    }
  }
  return(data)
}




