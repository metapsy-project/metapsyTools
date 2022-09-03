#' Check for (potential) data format conflicts
#'
#' This function checks for potential data formatting conflicts that may produce errors or incorrect
#' results when applying the [calculateEffectSizes()] or [runMetaAnalysis()] function later on.
#'
#' @param .data Meta-analysis data stored as a \code{data.frame}, to be checked by the function.
#' @param vars.for.id \code{character} vector, containing column names of all 
#' variables used to construct unique comparison IDs.
#' @param .condition \code{character}. The prefix of the two variables in \code{data} in 
#' which the conditions (e.g. "guided iCBT", "waitlist") of the trial arm comparison are stored.
#' @param .condition.specification \code{character}, name of the column containing the specific condition in each trial arm.
#' For multiarm trials, these conditions \emph{must} be distinct (e.g. \code{"cbt-guided"} and \code{"cbt-unguided"}).
#' @param .groups.column.indicator \code{character}. A character vector with two elements, 
#' representing the suffix used to differentiate between the first and second arm in a comparison.
#'
#' @usage checkConflicts(.data,
#'                vars.for.id = c("study", "outcome_type",
#'                                "instrument", "time",
#'                                "time_weeks",
#'                                "rating"),
#'                .condition = "condition",
#'                .condition.specification = "multi",
#'                .groups.column.indicator = c("_arm1", "_arm2"))
#'
#' @return The type of data returned by \code{checkConflicts} depends on the outcome of the evaluation. 
#' When no problems have been detected, the function simply returns the data set 
#' provided in \code{.data}.
#'
#' When (potential) formatting formatting issues have been detected, the function 
#' throws a message and returns the affected studies/\code{data.frame} columns. 
#' In particular, results are provided within a \code{list} of three objects:
#'
#' \itemize{
#' \item{\code{allConflicts}, a \code{data.frame} containing all affected rows, regardless of conflict type.}
#' \item{\code{idConflicts}, a \code{data.frame} containing rows with ID/number of arms conflicts.}
#' \item{\code{cgConflicts}, a \code{data.frame} containing rows with reference arm conflicts
#' (there must be a unique control/reference group for each comparison).}
#' }
#'
#' The returned list has class \code{checkConflicts}.
#'
#' @examples
#' \dontrun{
#' data("depressionPsyCtr")
#' 
#' # Example 1: Use defaults and simply run checks
#' depressionPsyCtr %>%
#'   checkDataFormat() %>%
#'   checkConflicts() -> res
#' 
#' # Example 2: Overrule defaults; this will produce a conflict
#' depressionPsyCtr %>%
#'   checkDataFormat() %>%
#'   checkConflicts(vars.for.id = "study") -> res
#' 
#' }
#' 
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com}, Paula Kuper \email{paula.r.kuper@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{checkDataFormat}}
#' @importFrom crayon green yellow
#' @importFrom dplyr select filter group_map group_by
#' @importFrom stringr str_replace_all str_remove_all
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export checkConflicts

checkConflicts = function(.data,
                          vars.for.id = c("study", "outcome_type",
                                          "instrument", "time",
                                          "time_weeks",
                                          "rating"),
                          .condition = "condition",
                          .condition.specification = "multi",
                          .groups.column.indicator = c("_arm1", "_arm2")){
  
  data.return = .data
  
  # Get id variables
  .data[colnames(.data) %in%
          c(vars.for.id, 
            paste0(.condition,
                   .groups.column.indicator),
            paste0(.condition.specification,
                   .groups.column.indicator))] %>%
    colnames() -> id.vars
  
  # Create comparison IDs
  apply(.data, 1,
        function(x){paste(as.character(x[id.vars]), 
                          collapse = "_")}) %>%
    stringr::str_replace_all(",", "_") %>% 
    stringr::str_remove_all(" ") -> .data$.id
  
  # Check for conflicts
  .data[.data$.id %in%
          names(table(.data$.id)
                [table(.data$.id) > 1]),] -> conflictIds
  
  # Set multiple CG IDs conflicts (not used)
  .data[NULL,] -> multipleCgIds
  
  list(allConflicts = rbind(conflictIds, multipleCgIds),
       idConflicts = conflictIds,
       cgConflicts = multipleCgIds,
       condition.specification = .condition.specification,
       vars.for.id = vars.for.id) -> returnlist
  
  
  class(returnlist) = c("checkConflicts", "list")
  
  if (nrow(returnlist$idConflicts) == 0){
    message("- ", crayon::green("[OK] "),  "No data format conflicts detected.")
    class(data.return) = c("wide", "data.frame")
    return(data.return)
  } else {
    message("- ", crayon::yellow("[!] "), "Data format conflicts detected!")
    return(returnlist)
  }
  
}

#' Print method for the 'checkConflicts' function
#'
#' This S3 method prints the studies with potential formatting conflicts.
#'
#' @param x A \code{list} object of class \code{checkConflicts}.
#' @param ... Additional arguments (not used).
#'
#'
#' @author Mathias Harrer & David Daniel Ebert
#'
#' @seealso \code{\link{checkConflicts}}
#' @export
#' @method print checkConflicts
#' @importFrom crayon green

print.checkConflicts = function(x, ...){
  if (nrow(x$idConflicts) != 0){
    message(paste0("ID conflicts \n",
                   "- check if variable(s) ", 
                   crayon::green(paste(x$vars.for.id, collapse = ", ")),
                   " create(s) unique assessment point IDs \n",
                   paste0("- check if ", crayon::green(x$condition.specification),
                          " uniquely identifies all trial arms in multiarm trials \n"),
                   "--------------------"))
    cat(paste(unique(x$idConflicts$study), collapse = "\n"))
  }
  
  if (nrow(x$cgConflicts) != 0){
    cat("\n")
    cat("\n")
    message(paste0("Trials with 2+ control groups \n",
                   "- NOTE: As of version 0.2.0, 'metapsyTools' can handle 2+ control groups! \n",
                   "--------------------"))
    cat(paste(unique(x$cgConflicts$study), collapse = "\n"))
  }
  
  if (nrow(x$allConflicts) == 0){
    message("Looks good!")
  }
  
}


