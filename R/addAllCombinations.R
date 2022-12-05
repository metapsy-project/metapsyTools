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
#' @importFrom dplyr select ends_with filter group_map group_by
#' @importFrom stringr str_replace_all str_remove_all
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export addAllCombinations



