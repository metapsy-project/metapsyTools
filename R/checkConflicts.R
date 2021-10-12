#' Check for (potential) data format conflicts
#'
#' This function checks for potential data formatting conflicts that may produce errors or incorrect
#' results when applying the \code{\link{expandMultiarmTrials}} function later on.
#'
#' @param .data Meta-analysis data stored as a \code{data.frame}, to be checked by the function.
#' @param vars.for.id \code{character} vector, containing column names of all variables used to construct unique comparison IDs.
#' @param .no.arms \code{character}, signifying the name of the variable containing the number of arms included in a study (typically 2).
#' @param .condition.specification \code{character}, name of the column containing the specific condition in each trial arm.
#' For multiarm trials, these conditions \emph{must} be distinct (e.g. \code{"cbt-guided"} and \code{"cbt-unguided"}).
#' @param .group.indicator \code{character}, variable encoding if a row represents an intervention or control/reference arm.
#' @param .group.names \code{list}, storing the name of the value in \code{.group.indicator} corresponding to the intervention group (\code{"ig"}) and control group (\code{"cg"}).
#'
#' @usage checkConflicts(.data,
#'                vars.for.id = c("study", "primary",
#'                                "Outc_measure",
#'                                "Time", "Time_weeks",
#'                                "sr_clinician"),
#'                .no.arms = "no.arms",
#'                .condition.specification = "Cond_spec",
#'                .group.indicator = "condition",
#'                .group.names = list("ig" = "ig",
#'                                    "cg" = "cg"))
#'
#' @return The type of data returned by \code{checkConflicts} depends on the outcome of the evaluation. When no problems
#' have been detected, the function simply returns the dataset provided in \code{.data}.
#'
#' When (potential) formatting formatting issues have been detected, the function throws a message and returns the affected
#' studies/\code{data.frame} columns. In particular, results are provided within a \code{list} of three objects:
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
#' data("inpatients")
#' data("psyCtrSubset")
#'
#' # Example 1: Use defaults and simply run checks
#' psyCtrSubset %>%
#'   checkDataFormat() %>%
#'   checkConflicts() -> result
#'
#' # Example 2: Use defaults, run check and directly
#' # proceed expanding multiarm trials/calculating effects
#' inpatients %>%
#'   checkDataFormat() %>%
#'   checkConflicts() %>%
#'   expandMultiarmTrials() %>%
#'   calculateEffectSizes() -> result
#'
#' # Example 3: Define custom IDs (example generates error)
#' vars.for.id = c("study", "primary",
#'                 "Time", "Time_weeks")
#' psyCtrSubset %>%
#'   checkDataFormat() %>%
#'   checkConflicts(vars.for.id) %>%
#'   expandMultiarmTrials(vars.for.id) -> res
#'  }
#'
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com}, Paula Kuper \email{paula.r.kuper@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{expandMultiarmTrials}}
#' @importFrom dplyr select filter group_map group_by
#' @importFrom stringr str_replace_all str_remove_all
#' @importFrom magrittr set_names
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export checkConflicts



checkConflicts = function(.data,
                          vars.for.id = c("study", "primary",
                                          "Outc_measure",
                                          "Time", "Time_weeks",
                                          "sr_clinician"),
                          .no.arms = "no.arms",
                          .condition.specification = "Cond_spec",
                          .group.indicator = "condition",
                          .group.names = list("ig" = "ig",
                                              "cg" = "cg")){

  data.return = .data

  # Extract no.arms variable and id.vars
  .data$no.arms = as.numeric(.data[[.no.arms]])
  id.vars = vars.for.id

  apply(.data, 1,
        function(x){paste(as.character(x[id.vars]), collapse = "_")}) %>%
    stringr::str_replace_all(",", "_") %>% stringr::str_remove_all(" ") -> .data$id

  # Check for number of ID - no.arms conflict
  table(.data$id, .data$no.arms) %>% as.matrix() -> M
  b = as.numeric(colnames(M))
  M[!(rowSums(t(t(M) %/% b)) == 1),] %>% rownames() -> conflictIds

  # Check for unspecified condition specifications
  .data$condition.expanded = paste0(.data[[.group.indicator]], "_",
                                    .data[[.condition.specification]])
  .data %>% dplyr::group_by(id) %>%
    dplyr::group_map(~ mean(.x$no.arms) - length(unique(.x$condition.expanded))) %>%
    do.call(c, .) %>% magrittr::set_names(unique(.data$id)) %>%
    {.[. != 0]} %>% names() %>% union(., conflictIds) -> conflictIds

  .data$condition.expanded = NULL

  # Check for double control groups
  table(.data$id, .data[[.group.indicator]] == .group.names[["cg"]])[,2] %>%
    {.[.>1]} %>% names() -> multipleCgIds

  list(allConflicts = .data %>%
         dplyr::filter(id %in% union(conflictIds, multipleCgIds)),
       idConflicts = .data %>%
         dplyr::filter(id %in% conflictIds),
       cgConflicts = .data %>%
         dplyr::filter(id %in% multipleCgIds),
       condition.specification = .condition.specification) -> returnlist


  class(returnlist) = c("checkConflicts", "list")

  if (nrow(returnlist$idConflicts) == 0){
    message("- [OK] No data format conflicts detected")
    return(data.return)
  } else {
    message("- [!] Data format conflicts detected!")
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

print.checkConflicts = function(x, ...){
  if (nrow(x$idConflicts) != 0){
    message(paste0("ID conflicts \n",
                   "- check if the specified number of arms is correct \n",
                   "- check if selected variables really create unique assessment point IDs \n",
                   paste0("- check if '", x$condition.specification,
                          "' uniquely identifies all trial arms in multiarm trials \n"),
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


