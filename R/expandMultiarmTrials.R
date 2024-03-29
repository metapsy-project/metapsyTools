#' Expand multiarm trials
#'
#'
#' This is a function to expand the format of meta-analysis data for its applicability to the \code{calculateEffectSizes} function of the \code{metapsyTools}-package.
#'
#' @param data Meta-analysis data stored as a \code{data.frame}.
#' @param vars.for.id \code{character} vector, containing column names of all variables used to construct unique comparison IDs.
#' @param study.indicator \code{character}, signifying the name of the variable containing the study name.
#' @param multiarm.indicator \code{numeric}, signifying if a row is part of a multiarm study (1) or not (0).
#' @param no.arms.indicator \code{character}, signifying the name of the variable containing the number of arms included in a study (typically 2).
#' @param group.indicator \code{character}, column name of the variable storing the study name.
#' @param condition.specification \code{character}, column name of the variable storing the trial condition name.
#' @param groups.column.indicator \code{character}. If the dataset is in wide format: a character vector with two elements, representing the suffix used to
#' differentiate between the first and second treatment in a comparison.
#' @param group.names \code{list}, storing the name of the value corresponding to the intervention group (\code{"ig"}) and control group (\code{"cg"}).
#' @param data.format \code{character}. Either \code{"long"} or \code{"wide"}, depending on the format of the dataset in
#' \code{data}. \code{NULL} by default, which lets the user define the format after the function has been called.
#'
#' @usage expandMultiarmTrials(data,
#'                      vars.for.id = c("study", "primary",
#'                                      "Outc_measure",
#'                                      "Time", "Time_weeks",
#'                                      "sr_clinician"),
#'                      study.indicator = "study",
#'                      multiarm.indicator = "is.multiarm",
#'                      no.arms.indicator = "no.arms",
#'                      group.indicator = "condition",
#'                      condition.specification = "Cond_spec",
#'                      groups.column.indicator = c("_trt1", "_trt2"),
#'                      group.names = list("ig" = "ig",
#'                                         "cg" = "cg"),
#'                      data.format = NULL)
#'
#' @return \code{expandMultiarmTrials} returns the meta-analysis data set as class \code{data.frame} (if results are saved to a variable). The rows of multiarm studies are expanded so that each intervention group has an unambiguously assigned control group. It also generates the following columns:
#' \itemize{
#' \item{\code{id} a \emph{comparison}-specific ID variable.}
#' \item{\code{study.id} a \emph{study}-specific ID variable.}
#' \item{\code{study} a study-specific variable containing the study name. For multiarm studies, this variable also specifies the active treatment indicated by \code{multiarm.group.indicator} behind the name of the study (e.g. \code{"Hauksson, 2017 -grp"}).}
#' }
#' @details This function expands multiarm studies in a meta-analysis data set, thereby ensuring that each comparison (intervention group vs. control group in \code{condition}) is unique for a specific outcome, and thus has two rows. For this purpose, it duplicates the corresponding row of the control group condition if required. A specific study indicator variable is created that enables further use, e.g. in 3-level models.
#'  For more details see the help vignette: \code{vignette("metapsyTools")}.
#' @examples
#' \dontrun{
#' data("inpatients")
#' expandMultiarmTrials(inpatients)
#' }
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com}, Paula Kuper \email{paula.r.kuper@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{calculateEffectSizes}}
#' @importFrom dplyr filter arrange pull all_of
#' @importFrom stringr str_replace_all str_remove_all
#' @importFrom purrr map_df
#' @importFrom crayon green yellow cyan bold
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export expandMultiarmTrials
#' @keywords internal

expandMultiarmTrials = function(data,
                                vars.for.id = c("study", "primary",
                                                "Outc_measure",
                                                "Time", "Time_weeks",
                                                "sr_clinician"),
                                study.indicator = "study",
                                multiarm.indicator = "is.multiarm",
                                no.arms.indicator = "no.arms",
                                group.indicator = "condition",
                                condition.specification = "Cond_spec",
                                groups.column.indicator = c("_trt1", "_trt2"),
                                group.names = list("ig" = "ig",
                                                   "cg" = "cg"),
                                data.format = NULL){


  warning("This function is deprecated as of version 1.0.0.")
  
  
  # Check for data conflicts flagged in checkConflicts
  if ("checkConflicts" %in% class(data)){
    stop("Data format conflicts were detected. Run 'checkConflicts' to diagnose potential problems.")
  }


  # Check data format
  if (!class(data)[1] %in% c("wide", "long")){
    if (is.null(data.format)){
      # Format switch
      input = readline("Enter: Is the data set in long [l] or wide [w] format? ")
      if (input[1] %in% c("l", "L", "long", "Long")){
        format = "long"
        data.format = "long"
      } else if (input[1] %in% c("w", "W", "wide", "Wide")){
        format = "wide"
        data.format = "wide"
      } else {
        stop("Data format must either be long [l] or wide [w].")
      }
    }

    if (!(data.format[1] %in% c("long", "wide"))){
      # Format switch
      input = readline("Enter: Is the data set in long [l] or wide [w] format? ")
      if (input[1] %in% c("l", "L", "long", "Long")){
        format = "long"
      } else if (input[1] %in% c("w", "W", "wide", "Wide")){
        format = "wide"
      } else {
        stop("Data format must either be long [l] or wide [w].")
      }
    } else {format = data.format}
  } else {
    format = class(data)[1]
  }


  # Wide format
  if (format[1] == "wide"){

    # Generate study ID (e.g. for 3L models)
    data$study.id = data[[study.indicator]]

    # Create comparison id
    data[colnames(data) %in%
            c(vars.for.id,
              paste0(vars.for.id, groups.column.indicator),
              paste0(condition.specification,
                     groups.column.indicator))] %>%
      colnames() -> id.vars

    apply(data, 1,
          function(x){paste(as.character(x[id.vars]), collapse = "_")}) %>%
      stringr::str_replace_all(",", "_") %>% str_remove_all(" ") -> data$id

    # Set class
    class(data) = c("data.frame", "expandMultiarmTrials", "wide")

    message("- ", crayon::green("[OK] "), "Data format is wide, so no multiarm expansion is required.")
    message("- ", crayon::green("[OK] "), "Added unique ID for each comparison.")
    return(data)
  }

  # Long format
  if (format[1] == "long"){

    # Check if all specified variables are in dataset
    allvars = c(vars.for.id, multiarm.indicator,
                no.arms.indicator, group.indicator, study.indicator)
    if (!(sum(allvars %in% colnames(data)) == length(allvars))){
      stop(paste0("variable ", paste(allvars[!allvars %in% colnames(data)],
                                     collapse = ", ")),
           " not found in data set.")
    }

    # Generate study ID (e.g. for 3L models)
    data$study.id = data[[study.indicator]]

    # Create comparison id (multiarm ids appear 3+ times)
    apply(data, 1,
          function(x){paste(as.character(x[vars.for.id]), collapse = "_")}) %>%
      stringr::str_replace_all(",", "_") %>% str_remove_all(" ") -> data$id

    # Filter out multiarm studies
    data[data[names(data) == multiarm.indicator] == 1,] %>%
      dplyr::arrange(id, group.indicator) -> data.multiarm

    # Filter out non-multiarm studies
    data[data[names(data) == multiarm.indicator] == 0,] %>%
      dplyr::arrange(id, group.indicator) -> data.no.multiarm

    # Split multiarm data by id; then expand ID-wise
    data.multiarm %>%
      split.data.frame(.$id) %>%
      purrr::map_df(function(x){
        multiarmExpander(x,
                         condition.specification = condition.specification,
                         group.indicator = group.indicator,
                         group.names = group.names)}) -> data.multiarm

    # Add trt column to data.no.multiarm
    data.no.multiarm$trt = ifelse(data.no.multiarm[[group.indicator]] == group.names[["ig"]],
                                  "trt1", "trt2")

    # Combine
    data = rbind(data.no.multiarm, data.multiarm)

    # Order dataset
    data[order(data$trt),] -> data
    data[order(data$id),] -> data

    # Set class
    class(data) = c("data.frame", "expandMultiarmTrials", "long")

    # Return
    if (sum(table(data$id) != 2) != 0){
      stop("'expandMultiarmTrials' could not create unique comparisons. Please check for residual formatting issues.")
    } else {
      message("- ", crayon::green("[OK] "), "Multiarm trial expansion successful.")
      return(data)
    }
  }



}
