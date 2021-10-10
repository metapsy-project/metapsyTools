#' Add information that varies between trial arms as extra columns to your meta-analysis dataset
#'
#' Creates two additional columns for each selected variable in which information is stored
#' separately for the intervention and control group. This is typically useful when trial-level variables
#' (i.e. variables that differ between trial arms) are to be included in the final meta-analysis dataset.
#'
#' @usage addTrialArmInfo(.data, ..., .arm.variable,
#'             .group.indicator = "condition",
#'             .name.intervention.group = "ig",
#'             .name.control.group = "cg",
#'             .vars.for.id = c("study", "primary",
#'                              "Outc_measure",
#'                              "Time", "Time_weeks"))
#'
#' @param .data A \code{data.frame} containing unique intervention-control group comparisons, as created by the \code{\link{expandMultiarmTrials}} function.
#' @param ... <\link[dplyr]{dplyr_data_masking}>. The name of several columns (included in \code{.data})
#' that are trial-level variables to be added as columns to \code{.data}. To add multiple variables, simply separate them using a comma.
#' @param .group.indicator \code{character}. Name of the column in \code{.data} which encodes the intervention/control group rows.
#' @param .name.intervention.group \code{character}. Name used in the \code{.group.indicator} variable to identify the intervention group rows.
#' @param .name.control.group \code{character}. Name used in the \code{.group.indicator} variable to identify the control group rows.
#' @param .vars.for.id \code{character} vector, containing column names of all variables used to construct unique comparison IDs.
#'
#' @return \code{addTrialArmInfo} returns a dataset as class \code{data.frame}. This dataset
#' contains all the information previously stored in \code{.data}, plus two columns for each selected
#' trial arm variable (one for the intervention and one for the control group).
#'
#' @examples
#' \dontrun{
#'
#' # Example 1: calculate effect sizes
#' # then add "Post_N" as trial arm variable
#' data("inpatients")
#' inpatients %>%
#'   checkDataFormat() %>%
#'   expandMultiarmTrials() %>%
#'   calculateEffectSizes() %>%
#'   addTrialArmInfo(Post_N) %>%
#'   filterPoolingData(primary == 1)
#'
#' # Example 2: add several trial arm variables simultaneously
#' inpatients %>%
#'   checkDataFormat() %>%
#'   expandMultiarmTrials() %>%
#'   calculateEffectSizes() %>%
#'   addTrialArmInfo(Post_N, Rand_N, Cond_spec) %>%
#'   filterPoolingData(primary == 1)
#' }
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com}, Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{expandMultiarmTrials}}
#'
#' @details Before running the meta-analysis, it is necessary to select only the rows containing calculated
#' effect sizes ('\code{es}'). This results in an information loss when data differs between trial
#' arms within one study (e.g. the sample size \emph{n} is often not identical in both arms of
#' a study); only the row of the "active"/intervention arm is selected, and the information of
#' the control group arm is discarded.
#'
#' \code{addTrialArmInfo} is a convenience function which allows to avoid this information loss
#' by adding trial-specific information as extra columns in the dataset. Two columns are created
#' for each feature: one containing the value of the intervention arm, and another containing the
#' information in the control arm.
#'
#' The function is only applicable to datasets with expanded multiarm trial; that is, the
#' output of \code{\link{expandMultiarmTrials}} (or \code{\link{expandMultiarmTrials}},
#' followed by \code{\link{calculateEffectSizes}}).
#'
#' For more details see the help vignette: \code{vignette("metapsyTools")}.
#'
#' @importFrom stringr str_replace_all
#' @importFrom dplyr filter select arrange all_of
#' @importFrom purrr map2
#' @export addTrialArmInfo

addTrialArmInfo = function(.data, ..., .arm.variable,
                           .group.indicator = "condition",
                           .name.intervention.group = "ig",
                           .name.control.group = "cg",
                           .vars.for.id = c("study", "primary",
                                           "Outc_measure",
                                           "Time", "Time_weeks")){

  if (class(.data) != "data.frame"){
    stop("object provided in '.data' must be of class 'data.frame'.")
  }

  # Re-Create unique comparison IDs
  apply(.data, 1,
        function(x){paste(as.character(x[.vars.for.id]), collapse = "_")}) %>%
    stringr::str_replace_all(",", "_") %>% str_remove_all(" ") -> .data$id

  # Check if comparisons are unique
  if (!(length(unique(.data$id)) == nrow(.data)/2)){
    stop("data.frame must consist of unique intervention vs. control group comparisons. Did you use expandMultiarmTrials()?")
  }

  # Extract all trial arm variables to be added to data set
  .data %>% dplyr::select(...) %>% colnames() -> trial.vars

  # Loop over all selected variables
  purrr::map2(list(.data), as.list(trial.vars), function(x,y){

    # Extract values for IG
    x %>%
      dplyr::filter(condition == .name.intervention.group) %>%
      dplyr::select(id, dplyr::all_of(y)) -> ig

    # Extract values
    x %>%
      dplyr::filter(condition == .name.control.group) %>%
      dplyr::select(id, dplyr::all_of(y)) -> cg

    # Combine
    merge(ig, cg, by = "id", suffixes = c(paste0("_", .name.intervention.group),
                                          paste0("_", .name.control.group))) %>%
      {rbind(.,.)} %>%
      dplyr::arrange(id)

  }) %>% do.call(cbind, .) %>%
    {cbind(id = .[,1], dplyr::select(., -id))} -> res

  message("- [OK] Attaching trial arm variable(s) ",
          paste(trial.vars, collapse = ", "))

  # Merge back with original .data
  cbind(dplyr::arrange(.data, id), res) %>%
    dplyr::select(-id) %>%
    return()

}


