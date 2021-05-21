#' Check data format
#'
#' This function checks if a `data.frame` object conforms with the 
#' [Metapsy data standard](https://docs.metapsy.org/data-preparation/format/).
#' 
#'
#' @param data A `data.frame` containint meta-analysis data.
#' @param must.contain \code{character} vector, containing all the variable names 
#' the data set should contain. Defaults correspond with the Metapsy data standard.
#' @param variable.class \code{list}, defining the required class 
#' for some or all variables. If the class differs in `data`, the function will try to 
#' convert the variable to the desired class. Defaults correspond with the 
#' Metapsy data standard.
#'
#' @usage checkDataFormat(data,
#'                 must.contain = c("study", "condition_arm1",
#'                                  "condition_arm2", 
#'                                  "multi_arm1", 
#'                                  "multi_arm2",
#'                                  "outcome_type", "instrument",
#'                                  "time", "time_weeks",
#'                                  "rating", "mean_arm1", "mean_arm2",
#'                                  "sd_arm1", "sd_arm2",
#'                                  "n_arm1", "n_arm2",
#'                                  "event_arm1", "event_arm2",
#'                                  "totaln_arm1", "totaln_arm2"),
#'                 variable.class = list("study" = "character", 
#'                                       "condition_arm1" = "character",
#'                                       "condition_arm2" = "character", 
#'                                       "multi_arm1" = "character", 
#'                                       "multi_arm2" = "character",
#'                                       "outcome_type" = "character", 
#'                                       "instrument" = "character",
#'                                       "time" = "character", 
#'                                       "time_weeks" = "numeric",
#'                                       "rating" = "character", 
#'                                       "mean_arm1" = "numeric", 
#'                                       "mean_arm2" = "numeric",
#'                                       "sd_arm1" = "numeric", 
#'                                       "sd_arm2" = "numeric",
#'                                       "n_arm1" = "numeric", 
#'                                       "n_arm2" = "numeric",
#'                                       "event_arm1" = "numeric", 
#'                                       "event_arm2" = "numeric",
#'                                       "totaln_arm1" = "numeric", 
#'                                       "totaln_arm2" = "numeric"))
#'
#' @return \code{checkDataFormat} returns messages that specify if input variables, 
#' values and classes of the variables are as defined. The output should then be passed to
#' \code{\link[metapsyTools]{checkConflicts}} and then to \code{\link[metapsyTools]{runMetaAnalysis}}.
#' 
#' If default settings are used, `checkDataFormat()` can be used in combination
#' with \code{\link[metapsyTools]{checkConflicts}} to determine if a dataset follows the [Metapsy data standard](https://docs.metapsy.org/data-preparation/format/).
#' Datasets that are formatted using this standard can be directly used in 
#' the [analysis module](https://tools.metapsy.org/articles/metapsytools#the-analysis-module) in `metapsyTools`; 
#' for example the \code{\link[metapsyTools]{runMetaAnalysis}} function.
#' 
#' @details The function checks if:
#'  \itemize{
#'  \item{the data set contains all relevant variables and}
#'  \item{variables have the desired class (if not, it tries to convert).}
#'  }

#'
#' @examples
#' \dontrun{
#' data("depressionPsyCtr")
#'
#' # Example 1: Check with default arguments
#' checkDataFormat(depressionPsyCtr)
#'
#' #Example 2: Check for non-default arguments
#' checkDataFormat(depressionPsyCtr,
#'                 must.contain = c("study", "condition",
#'                                  "primary", "year"),
#'                 variable.class = list(study = "character",
#'                                       no.arms = "numeric"))
#' }
#'
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com}, 
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, 
#' Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{expandMultiarmTrials}}
#' @export checkDataFormat
#' @importFrom crayon green yellow cyan bold
#' @importFrom methods as
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @importFrom stringr str_remove_all

checkDataFormat = function(data,
                           must.contain = c("study", "condition_arm1",
                                            "condition_arm2", 
                                            "multi_arm1", 
                                            "multi_arm2",
                                            "outcome_type", "instrument",
                                            "time", "time_weeks",
                                            "rating", "mean_arm1", "mean_arm2",
                                            "sd_arm1", "sd_arm2",
                                            "n_arm1", "n_arm2",
                                            "event_arm1", "event_arm2",
                                            "totaln_arm1", "totaln_arm2"),
                           variable.class = list("study" = "character", 
                                                 "condition_arm1" = "character",
                                                 "condition_arm2" = "character", 
                                                 "multi_arm1" = "character", 
                                                 "multi_arm2" = "character",
                                                 "outcome_type" = "character", 
                                                 "instrument" = "character",
                                                 "time" = "character", 
                                                 "time_weeks" = "numeric",
                                                 "rating" = "character", 
                                                 "mean_arm1" = "numeric", 
                                                 "mean_arm2" = "numeric",
                                                 "sd_arm1" = "numeric", 
                                                 "sd_arm2" = "numeric",
                                                 "n_arm1" = "numeric", 
                                                 "n_arm2" = "numeric",
                                                 "event_arm1" = "numeric", 
                                                 "event_arm2" = "numeric",
                                                 "totaln_arm1" = "numeric", 
                                                 "totaln_arm2" = "numeric")){

  # Validate input
  validate_data_frame(data, arg_name = "data")
  
  # 1. Check if data set contains all relevant variables.
  if (length(must.contain) > 0){
    missing <- validate_required_columns(data, must.contain, arg_name = "data")
    
    if (length(missing) > 0){
      message("- ", crayon::yellow("[!] "), 
              "Data set does not contain variable(s) ",
              crayon::green(paste(missing, collapse = ", ")), ".")
    } else {
      message("- ", crayon::green("[OK] "), 
              "Data set contains all variables in 'must.contain'.")
    }
  }
  
  # 2. Check if variables are of desired type; if not, try to convert
  if (length(variable.class) > 0){
    for (var_name in names(variable.class)){
      target_class <- variable.class[[var_name]]
      
      # Skip if variable doesn't exist in data
      if (!var_name %in% colnames(data)) {
        next
      }
      
      current_class <- class(data[[var_name]])[1]
      
      if (current_class == target_class){
        message(paste0("- ", crayon::green("[OK] "), 
                       "'", var_name, 
                       "' has desired class ",
                       target_class, "."))
      } else {
        # Try to convert
        data[[var_name]] <- safe_convert_class(
          data[[var_name]], 
          target_class, 
          var_name
        )
        
        # Check if conversion was successful
        new_class <- class(data[[var_name]])[1]
        if (new_class == target_class) {
          message(paste0("- ", crayon::green("[OK] "), "'",
                         var_name, 
                         "' has been converted to class ",
                         target_class, "."))
        } else {
          message(paste0("- ", crayon::yellow("[!] "), "'",
                         var_name, 
                         "' could not be converted from ",
                         current_class, " to ", target_class, "."))
        }
      }
    }
  }
  
  # 3. Check if dataset contains 'subset' or 'exclude', rename if necessary
  rename_result <- rename_reserved_columns(data)
  data <- rename_result$data
  
  if (rename_result$renamed) {
    message(paste0("- ", crayon::green("[OK] "), 
                   "Reserved variable name(s) have been renamed."))
  }
  
  return(data)
}
