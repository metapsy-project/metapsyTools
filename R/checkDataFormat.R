# data: meta-analysis data stored as a data.frame
# must.contain: character vector, containing all the variable names the
# data set should contain.
# variable.contains: list, defining which values should only be contained in a variable.
# variable.class: list, defining the required data class for some or all variables. If the
# class differs, the function will try to convert the variable to the desired class.

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
      message(paste0("✗ data set does not contain variable(s) ",
                     paste(miss, collapse = ", "), "."))
    } else {
      message("✓ data set contains all variables in 'must.contain'.")
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
      message(paste0("✗ ",
                     paste(names(variable.contains)[issueList == 1],
                           collapse = ", "),
                     " not (only) contains the values specified in 'variable.contains'."))
    } else {
      message("✓ variables contain only the values specified in 'variable.contains'.")
    }
  }

  # 3. Check if variables are of desired type; if not, try to convert
  if (length(variable.class) > 0){
    for (i in 1:length(variable.class)){
      if (class(data[[names(variable.class)[i]]]) == variable.class[[i]]){
        message(paste0("✓ '", names(variable.class)[i], "' has desired class ",
                       variable.class[[i]], "."))
      } else {
        data[[names(variable.class)[i]]] = as(data[[names(variable.class)[i]]],
                                              variable.class[[i]])
        message(paste0("✓ '", names(variable.class)[i], "' has been converted to class ",
                       variable.class[[i]], "."))
      }
    }
  }
  return(data)
}

