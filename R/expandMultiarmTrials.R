# data: meta-analysis data stored as a data.frame
# vars.for.id: character vector, containing column names of all variables used to
# construct unique comparison IDs.
# study.indicator: character, signifying the name of the variable containing the study name.
# multiarm.indicator: numeric, signifying if a row is part of a multiarm study (1) or not (0).
# multiarm.group.indicator: character, signifying the name of the variable containing the name of the (active) treatment arm.
# should be NA when study is not a multiarm trial/the row is part of the control group.
# no.arms.indicator: character, signifying the name of the variable containing the number of arms included in a study (typically 2).
# group.indicator: character, column name of the variable storing the study name
# group.names: list, storing the name of the value corresponding to the intervention group ("ig") and control group ("cg").

# dplyr `%>%` filter arrange pull
# stringr str_replace_all str_remove_all
# magrittr set_rownames

expandMultiarmTrials = function(data,
                                vars.for.id = c("study", "primary",
                                                "Outc_measure",
                                                "Time", "Time_weeks"),
                                study.indicator = "study",
                                multiarm.indicator = "is.multiarm",
                                multiarm.group.indicator = "multiple_arms",
                                no.arms.indicator = "no.arms",
                                group.indicator = "condition",
                                group.names = list("ig" = "ig",
                                                   "cg" = "cg")){

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
    str_replace_all(",", "_") %>% str_remove_all(" ") -> data$id

  # Filter out multiarm studies
  data[data[names(data) == multiarm.indicator] == 1,] %>%
    dplyr::arrange(id, group.indicator) -> data.multiarm

  # Create IG arm names in multiarm studies
  data.multiarm[
    data.multiarm[names(data.multiarm) == group.indicator] ==
      group.names[["ig"]],] %>%
    pull(multiarm.group.indicator) -> arm.names

  # Expand data set so that each comparison has two lines (unique comparisons,
  # including multi-arm trials)
  i.ma = which(data.multiarm[,names(data.multiarm) == group.indicator] == group.names[["cg"]])

  # Create CG arms via duplication
  n.reps = data.multiarm[i.ma,
                         names(data.multiarm)[names(data.multiarm) ==
                                                no.arms.indicator]] %>%
                                                          pull(1) %>% {. -1}
  data.multiarm[rep(i.ma, n.reps),] %>%
    set_rownames(NULL) -> data.cg
  data.cg[[study.indicator]] = paste0(data.cg[[study.indicator]], " -", arm.names)

  # Prepare non-CG arms
  data.multiarm[data.multiarm[,names(data.multiarm) ==
                                group.indicator] ==
                  group.names[["ig"]],] -> data.igs

  data.igs[[study.indicator]] = paste0(data.igs[[study.indicator]],
                                      " -", data.igs[[multiarm.group.indicator]])


  # Combine all multiarm data
  data.multiarm = rbind(data.igs, data.cg)
  apply(data.multiarm, 1,
        function(x){paste(as.character(x[vars.for.id]), collapse = "_")}) %>%
    str_replace_all(",", "_") %>% str_remove_all(" ") -> data.multiarm$id
  data.multiarm = dplyr::arrange(data.multiarm, id)

  # Add multiarm and non-multiarm data
  rbind(data[data[[multiarm.indicator]] == 0,], data.multiarm) -> data
  data = dplyr::arrange(data, id)

  # Set class
  class(data) = c("data.frame", "expandMultiarmTrials")

  # Return
  return(data)

}


