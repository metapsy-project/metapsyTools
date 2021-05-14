# data: meta-analysis dataset of class expandMultiarmTrials, created using expanMultiarmTrials
# comp.id.indicator: character, column name of the variable storing the comparison ID
# group.indicator: character, column name of the variable storing the study name
# group.indicator: character, column name of the variable storing the condition
# (intervention or control group)
# pivot.vars: character vector, contains the name of the comparison ID variable, the condition variable,
# and all raw effect size data variables from which effect sizes should be calculated
# funcs: list of functions. These functions will be used to calculate the effect sizes (Hedges' g) based on the raw data.
# By default, g is calculated from (i) the mean, SD and N; (ii) binary outcome data (i.e. response counts),
# and (iii) change scores. Other functions can be added to the list. Results of the function must result in a data.frame
# with the same number of rows as in data, and two columns: one for the calculated g value (named es) and its standard error (named se).
# In rows for which no fitting raw data was supplied, the resulting data.frame should contain NA.

# tidyr pivot_wider pivot_longer
# dplyr `%>%` all_of select filter mutate arrange
# purrr pmap_dfr
# esc esc_mean_sd

calculateEffectSizes = function(data,
                                comp.id.indicator = "id",
                                study.indicator = "study",
                                group.indicator = "condition",
                                pivot.vars = c("id", "condition",
                                               "Post_M", "Post_SD", "Post_N",
                                               "Rand_N", "Improved_N",
                                               "Change_m", "Change_SD",
                                               "Change_N"),
                                funcs = list(
                                  mean.sd = function(x, ...){
                                    x %>%
                                      purrr::pmap_dfr(function(Post_M_ig, Post_M_cg, Post_SD_ig,
                                                               Post_SD_cg, Post_N_ig, Post_N_cg, ...)
                                      {esc::esc_mean_sd(Post_M_ig, Post_SD_ig, Post_N_ig,
                                                        Post_M_cg, Post_SD_cg, Post_N_cg, es.type = "g") %>%
                                          as.data.frame() %>% dplyr::select(es, se) %>%
                                          suppressWarnings()})
                                  },
                                  binary = function(x, ...){
                                    x %>%
                                      purrr::pmap_dfr(function(Improved_N_ig, Improved_N_cg, Rand_N_ig, Rand_N_cg, ...)
                                      {esc::esc_2x2(Improved_N_ig,
                                                    Rand_N_ig - Improved_N_ig,
                                                    Improved_N_cg,
                                                    Rand_N_cg - Improved_N_cg,
                                                    es.type = "g") %>%
                                          as.data.frame() %>% dplyr::select(es, se) %>%
                                          suppressWarnings() %>%
                                          mutate(es = es*-1)})
                                  },
                                  change = function(x, ...){
                                    x %>%
                                      purrr::pmap_dfr(function(Change_m_ig, Change_m_cg, Change_SD_ig,
                                                               Change_SD_cg, Change_N_ig, Change_N_cg, ...)
                                      {esc::esc_mean_sd(Change_m_ig, Change_SD_ig, Change_N_ig,
                                                        Change_m_cg, Change_SD_cg, Change_N_cg, es.type = "g") %>%
                                          as.data.frame() %>% dplyr::select(es, se) %>%
                                          suppressWarnings()})})
                                ){

  # check class
  if (class(data)[2] != "expandMultiarmTrials"){
    message(paste("class of 'data' is not 'expandMultiarmTrials'.",
            "Sure that the data has the right format?"))
  }

  # Convert to data.frame (conventional)
  data = data.frame(data)

  # Define id
  data$id = data[[comp.id.indicator]]

  # Convert to wide
  data %>%
    dplyr::select(dplyr::all_of(pivot.vars)) %>%
    tidyr::pivot_wider(names_from = dplyr::all_of(group.indicator),
                values_from = c(-id, -dplyr::all_of(group.indicator))) %>%
    {.[,colSums(is.na(.))<nrow(.)]} -> data.wide

  # Apply funcs
  message("calculating effect sizes...")
  es.res = list()
  for (i in 1:length(funcs)){
    es.res[[i]] = funcs[[i]](data.wide)
  }
  es.res = do.call(cbind, es.res)
  message("SUCCESS")

  # Now, bind all calculated ES together,
  # then bind together with wide dataset
  es.res %>%
    apply(., 1, function(x) x[!is.na(x)]) %>% t() %>%
    cbind(data.wide, .) -> data.wide.es

  # Now, transform the wide format data set with calculated ES
  # back to long and merge back with original version
  data.wide.es %>%
    dplyr::select(id, es, se) %>%
    dplyr::mutate(trt1 = "t.1", trt2 = "t.2") %>%
    tidyr::pivot_longer(-id,
                 names_to = c(".value"),
                 names_pattern = "(..)") %>%
    dplyr::arrange(id) %>%
    dplyr::select(2:3) %>%
    cbind(data %>% arrange(id) %>% select(-id), .) -> dat.final

  # Return
  return(dat.final)
}








