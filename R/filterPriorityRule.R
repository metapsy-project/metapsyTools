psyCtrSubset %>%
  checkDataFormat() %>%
  checkConflicts() %>%
  expandMultiarmTrials() %>%
  calculateEffectSizes() -> data

filterPriorityRule(data,
                   Cond_spec_trt1 = c("cbt", "pst"),
                   Cond_spec_trt2 = c("cau", "wl", "cbt"),
                   Outc_measure = c("cesd", "phq-9", "scl", "hdrs"),
                   Time = c("post", "fu")) -> res

filterPriorityRule = function(data, ..., study.indicator = "study"){

  rules = enquos(...)
  vars = names(rules)

  for (i in 1:length(rules)){

    data %>%
      split(.[[study.indicator]]) %>%
      purrr::map_df(function(x){

        data.frame(factor = rlang::eval_tidy(rules[[i]]),
                   weight = length(rlang::eval_tidy(rules[[i]])):1) -> weight

        data.frame(factor = x[[vars[i]]],
                   dat = as.numeric(x[[vars[i]]] %in%
                                      rlang::eval_tidy(rules[[i]]))) -> dat

        merge(dat, weight, by = "factor", all.x = TRUE) -> tab
        tab[is.na(tab)] = 0

        if (sum(with(tab, {dat*weight})) > 0){
          with(tab, {dat*weight == max(dat*weight)}) -> mask
          x[x[[vars[i]]] %in% unique(tab[mask, "factor"]),]
        } else {
          x = x[NULL,]
        }

      }) -> data
  }

  return(data)
}

filterPriorityRule(data,
                   Cond_spec_trt2 = c("cau", "wl", "cbt"),
                   Outc_measure = c("cesd", "bdi-ii", "hamd"),
                   Time = c("FU", "post")) -> res

runMetaAnalysis(res) %>% plot()

