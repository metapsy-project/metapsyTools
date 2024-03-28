#' Explore included treatments and comparisons
#'
#' @description This function allows to summarize included treatments and 
#' treatment comparisons in a data set.
#'
#' @usage exploreStudies(data,
#'                which = c("treatments", "comparisons"),
#'                
#'                # Metapsy standard variables
#'                .study.var = "study",
#'                .condition = "condition",
#'                .condition.specification = "multi",
#'                .groups.column.indicator = c("_arm1", "_arm2"),
#'                .trt.indicator = "arm",
#'                .n.vars = c("n", "n_change", "totaln", "N"),
#'                
#'                # Output
#'                html = TRUE)
#'
#' @param data \code{data.frame}. Effect size data in the wide format, as 
#' created by \code{\link{calculateEffectSizes}}. For the other default settings to be
#' applicable, the data set should follow the [Metapsy data standard](https://docs.metapsy.org/data-preparation/format/).
#' Alternatively, one can also provide an `metapsyDatabase` object as returned by [metapsyData::getData()],
#' or a meta-analysis object returned by [runMetaAnalysis()].
#' @param which Should the data set be summarized with respect to the included
#' treatments (`"treatments"`) or treatment comparisons (`"comparisons"`)? Defaults to
#' `"treatments"`.
#' @param .study.var `character`. The name of the variable in the data set in which
#' the study labels are stored. 
#' @param .condition \code{character}. The prefix of the two variables in \code{data} in 
#' which the conditions (e.g. "guided iCBT", "waitlist") of the trial arm comparison are stored.
#' @param .condition.specification \code{character}. The prefix of the two variables in the dataset
#' which provide a "specification" of the trial arm condition in multiarm trials.
#' @param .groups.column.indicator \code{character}. A character vector with two elements, 
#' representing the suffix used to differentiate between the first and second arm in a comparison.
#' @param .trt.indicator \code{character}. A character specifying the name used to indicate the treatment arm.
#' @param .n.vars `character`. A character vector which includes the names of all variables in the data set in
#' which sample size information is stored. Only the prefix is needed, where `.groups.column.indicator`
#' provides the suffixes.
#' @param html \code{logical}. Should an HTML table be created for the results? Default is \code{TRUE}.
#'
#' @return Returns an object of class \code{"exploreStudies"}. This object includes a list object called `summary`
#' in which the counts for distinct treatments (`conditions`) and comparisons (`comparisons`) are summarized, as well as
#' a `data.frame`  data. This data frame includes the initially provided data set collapsed by study 
#' (so that each row represents one study). To this data set, variables are added that encode how many arms with a specific condition
#' are included in the trial (e.g. if `cbt=2`, this means that two CBT groups are included in the trial),
#' as well as the number of distinct comparisons, and the sample size of both (these columns all start with `n.`).
#' This can be helpful to perform further descriptive analyses.
#' 
#' @details Using the variables provided in the `.n.vars` argument, `exploreStudies` calculates the arm- and study-specific
#' sample sizes. If no adequate information is provided, sample sizes cannot be calculated for a study. If this
#' is the case, a warning is printed, pointing to the studies with missing sample size information.
#'
#' @examples
#' \dontrun{
#' # Explore studies in built-in dataset
#' data("depressionPsyCtr")
#' exploreStudies(depressionPsyCtr, "treatments") 
#' exploreStudies(depressionPsyCtr, "comparisons") 
#' 
#' # - Extract metapsy database using metapsyData
#' # - Filter CBT and PST studies
#' # - Run a meta-analysis and explore synthesize studies
#' library(metapsyData)
#' getData("depression-psyctr", version="22.0.2") %>% 
#'   filterPoolingData(condition_arm1 %in% c("cbt", "pst")) %>% 
#'   runMetaAnalysis(which.run = c("combined")) -> res
#' 
#' exploreStudies(res)
#' exploreStudies(res, "comparisons")
#' }
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{createStudyTable}}, \code{\link{calculateEffectSizes}}, 
#' \code{\link{subgroupAnalysis}}, \code{\link{correctPublicationBias}},
#' \code{\link{metaRegression}}, \code{\link{runMetaAnalysis}}.
#'
#' @importFrom crayon green blue bold
#' @importFrom scales pvalue
#' @importFrom stats dffits model.matrix rnorm rstudent complete.cases median quantile aggregate
#' @export exploreStudies

exploreStudies = function(data,
                          which = c("treatments", "comparisons"),
                          .study.var = "study",
                          .condition = "condition",
                          .condition.specification = "multi",
                          .groups.column.indicator = c("_arm1", "_arm2"),
                          .trt.indicator = "arm",
                          .n.vars = c("n", "n_change", "totaln", "N"),
                          html = TRUE){
  
  if (identical(class(data)[1], "runMetaAnalysis")){
    data$data -> dat
  } else if (identical(class(data)[1], "metapsyDatabase")) {
    data$data -> dat
  } else {
    dat = data
  }
  
  message("- ", crayon::green("[OK] "), 
          "Generating summary table for '", which[1],"'... ", 
          appendLF = FALSE)
  
  # # # # # # # # # # # # # # # # # # # # # #
  #  CONDITIONS                             #
  # # # # # # # # # # # # # # # # # # # # # #
  
  # Extract number of studies per condition & conditions per study
  split(dat, dat[[.study.var]]) %>% 
    lapply(function(x){
      v1 = x[[paste0(.condition.specification, .groups.column.indicator[1])]]
      v2 = x[[paste0(.condition.specification, .groups.column.indicator[2])]]
      names(v1) = x[[paste0(.condition, .groups.column.indicator[1])]]
      names(v2) = x[[paste0(.condition, .groups.column.indicator[2])]]
      vec = c(v1, v2)
      if (all(is.na(vec))){ unique(names(vec)) } else { 
        names(vec[!duplicated(vec)]) }
    }) -> condition.list.ma
  
  all.conditions.ma = unique(unlist(condition.list.ma))
  
  sapply(all.conditions.ma, function(x){
    lapply(condition.list.ma, function(y){sum(y %in% x)})
  }) %>% as.matrix() -> condition.mat.ma
  
  # Get studies per condition
  condition.mat.ma %>% 
    apply(2, function(z){
      sum(unlist(z>0))
    }) -> condition.k
  
  # Get condition per study
  condition.mat.ma %>% 
    apply(2, function(z){
      sum(unlist(z))
    }) -> condition.n
  
  
  # Get sample size per condition
  n.names.1 = c(paste0(.n.vars, .groups.column.indicator[1]))
  n.names.2 = c(paste0(.n.vars, .groups.column.indicator[2]))
  
  split(dat, dat[[.study.var]]) %>% 
    lapply(function(B){
      data.frame(
        condition_arm1 = B[[paste0(.condition, .groups.column.indicator[1])]],
        condition_arm2 = B[[paste0(.condition, .groups.column.indicator[2])]],
        multi_arm1 = B[[paste0(.condition.specification, .groups.column.indicator[1])]],
        multi_arm2 = B[[paste0(.condition.specification, .groups.column.indicator[2])]],
        N1 = apply(B[colnames(B) %in% n.names.1], 1, function(x){
          x[!is.na(x)][1]}) %>% as.numeric(),
        N2 = apply(B[colnames(B) %in% n.names.2], 1, function(x){
          x[!is.na(x)][1]}) %>% as.numeric()) %>% 
        within({
          .condition_arm1 = paste0(condition_arm1, multi_arm1)
          .condition_arm2 = paste0(condition_arm2, multi_arm2)
        }) %>% 
        with(data.frame(condition = c(condition_arm1, condition_arm2),
                        .condition = c(.condition_arm1, .condition_arm2),
                        N = c(N1, N2))) %>% 
        with({
          aggregate(N, by=list(g1 = condition, g2 = .condition), FUN=mean)}) %>% 
        with({aggregate(x, by=list(g=g1), FUN=sum)}) %>% 
        {.$x -> res; names(res) = .$g; res}
    }) -> n.res
  
  sapply(all.conditions.ma, function(x){
    lapply(n.res, function(y){y[x]})
  }) %>% as.matrix() %>% 
    {.[is.na(.)] = 0;.} -> condition.mat.n
  
  condition.mat.n %>% 
    apply(2, function(z){
      sum(unlist(z), na.rm = TRUE)
    }) %>% round() -> sample.n
  
  
  # # # # # # # # # # # # # # # # # # # # # #
  #  COMPARISON                             #
  # # # # # # # # # # # # # # # # # # # # # #
  
  # Extract number of studies per comparison & comparisons per study
  split(dat, dat[[.study.var]]) %>% 
    lapply(function(B){
      data.frame(
        condition = paste(B[[paste0(.condition, .groups.column.indicator[1])]], "-",
                          B[[paste0(.condition, .groups.column.indicator[2])]]),
        multi = paste(B[[paste0(.condition.specification, .groups.column.indicator[1])]], 
                      B[[paste0(.condition, .groups.column.indicator[1])]], "-",
                      B[[paste0(.condition.specification, .groups.column.indicator[2])]],
                      B[[paste0(.condition, .groups.column.indicator[2])]]),
        count = 1) %>% 
        with(aggregate(count, by=list(g1 = condition, g2 = multi), FUN = mean)) %>% 
        with(aggregate(x, by=list(g = g1), FUN = sum)) %>% 
        {.$x -> res; names(res) = .$g; res}
    }) -> comp.list
  
  all.comps = lapply(comp.list, function(z) names(z)) %>% unlist() %>% unique()
  
  sapply(all.comps, function(x){
    lapply(comp.list, function(y){sum(y[x])})
  }) %>% as.matrix() %>% {.[is.na(.)] = 0;.} -> comp.mat
  
  # Get studies per comparison
  comp.mat %>% 
    apply(2, function(z){
      sum(unlist(z>0))
    }) -> comp.k
  
  # Get comparisons per study
  comp.mat %>% 
    apply(2, function(z){
      sum(unlist(z))
    }) -> comp.n
  
  # Get sample size per comparison
  n.names.1 = c(paste0(.n.vars, .groups.column.indicator[1]))
  n.names.2 = c(paste0(.n.vars, .groups.column.indicator[2]))
  
  split(dat, dat[[.study.var]]) %>% 
    lapply(function(B){
      data.frame(
        condition = paste(B[[paste0(.condition, .groups.column.indicator[1])]], "-",
                          B[[paste0(.condition, .groups.column.indicator[2])]]),
        multi = paste(B[[paste0(.condition.specification, .groups.column.indicator[1])]], 
                      B[[paste0(.condition, .groups.column.indicator[1])]], "-",
                      B[[paste0(.condition.specification, .groups.column.indicator[2])]],
                      B[[paste0(.condition, .groups.column.indicator[2])]]),
        N = apply(B[colnames(B) %in% n.names.1], 1, function(x){x[!is.na(x)][1]}) %>% as.numeric() + 
            apply(B[colnames(B) %in% n.names.2], 1, function(x){x[!is.na(x)][1]}) %>% as.numeric()) %>% 
        with({
          aggregate(N, by=list(g1 = condition, g2 = multi), FUN=mean)}) %>% 
        with({aggregate(x, by=list(g=g1), FUN=sum)}) %>% 
        {.$x -> res; names(res) = .$g; res}
    }) -> n.comp.res
  
  sapply(all.comps, function(x){
    lapply(n.comp.res, function(y){y[x]})
  }) %>% as.matrix() -> comp.mat.n
  
  comp.mat.n %>% 
    apply(2, function(z){
      sum(unlist(z), na.rm = TRUE)
    }) %>% round() -> comp.sample.n
  
  
  # # # # # # # # # # # # # # # # # # # # # #
  #  COMPILE FOR RETURN                     #
  # # # # # # # # # # # # # # # # # # # # # #
  
  # Summary table
  summary.comp = data.frame(
    studies = comp.k,
    conditions = comp.n[names(comp.k)],
    patients = comp.sample.n[names(comp.k)]) %>% 
    {.[order(rownames(.)),]}
  
  summary.condition = data.frame(
    studies = condition.k,
    conditions = condition.n[names(condition.k)],
    patients = sample.n[names(condition.k)])
  
  # Create collapsed data set
  colnames(comp.mat.n) = paste0("n.", colnames(comp.mat.n))
  comp.mat.n[is.na(comp.mat.n)] = 0
  colnames(condition.mat.n) = paste0("n.", colnames(condition.mat.n))
  condition.mat.n[is.na(condition.mat.n)] = 0
  
  data = cbind(study = dat[!duplicated(dat$study),][[.study.var]],
               condition.mat.ma, condition.mat.n,
               comp.mat, comp.mat.n,
               dat[!duplicated(dat$study), setdiff(colnames(dat), "study")]) %>% 
    as.data.frame() %>% 
    {rownames(.) = NULL;.}
  
  
  # Check for missing sample size info
  condition.mat.n %>% 
    apply(1, function(z) 
      sum(unlist(z!=0))) < 2 -> incomplete
  
 if (!all(!incomplete)){
    warning("No sample size information found for studies ",
            paste(names(incomplete)[incomplete], collapse = "; "), ".")
 }
  
  # Return
  retlist = list(summary = 
                   list(conditions = summary.condition,
                        comparisons = summary.comp),
                 data = data, html = html[1],
                 which = which[1])
  
  
  cat(crayon::green("DONE \n"))
  class(retlist) = c("exploreStudies", "list")
  return(retlist)
  
}



#' Print method for objects of class 'exploreStudies'
#'
#' Print S3 method for objects of class \code{exploreStudies}.
#'
#' @param x An object of class \code{exploreStudies}.
#' @param ... Additional arguments.
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers 
#' \email{p.cuijpers@@vu.nl}
#'
#' @importFrom knitr kable
#' @importFrom dplyr as_tibble
#' @importFrom kableExtra kable_styling column_spec footnote
#' @importFrom crayon green blue magenta bold
#' @importFrom stringr str_sub
#'
#' @export
#' @method print exploreStudies

print.exploreStudies = function(x, ...){
  
  studies.total = length(x$data$study)
  conditions.total = sum(x$summary$conditions$conditions)
  patients.total = sum(x$summary$conditions$patients)
  
  if (identical(x$which[1], "comparison")||identical(x$which[1], "comparisons")){
    
    x$summary$comparisons$studies = 
      paste0(x$summary$comparisons$studies, " (",
             round(x$summary$comparisons$studies/studies.total,4)*100, "%)")
    x$summary$comparisons$conditions = 
      paste0(x$summary$comparisons$conditions, " (",
             round(x$summary$comparisons$conditions/conditions.total,4)*100, "%)")
    x$summary$comparisons$patients = 
      paste0(x$summary$comparisons$patients, " (",
             round(x$summary$comparisons$patients/patients.total,4)*100, "%)")
    
    cat(crayon::blue$bold("Summary "))
    cat(crayon::blue$bold("------------------------------------------------ \n"))
    cat("- Total number of studies: k=", studies.total, "\n", sep="")
    cat("- Total number of conditions/trial arms: n=", conditions.total, "\n", sep="")
    cat("- Total number of patients: n=", patients.total, "\n", sep=""); cat("\n")
    cat(crayon::blue$bold("Comparisons "))
    cat(crayon::blue$bold("-------------------------------------------- \n"))
    print(x$summary$comparisons)
    
    if (x$html[1]){
      footnote = c(
        paste("Total number of studies: k=", studies.total, sep=""),
        paste("Total number of conditions/trial arms: n=", conditions.total, sep=""),
        paste("Total number of patients: n=", patients.total, sep=""))
      
      x$summary$comparisons %>% 
        knitr::kable(escape = FALSE, format = "html") %>%
        kableExtra::kable_styling(font_size = 8, full_width = FALSE) %>%
        kableExtra::column_spec(1, bold = TRUE, width_min = "13em") %>%
        kableExtra::footnote(alphabet = footnote) %>% 
        print()
    }
    
  } else {
    
    x$summary$conditions$studies = 
      paste0(x$summary$conditions$studies, " (",
             round(x$summary$conditions$studies/studies.total,4)*100, "%)")
    x$summary$conditions$conditions = 
      paste0(x$summary$conditions$conditions, " (",
             round(x$summary$conditions$conditions/conditions.total,4)*100, "%)")
    x$summary$conditions$patients = 
      paste0(x$summary$conditions$patients, " (",
             round(x$summary$conditions$patients/patients.total,4)*100, "%)")
    
    cat(crayon::blue$bold("Summary "))
    cat(crayon::blue$bold("------------------------------------------------ \n"))
    cat("- Total number of studies: k=", studies.total, "\n", sep="")
    cat("- Total number of conditions/trial arms: n=", conditions.total, "\n", sep="")
    cat("- Total number of patients: n=", patients.total, "\n", sep=""); cat("\n")
    cat(crayon::blue$bold("Treatments "))
    cat(crayon::blue$bold("--------------------------------------------- \n"))
    print(x$summary$conditions)
    
    if (x$html[1]){
      footnote = c(
        paste("Total number of studies: k=", studies.total, sep=""),
        paste("Total number of conditions/trial arms: n=", conditions.total, sep=""),
        paste("Total number of patients: n=", patients.total, sep=""))
      
      x$summary$conditions %>% 
        knitr::kable(escape = FALSE, format = "html") %>%
        kableExtra::kable_styling(font_size = 8, full_width = FALSE) %>%
        kableExtra::column_spec(1, bold = TRUE, width_min = "13em") %>%
        kableExtra::footnote(alphabet = footnote) %>% 
        print()
    }
  }
}

