#' Run subgroup analyses
#'
#' This function allows to simultaneously conduct different subgroup analyses using
#' \code{runMetaAnalysis} objects.
#'
#' @usage subgroupAnalysis(.model,
#'                  ...,
#'                  .which.run = .model$which.run[1],
#'                  .round.digits = 2,
#'                  .nnt.cer = NULL,
#'                  .tau.common = FALSE,
#'                  .html = TRUE)
#'
#' @param .model An object of class \code{"runMetaAnalysis"}, created by \code{\link{runMetaAnalysis}}.
#' @param ... <\link[dplyr]{dplyr_data_masking}>. A number of subgroup variables included in the original dataset
#' provided to \code{\link{runMetaAnalysis}}, separated by commas.
#' @param .which.run The model in \code{.model} that should be used for the subgroup analyses. Uses the default
#' analysis in \code{.model} if no value is specified by the user.
#' @param .round.digits \code{numeric}. Number of digits to round the (presented) results by. Default is \code{2}.
#' @param .nnt.cer \code{numeric}. Value between 0 and 1, indicating the assumed control group event rate to be used
#' for calculating NNTs via the Furukawa-Leucht method. If set to \code{NULL} (default),
#' the value saved in \code{.model} is (re-)used.
#' @param .tau.common \code{logical}. Should a common (\code{TRUE}) or subgroup-specific (\code{FALSE}) estimate
#' of the between-study heterogeneity be calculated when analyzing the subgroups? \code{FALSE} by default. Note that subgroup
#' analyses based on "multilevel" models automatically assume common heterogeneity estimates.
#' @param .html \code{logical}. Should an HTML table be created for the results? Default is \code{TRUE}.
#'
#' @return Returns an object of class \code{"subgroupAnalysis"}. This object includes, among other things,
#' a \code{data.frame} with the name \code{summary}, in which all subgroup analysis results are summarized.
#' Other objects are the "raw" subgroup analysis model objects returned.
#' This allows to conduct further operations on some subgroup analysis specifically.
#'
#' @examples
#' \dontrun{
#' data("depressionPsyCtr")
#' 
#' depressionPsyCtr %>%
#'   checkDataFormat() %>%
#'   checkConflicts() %>%
#'   calculateEffectSizes() %>% 
#'   filterPoolingData(condition_arm2 %in% 
#'                       c("wl", "other ctr")) -> data
#'
#' # Run the meta-analyses
#' runMetaAnalysis(data) -> res
#'
#' # Subgroup analysis
#' subgroupAnalysis(res, condition_arm2, country,
#'                  .which.run = "combined",
#'                  .tau.common = TRUE) -> sg
#' plot(sg, "condition_arm2")
#' plot(sg, "country")
#' }
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, 
#' Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{runMetaAnalysis}}
#'
#' @details For more details see the [Get Started](https://tools.metapsy.org/articles/metapsytools) vignette.
#'
#' @import dplyr
#' @importFrom crayon green yellow cyan bold
#' @importFrom scales pvalue
#' @importFrom purrr map
#' @importFrom meta metagen
#' @importFrom metafor escalc aggregate.escalc rma.mv
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @export subgroupAnalysis

subgroupAnalysis = function(.model, ...,
                            .which.run = .model$which.run[1],
                            .round.digits = 2,
                            .nnt.cer = NULL,
                            .tau.common = FALSE,
                            .html = TRUE){
  
  if (class(.model)[1] != "runMetaAnalysis"){
    stop("Input must be of class 'runMetaAnalysis'. Did you apply 'runMetaAnalysis' first?")
  }

  variables = .model$data %>% dplyr::select(...) %>% colnames()
  .round.digits = abs(.round.digits)
  .has.event.rate = (.model$.type.es=="CER"||.model$.type.es=="EER")[1]

  # Get model type
  model.type = paste0("model.", .which.run)
  M = .model[[model.type]]
  message("- ", crayon::green("[OK] "), "'", 
          model.type, "' used for subgroup analyses.")

  # Throw error if which.run is not included
  if (.which.run[1] == "combined" & !is.null(.model$model.combined$k)){
    if (.model$model.combined$k == 1){
      stop("'which.run' is set to 'combined', but there is only k=1 study/ES.")
    }
  }
  if (.which.run[1] == "lowest"){
    if (.model$model.lowest$k == 1){
      stop("'which.run' is set to 'lowest', but there is only k=1 study/ES.") 
    }
  }
  if (.which.run[1] == "highest"){
    if (.model$model.highest$k == 1){
      stop("'which.run' is set to 'highest', but there is only k=1 study/ES.")
    }
  }
  if (.which.run[1] == "influence"){
    if (.model$model.influence$k == 1){
      stop("'which.run' is set to 'influence', but there is only k=1 study/ES.")
    }
  }
  if (.which.run[1] == "rob"){
    if (.model$model.rob$k == 1){
      stop("'which.run' is set to 'rob', but there is only k=1 study/ES.")
    }
  }

  # Run all subgroup analyses
  if (class(M)[1] == "metagen" ||
      class(M)[1] == "metabin"){
    
    if (.tau.common)
      message("- ", crayon::green("[OK] "), 
              "Subgroup analyses conducted using a common heterogeneity variance estimate.")
    
    purrr::map(as.list(variables), function(x){
      
      na.mask = !is.na(M$data[[x]])
      meta::metagen(TE = M$TE[na.mask],
                    seTE = M$seTE[na.mask],
                    studlab = M$studlab[na.mask],
                    method.tau = M$method.tau,
                    method.tau.ci = M$method.tau.ci,
                    method.random.ci = ifelse(M$hakn, "HK", "classic"),
                    data = M$data[na.mask,],
                    subgroup = M$data[[x]][na.mask],
                    common = M$common,
                    random = M$random,
                    sm = switch(.model$.type.es, 
                                "CER" = "PLOGIT", 
                                "EER" = "PLOGIT", NULL),
                    tau.common = .tau.common)

    }) -> subgroup.analysis.list
    names(subgroup.analysis.list) = variables


    if (sum(is.na(M$data[variables])) > 0){
      warning("Some subgroup variables contained NA. These entries were omitted from model fitting.")
    }

    # Extract information
    subgroup.analysis.list %>%
      purrr::map2(as.list(variables), function(x,y){

        # Effect size in each group
        if (x$common == TRUE){
          g = round(
                if (identical(.model$.type.es, "RR")){
                   exp(x$TE.common.w)
                } else {
                  if (.has.event.rate) {plogis(x$TE.common.w)} else { x$TE.common.w }
                }, .round.digits)
          g.full = 
            if (identical(.model$.type.es, "RR")){
              exp(x$TE.common.w)
            } else {
              if (.has.event.rate) {plogis(x$TE.common.w)} else {x$TE.common.w}
            }
        } else {
          g = round(
                if (identical(.model$.type.es, "RR")){
                  exp(x$TE.random.w)
                } else {
                  if (.has.event.rate) {plogis(x$TE.random.w)} else {x$TE.random.w}
                }, .round.digits)
          g.full = 
            if (identical(.model$.type.es, "RR")){
              exp(x$TE.random.w)
            } else {
              if (.has.event.rate) {plogis(x$TE.random.w)} else {x$TE.random.w}
            }
        }

        # Confidence interval for g
        if (x$common == TRUE){
          g.ci = paste0("[", 
                        round(
                          if (identical(.model$.type.es, "RR")){
                            exp(x$lower.common.w)
                          } else {
                            if (.has.event.rate) {plogis(x$lower.common.w)} else {x$lower.common.w}
                          }, .round.digits), "; ",
                        round(
                          if (identical(.model$.type.es, "RR")){
                            exp(x$upper.common.w)
                          } else {
                            if (.has.event.rate) {plogis(x$upper.common.w)} else {x$upper.common.w}
                          }, .round.digits), "]")
        } else {
          g.ci = paste0("[", 
                        round(
                          if (identical(.model$.type.es, "RR")){
                            exp(x$lower.random.w)
                          } else {
                            if (.has.event.rate) {plogis(x$lower.random.w)} else {x$lower.random.w}
                          }, .round.digits), "; ",
                        round(
                          if (identical(.model$.type.es, "RR")){
                            exp(x$upper.random.w)
                          } else {
                            if (.has.event.rate) {plogis(x$upper.random.w)} else {x$upper.random.w}
                          }, .round.digits), "]")
        }

        # I-squared
        i2 = ifelse(is.na(x$I2.w), "-", round(x$I2.w*100, .round.digits-1))
        i2.ci = paste0("[", ifelse(is.na(x$lower.I2.w), "-",
                                   round(x$lower.I2.w*100, .round.digits-1)),
                       "; ", ifelse(is.na(x$upper.I2.w), "-",
                                    round(x$upper.I2.w*100, .round.digits-1)),"]") %>%
          ifelse(.=="[-; -]", "-", .)

        # NNT
        if (is.null(.nnt.cer)){
          metapsyNNT(abs(g.full), .model$nnt.cer) %>%
            round(.round.digits) %>% abs() -> nnt
        } else {
          metapsyNNT(abs(g.full), .nnt.cer) %>%
            round(.round.digits) %>% abs() -> nnt
        }

        
        if (M$common[1] & !M$random[1]) {
          pval.q.test = ifelse(is.na(x$pval.Q.b.common), NA,
                 scales::pvalue(x$pval.Q.b.common))
        } else {
          pval.q.test = ifelse(is.na(x$pval.Q.b.random), NA,
                 scales::pvalue(x$pval.Q.b.random))
        }

        data.frame(variable = y,
                   group = x$bylevs,
                   n.comp = x$k.w, g = g, g.ci = g.ci,
                   i2 = i2, i2.ci = i2.ci,
                   nnt = nnt,
                   p = pval.q.test)
        
        
      }) %>% do.call(rbind, .) %>% {rownames(.) = NULL;.} -> summary
    
    if (identical(.model$.type.es, "RR")){
      summary["nnt"] = "-"
      colnames(summary)[4:5] = c("rr", "rr.ci")
    }
    if (.has.event.rate){
      summary["nnt"] = "-"
      colnames(summary)[4:5] = c(tolower(.model$.type.es), 
                                 paste0(tolower(.model$.type.es),".ci"))
    }
  }

  if (class(M)[1] == "rma.mv"){

    stringr::str_replace_all(as.character(M$random[[1]]), "1 \\| ", "")[2] %>%
      strsplit("/") %>% {.[[1]]} %>% {.[1]} -> study.id

    dat.mv = data.frame(yi = M$yi,
                        vi = M$vi,
                        slab = M$slab,
                        study = M$data[[study.id]],
                        es.id = 1:length(M$yi))
    dat.mv = cbind(dat.mv, .model$data[!is.na(M$yi),] %>%
                     dplyr::select(dplyr::all_of(variables)) %>%
                     purrr::map_dfr(~as.factor(.)))

    purrr::map(as.list(variables), function(x){
      form = as.formula(paste0("~ 0 +", x))
      tryCatch2({
        metafor::rma.mv(yi = yi,
                        V = vi,
                        data = dat.mv,
                        random = M$random[[1]],
                        test = M$test,
                        method = "REML",
                        mods = form)}) %>% 
        {.$value} -> res.mv
      }) -> subgroup.analysis.list
    names(subgroup.analysis.list) = variables


    purrr::map2(subgroup.analysis.list, variables, 
                function(x, y){
      
      if (identical(.model$.type.es, "RR")){
        g = c(exp(as.numeric(x$b))) %>% round(.round.digits)
        g.full = c(exp(as.numeric(x$b)))
        g.lower = {exp(x$ci.lb)} %>% round(.round.digits)
        g.upper = {exp(x$ci.ub)} %>% round(.round.digits)
        g.ci = paste0("[", g.lower, "; ", g.upper, "]")
      } else {
        if (.has.event.rate) {
          g = as.numeric(x$b) %>% plogis() %>% round(.round.digits)
          g.full = as.numeric(x$b) %>% plogis()
          g.lower = {x$ci.lb} %>% plogis() %>% round(.round.digits)
          g.upper = {x$ci.ub} %>% plogis() %>% round(.round.digits)
          g.ci = paste0("[", g.lower, "; ", g.upper, "]")
        } else {
          g = as.numeric(x$b) %>% round(.round.digits)
          g.full = as.numeric(x$b)
          g.lower = {as.numeric(x$ci.lb)} %>% round(.round.digits)
          g.upper = {as.numeric(x$ci.ub)} %>% round(.round.digits)
          g.ci = paste0("[", g.lower, "; ", g.upper, "]")
        }
      }

      if (is.null(.nnt.cer)){
        metapsyNNT(abs(g.full), .model$nnt.cer) %>%
          round(.round.digits) %>% abs() -> nnt
      } else {
        metapsyNNT(abs(g.full), .nnt.cer) %>%
          round(.round.digits) %>% abs() -> nnt
      }

      data.frame(variable = y,
                 group = levels(dat.mv[[y]]),
                 n.comp = table(dat.mv[[y]]) %>% as.numeric(),
                 g = g, g.ci = g.ci, i2 = "-", i2.ci = "-",
                 nnt = nnt, p = x$QMp %>% scales::pvalue())
    }) %>% do.call(rbind, .) %>% {rownames(.) = NULL;.} -> summary
    
    if (identical(.model$.type.es, "RR")){
      colnames(summary)[4:5] = c("rr", "rr.ci")
    }
    if (.has.event.rate) {
      colnames(summary)[4:5] = c(tolower(.model$.type.es), paste0(tolower(.model$.type.es), ".ci"))
    }

  }


  # Replace NAs
  summary[is.na(summary)] = "-"
  summary[summary == "[; ]" &
            !is.na(summary)] = "[-;-]"
  summary[summary == "" &
            !is.na(summary)] = "-"
  
  if (identical(.model$.type.es, "RR")){
    summary$nnt = "-"
  }
  if (.has.event.rate) {
    summary$nnt = "-"
  }
  
  # Return
  returnlist = list(summary = summary,
                    subgroup.analysis.list = subgroup.analysis.list,
                    html = .html, 
                    model.type = class(M)[1],
                    .type.es = .model$.type.es)
  class(returnlist) = c("subgroupAnalysis", "list")
  return(returnlist)

}




#' Print method for objects of class 'subgroupAnalysis'
#'
#' Print S3 method for objects of class \code{subgroupAnalysis}.
#'
#' @param x An object of class \code{subgroupAnalysis}.
#' @param ... Additional arguments.
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @importFrom knitr kable
#' @importFrom dplyr as_tibble
#' @importFrom kableExtra kable_styling column_spec footnote
#' @importFrom crayon green blue magenta bold
#' @importFrom stringr str_sub
#'
#' @export
#' @method print subgroupAnalysis

print.subgroupAnalysis = function(x, ...){

  unique(with(x$summary, 
         unlist(lapply(variable, 
                       function(i) grep(i, variable)[1])))
    ) -> first.var.mask
  x$summary$variable[!1:nrow(x$summary) %in% first.var.mask] = "."
  x$summary$p[!1:nrow(x$summary) %in% first.var.mask] = "."
  
  cat(crayon::blue$bold("Subgroup analysis results "))
  cat(crayon::blue$bold(
    "---------------------- \n"))
  dat = x$summary
  tbl = dplyr::as_tibble(dat)
  old = options(pillar.bold=TRUE)
  tbl.format = format(tbl)[-c(1,3)]
  tbl.format[-1] = lapply(tbl.format[-1], 
                          function(x) stringr::str_sub(x, 19))
  tbl.format[1] = stringr::str_sub(tbl.format[1], 3)
  cat(do.call(c, tbl.format), sep="\n")
  options(old)
  
  if (x$html == TRUE){
    
    if (identical(x$.type.es, "RR")){
      colNames = c("Variable", "Level", "<i>n</i><sub>comp</sub>",
                   "<i>RR</i>", "CI", "<i>I</i><sup>2</sup>",
                   "CI", "NNT", "<i>p</i>")
    } else {
      if (identical(x$.type.es, "EER")||identical(x$.type.es, "CER")) {
        colNames = c("Variable", "Level", "<i>n</i><sub>comp</sub>",
                     paste0("<i>", x$.type.es, "</i>"), "CI", "<i>I</i><sup>2</sup>",
                     "CI", "NNT", "<i>p</i>")
      } else {
        colNames = c("Variable", "Level", "<i>n</i><sub>comp</sub>",
                     "<i>g</i>", "CI", "<i>I</i><sup>2</sup>",
                     "CI", "NNT", "<i>p</i>") 
      }
    }

    x$summary %>%
      {.$p = ifelse(.$p == "<0.001", "&lt;0.001", .$p);.} %>% 
      {colnames(.) = colNames;.} %>%
      knitr::kable(escape = FALSE, format = "html") %>%
      kableExtra::kable_styling(font_size = 8, full_width = FALSE) %>%
      kableExtra::column_spec(1, bold = TRUE) %>%
      kableExtra::collapse_rows(columns = 1, valign = "top") %>%
      print()
  }
}


#' Plot method for objects of class 'runMetaAnalysis'
#'
#' Plot S3 method for objects of class \code{runMetaAnalysis}.
#'
#' @param x An object of class \code{runMetaAnalysis}.
#' @param which \code{character}. Subgroup analysis to be plotted (variable name).
#' @param ... Additional arguments.
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @importFrom meta forest
#'
#' @export
#' @method plot subgroupAnalysis

plot.subgroupAnalysis = function(x, which = NULL, ...){

  if (x$model.type[1] == "rma.mv"){
    stop("Cannot generate subgroup analysis forest plots for 'threelevel' models.")
  }

  if (is.null(which)){
    message("- ", crayon::green("[OK] "), "'", 
            names(x$subgroup.analysis.list)[1], "' used for forest plot.")
    if (identical(x$.type.es, "RR")){
      x$subgroup.analysis.list[[1]] = updateMeta(
        x$subgroup.analysis.list[[1]], sm = "RR")
    }
    meta::forest(x$subgroup.analysis.list[[1]], layout = "JAMA")
  } else {
    if (identical(x$.type.es, "RR")){
      x$subgroup.analysis.list[[which[1]]] = updateMeta(
        x$subgroup.analysis.list[[which[1]]], sm = "RR")
    }
    meta::forest(x$subgroup.analysis.list[[which[1]]], layout = "JAMA")
  }
}
