# subgroupAnalysis = function(.model, ...,
#                  .which.run = .model$which.run[1],
#                  .round.digits = 2,
#                  .nnt.cer = NULL,
#                  .html = TRUE){
#
#   if (class(.model)[1] != "runMetaAnalysis"){
#     stop("Input must be of class 'runMetaAnalysis'. Did you apply 'runMetaAnalysis' first?")
#   }
#
#   variables = .model$data %>% dplyr::select(...) %>% colnames()
#
#   # Get model type
#   model.type = paste0("model.", .which.run)
#   M = .model[[model.type]]
#   message("- [OK] '", model.type, "' used for subgroup analyses.")
#
#   # Run all subgroup analyses
#   if (class(M)[1] == "metagen"){
#
#     purrr::map(as.list(variables), function(x){
#       meta::update.meta(M, byvar = M$data[[x]])
#     }) -> subgroup.analysis.list
#     names(subgroup.analysis.list) = variables
#
#     # Extract information
#     subgroup.analysis.list %>%
#       purrr::map2(as.list(variables), function(x,y){
#
#         # Effect size in each group
#         if (x$comb.fixed == TRUE){g = round(x$TE.fixed.w, .round.digits)} else {
#           g = round(x$TE.random.w, .round.digits)}
#
#         # Confidence interval for g
#         if (x$comb.fixed == TRUE){
#           g.ci = paste0("[", round(x$lower.fixed.w, .round.digits), "; ",
#                         round(x$upper.fixed.w, .round.digits), "]")
#         } else {
#           g.ci = paste0("[", round(x$lower.random.w, .round.digits), "; ",
#                         round(x$upper.random.w, .round.digits), "]")}
#
#         # I-squared
#         i2 = ifelse(is.na(x$I2.w), "-", round(x$I2.w*100, .round.digits-1))
#         i2.ci = paste0("[", ifelse(is.na(x$lower.I2.w), "-",
#                                    round(x$lower.I2.w*100, .round.digits-1)),
#                        "; ", ifelse(is.na(x$upper.I2.w), "-",
#                                     round(x$upper.I2.w*100, .round.digits-1)),"]") %>%
#           ifelse(.=="[-; -]", "-", .)
#
#         # NNT
#         if (is.null(.nnt.cer)){
#           metapsyNNT(g, .model$nnt.cer) %>%
#             round(.round.digits) %>% abs() -> nnt
#         } else {
#           metapsyNNT(g, .nnt.cer) %>%
#             round(.round.digits) %>% abs() -> nnt
#         }
#
#
#         data.frame(variable = y,
#                    group = x$bylevs,
#                    n.comp = x$k.w, g = g, g.ci = g.ci,
#                    i2 = i2, i2.ci = i2.ci,
#                    nnt = nnt,
#                    p = x$pval.Q.b.random %>% scales::pvalue())
#       }) %>% do.call(rbind, .) %>% {rownames(.) = NULL;.} -> summary
#   }
#
#
#   if (class(M)[1] == "rma.mv"){
#
#     stringr::str_replace_all(as.character(M$random[[1]]), "1 \\| ", "")[2] %>%
#       strsplit("/") %>% {.[[1]]} %>% {.[1]} -> study.id
#
#     dat.mv = data.frame(yi = M$yi,
#                         vi = M$vi,
#                         slab = M$slab,
#                         study.id = .model$data[!is.na(.model$data$es), study.id],
#                         es.id = 1:length(M$yi))
#     dat.mv = cbind(dat.mv, .model$data[!is.na(.model$data$es), variables] %>%
#                      purrr::map_dfr(~as.factor(.)))
#
#     purrr::map(as.list(variables), function(x){
#       form = as.formula(paste0("~", x))
#       metafor::rma.mv(yi = yi,
#                       V = vi,
#                       data = dat.mv,
#                       random = M$random[[1]],
#                       test = M$test,
#                       method = "REML",
#                       mods = form) -> res.mv
#       }) -> subgroup.analysis.list
#     names(subgroup.analysis.list) = variables
#
#
#     purrr::map2(subgroup.analysis.list, variables, function(x, y){
#
#       g = c(as.numeric(x$b)[1],
#             as.numeric(x$b)[-1] + as.numeric(x$b)[1]) %>% round(.round.digits)
#       g.lower = {as.numeric(x$b)[1] + x$ci.lb} %>% round(.round.digits)
#       g.upper = {as.numeric(x$b)[1] + x$ci.ub} %>% round(.round.digits)
#       g.ci = paste0("[", g.lower, "; ", g.upper, "]")
#
#       if (is.null(.nnt.cer)){
#         metapsyNNT(g, .model$nnt.cer) %>%
#           round(.round.digits) %>% abs() -> nnt
#       } else {
#         metapsyNNT(g, .nnt.cer) %>%
#           round(.round.digits) %>% abs() -> nnt
#       }
#
#       data.frame(variable = y,
#                  group = levels(dat.mv[[y]]),
#                  n.comp = table(dat.mv[[y]]) %>% as.numeric(),
#                  g = g, g.ci = g.ci, i2 = "-", i2.ci = "-",
#                  nnt = nnt, p = x$QMp %>% scales::pvalue())
#     }) %>% do.call(rbind, .) %>% {rownames(.) = NULL;.} -> summary
#
#   }
#
#
#   # Return
#   returnlist = list(summary = summary,
#                     subgroup.analysis.list = subgroup.analysis.list,
#                     html = .html)
#   class(returnlist) = c("subgroupAnalysis", "list")
#   return(returnlist)
#
# }
#
#
#
#
# psyCtrSubset %>%
#   expandMultiarmTrials() %>%
#   calculateEffectSizes() %>%
#   filterPoolingData(Detect(Cond_spec_trt1 , "cbt|sup"),
#                     Detect(Cond_spec_trt2 , "wl|cau"),
#                     !study %in% "Heckman, 2011") -> data
#
# runMetaAnalysis(data) -> res
#
# subgroupAnalysis(res, Outc_measure, sr_clinician, country) -> result
#
#
# print.subgroupAnalysis = function(x, ...){
#
#   cat("Subgroup analysis results ")
#   cat("---------------------- \n")
#   print(dplyr::as_tibble(x$summary), n = nrow(x$summary))
#
#   if (html == TRUE){
#
#     x$summary %>%
#       magrittr::set_colnames(c("Variable", "Level", "<i>n</i><sub>comp</sub>",
#                                "<i>g</i>", "CI",
#                                "<i>I</i><sup>2</sup>",
#                                "CI", "NNT", "<i>p</i>")) %>%
#       knitr::kable(escape = FALSE) %>%
#       kableExtra::collapse_rows(columns = 1:9, valign = "top") %>%
#       kableExtra::kable_styling(font_size = 8, full_width = FALSE) %>%
#       kableExtra::column_spec(1, bold = TRUE) %>%
#       kableExtra::collapse_rows(columns = 1, valign = "top") %>%
#       print()
#   }
# }
#
#
# plot.subgroupAnalysis = function(x, which = NULL, ...){
#
#   if (is.null(which)){
#     message("- [OK] '", names(x$subgroup.analysis.list)[1], "' used for forest plot.")
#     forest.meta(x$subgroup.analysis.list[[1]])
#   } else {
#     forest.meta(x$subgroup.analysis.list[[which[1]]])
#   }
# }
#
#
#
#
#
#
#
#
#
