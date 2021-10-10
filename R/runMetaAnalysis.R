psyCtrSubset %>%
  checkDataFormat() %>%
  checkConflicts() %>%
  expandMultiarmTrials() %>%
  calculateEffectSizes() -> data

runMetaAnalysis(data) -> res


# scales pvalue

runMetaAnalysis = function(data,
                           which.run = c("overall", "combined",
                                         "lowest.highest", "outliers",
                                         "influence", "rob", "threelevel"),
                           method.tau = "REML",
                           hakn = TRUE,
                           study.var = "study",
                           arm.var.1 = "Cond_spec_trt1",
                           arm.var.2 = "Cond_spec_trt2",
                           measure.var = "Outc_measure",
                           es.var = "es",
                           se.var = "se",
                           low.rob.filter = "rob > 2",
                           method.tau.ci = "Q-Profile",
                           round.digits = 2,
                           which.outliers = c("general", "combined"),
                           which.influence = c("general", "combined"),
                           which.rob = c("general", "combined"),
                           nnt.cer = 0.2,
                           rho.within.study = 0.5,
                           html = TRUE){


  message("- Pooling the data...")

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  1. Overall Model                                                         #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  data.original = data

  # Create comparison variable
  paste0(data[[study.var]], " (",
         data[[arm.var.1]], " vs. ",
         data[[arm.var.2]], "; ",
         data[[measure.var]], ")") -> data$comparison

  paste0(data[[arm.var.1]], " vs. ",
         data[[arm.var.2]]) -> data$comparison.only

  data[[measure.var]] -> data$instrument


  data = data[!is.na(data[[es.var]]),]
  method.tau.ci = ifelse(method.tau.ci == "Q-Profile", "QP", method.tau.ci)

  mGeneral = meta::metagen(TE = data[[es.var]],
                           seTE = data[[se.var]],
                           studlab = data[[study.var]],
                           data = data,
                           sm = "g",
                           hakn = hakn,
                           method.tau = method.tau,
                           method.tau.ci = method.tau.ci,
                           prediction = TRUE,
                           comb.fixed = ifelse(method.tau == "FE", TRUE, FALSE),
                           comb.random = ifelse(method.tau == "FE", FALSE, TRUE))

  mGeneralRes = with(mGeneral,{

    data.frame(
      k = k,
      g = ifelse(method.tau == "FE", TE.fixed, TE.random) %>% round(round.digits),
      g.ci = paste0("[", ifelse(method.tau == "FE",
                                lower.fixed, lower.random) %>% round(round.digits),"; ",
                    ifelse(method.tau == "FE",
                           upper.fixed, upper.random) %>% round(round.digits), "]"),
      p = ifelse(method.tau == "FE", pval.fixed, pval.random) %>% scales::pvalue(),
      i2 = round(I2*100, 2),
      i2.ci = paste0("[", round(lower.I2*100, 2), "; ", round(upper.I2*100, 2), "]"),
      prediction.ci = paste0("[", round(lower.predict, round.digits), "; ",
                             round(upper.predict, round.digits), "]"),
      nnt = metapsyNNT(ifelse(method.tau == "FE", TE.fixed, TE.random), nnt.cer) %>%
        round(round.digits) %>% abs(),
      excluded = "none"
    )
  })
  rownames(mGeneralRes) = "Overall"


  if ("overall" %in% which.run){
    message("- [OK] Calculated overall effect size")
  }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  2. Lowest Only                                                           #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  multi.study = names(table(data[[study.var]])[table(data[[study.var]]) > 1])

  if (length(multi.study) > 0){
    data %>%
      split(.[[study.var]]) %>%
      purrr::map(function(x){
        tiebreaker = rnorm(nrow(x), sd=1e-10)
        x[[es.var]] + tiebreaker == min(x[[es.var]] + tiebreaker)}) %>%
      do.call(c,.) -> lowest
  } else {
    lowest = NULL
  }

  mLowest = update.meta(mGeneral, exclude = !lowest)
  mLowestRes = with(mLowest,{

    data.frame(
      k = k,
      g = ifelse(method.tau == "FE", TE.fixed, TE.random) %>% round(round.digits),
      g.ci = paste0("[", ifelse(method.tau == "FE",
                                lower.fixed, lower.random) %>% round(round.digits),"; ",
                    ifelse(method.tau == "FE",
                           upper.fixed, upper.random) %>% round(round.digits), "]"),
      p = ifelse(method.tau == "FE", pval.fixed, pval.random) %>% scales::pvalue(),
      i2 = round(I2*100, 2),
      i2.ci = paste0("[", round(lower.I2*100, 2), "; ", round(upper.I2*100, 2), "]"),
      prediction.ci = paste0("[", round(lower.predict, round.digits), "; ",
                             round(upper.predict, round.digits), "]"),
      nnt = metapsyNNT(ifelse(method.tau == "FE", TE.fixed, TE.random), nnt.cer) %>%
        round(round.digits) %>% abs()
    )
  })

  if (is.null(mLowest$exclude[1])){
    mLowestRes$excluded = "none"
  } else {
    mLowestRes$excluded = paste(data$comparison[mLowest$exclude],
                                collapse = ", ")}
  rownames(mLowestRes) = "One ES/study (lowest)"


  if ("lowest.highest" %in% which.run){
    message("- [OK] Calculated effect size using only lowest effect")
  }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  3. Highest Only                                                          #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  multi.study = names(table(data[[study.var]])[table(data[[study.var]]) > 1])

  if (length(multi.study) > 0){
    data %>%
      split(.[[study.var]]) %>%
      purrr::map(function(x){
        tiebreaker = rnorm(nrow(x), sd=1e-10)
        x[[es.var]] + tiebreaker == max(x[[es.var]] + tiebreaker)}) %>%
      do.call(c,.) -> highest
  } else {
    highest = NULL
  }

  mHighest = update.meta(mGeneral, exclude = !highest)
  mHighestRes = with(mHighest,{

    data.frame(
      k = k,
      g = ifelse(method.tau == "FE", TE.fixed, TE.random) %>% round(round.digits),
      g.ci = paste0("[", ifelse(method.tau == "FE",
                                lower.fixed, lower.random) %>% round(round.digits),"; ",
                    ifelse(method.tau == "FE",
                           upper.fixed, upper.random) %>% round(round.digits), "]"),
      p = ifelse(method.tau == "FE", pval.fixed, pval.random) %>% scales::pvalue(),
      i2 = round(I2*100, 2),
      i2.ci = paste0("[", round(lower.I2*100, 2), "; ", round(upper.I2*100, 2), "]"),
      prediction.ci = paste0("[", round(lower.predict, round.digits), "; ",
                             round(upper.predict, round.digits), "]"),
      nnt = metapsyNNT(ifelse(method.tau == "FE", TE.fixed, TE.random), nnt.cer) %>%
        round(round.digits) %>% abs()
    )
  })

  if (is.null(mHighest$exclude[1])){
    mHighestRes$excluded = "none"
  } else {
    mHighestRes$excluded = paste(data$comparison[mHighest$exclude],
                                 collapse = ", ")}
  rownames(mHighestRes) = "One ES/study (highest)"

  if ("lowest.highest" %in% which.run){
    message("- [OK] Calculated effect size using only highest effect")
  }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  4. Combined Effect Sizes                                                 #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  data = metafor::escalc("SMD", yi = data[[es.var]],
                         sei = data[[se.var]], data = data)
  data.comb = metafor::aggregate.escalc(data, cluster = data[[study.var]],
                                        rho = rho.within.study)
  mComb = meta::metagen(TE = data.comb$yi,
                        seTE = sqrt(data.comb$vi),
                        studlab = data.comb[[study.var]],
                        data = data.comb,
                        sm = "g",
                        hakn = hakn,
                        method.tau = method.tau,
                        method.tau.ci = method.tau.ci,
                        prediction = TRUE,
                        comb.fixed = ifelse(method.tau == "FE", TRUE, FALSE),
                        comb.random = ifelse(method.tau == "FE", FALSE, TRUE))

  mCombRes = with(mComb,{

    data.frame(
      k = k,
      g = ifelse(method.tau == "FE", TE.fixed, TE.random) %>% round(round.digits),
      g.ci = paste0("[", ifelse(method.tau == "FE",
                                lower.fixed, lower.random) %>% round(round.digits),"; ",
                    ifelse(method.tau == "FE",
                           upper.fixed, upper.random) %>% round(round.digits), "]"),
      p = ifelse(method.tau == "FE", pval.fixed, pval.random) %>% scales::pvalue(),
      i2 = round(I2*100, 2),
      i2.ci = paste0("[", round(lower.I2*100, 2), "; ", round(upper.I2*100, 2), "]"),
      prediction.ci = paste0("[", round(lower.predict, round.digits), "; ",
                             round(upper.predict, round.digits), "]"),
      nnt = metapsyNNT(ifelse(method.tau == "FE", TE.fixed, TE.random), nnt.cer) %>%
        round(round.digits) %>% abs()
    )
  })
  rownames(mCombRes) = "Combined"
  mCombRes$excluded = paste("combined:", ifelse(length(multi.study) > 0,
                                                paste(multi.study, collapse = ", "), "none"))

  class(data) = "data.frame"

  if ("combined" %in% which.run){
    message("- [OK] Calculated effect size using combined effects (rho=", rho.within.study, ")")
  }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  5. Outliers Removed                                                      #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  if (which.outliers[1] == "general"){
    m.for.outliers = mGeneral
  }

  if (which.outliers[1] == "combined"){
    m.for.outliers = mComb
  }

  if (!(which.outliers[1] %in% c("general", "combined"))){
    stop("'which.outliers' must either be 'general' or 'combined'.")
  }

  if (method.tau == "FE"){
    mOutliers = metapsyFindOutliers(m.for.outliers)$m.fixed
  } else {
    mOutliers = metapsyFindOutliers(m.for.outliers)$m.random
  }

  mOutliersRes = with(mOutliers,{

    data.frame(
      k = k,
      g = ifelse(method.tau == "FE", TE.fixed, TE.random) %>% round(round.digits),
      g.ci = paste0("[", ifelse(method.tau == "FE",
                                lower.fixed, lower.random) %>% round(round.digits),"; ",
                    ifelse(method.tau == "FE",
                           upper.fixed, upper.random) %>% round(round.digits), "]"),
      p = ifelse(method.tau == "FE", pval.fixed, pval.random) %>% scales::pvalue(),
      i2 = round(I2*100, 2),
      i2.ci = paste0("[", round(lower.I2*100, 2), "; ", round(upper.I2*100, 2), "]"),
      prediction.ci = paste0("[", round(lower.predict, round.digits), "; ",
                             round(upper.predict, round.digits), "]"),
      nnt = metapsyNNT(ifelse(method.tau == "FE", TE.fixed, TE.random), nnt.cer) %>%
        round(round.digits) %>% abs()
    )
  })
  rownames(mOutliersRes) = "Outliers removed"
  mOutliersRes$excluded = paste(mOutliers$data$comparison[mOutliers$exclude], collapse = ", ")

  if ("outliers" %in% which.run){
    message("- [OK] Calculated effect size with outliers removed")
  }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  6. Influential Cases                                                     #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  if (which.influence[1] == "general"){
    m.for.influence = mGeneral
  }

  if (which.outliers[1] == "combined"){
    m.for.influence = mComb
  }

  if (!(which.influence[1] %in% c("general", "combined"))){
    stop("'which.influence' must either be 'general' or 'combined'.")
  }

  influenceRes = metapsyInfluenceAnalysis(m.for.influence,
                                          random = ifelse(method.tau == "FE",
                                                          FALSE, TRUE))
  mInfluence = update.meta(m.for.influence,
                           exclude = influenceRes$Data$is.infl == "yes")

  mInfluenceRes = with(mInfluence,{

    data.frame(
      k = k,
      g = ifelse(method.tau == "FE", TE.fixed, TE.random) %>% round(round.digits),
      g.ci = paste0("[", ifelse(method.tau == "FE",
                                lower.fixed, lower.random) %>% round(round.digits),"; ",
                    ifelse(method.tau == "FE",
                           upper.fixed, upper.random) %>% round(round.digits), "]"),
      p = ifelse(method.tau == "FE", pval.fixed, pval.random) %>% scales::pvalue(),
      i2 = round(I2*100, 2),
      i2.ci = paste0("[", round(lower.I2*100, 2), "; ", round(upper.I2*100, 2), "]"),
      prediction.ci = paste0("[", round(lower.predict, round.digits), "; ",
                             round(upper.predict, round.digits), "]"),
      nnt = metapsyNNT(ifelse(method.tau == "FE", TE.fixed, TE.random), nnt.cer) %>%
        round(round.digits) %>% abs()
    )
  })
  rownames(mInfluenceRes) = "Influence Analysis"
  mInfluenceRes$excluded = paste("removed as influential cases:",
                                 ifelse(sum(influenceRes$Data$is.infl == "yes") > 0,
                                        paste(mInfluence$data$comparison[influenceRes$Data$is.infl == "yes"],
                                              collapse = ", "), "none"))

  if ("influence" %in% which.run){
    message("- [OK] Calculated effect size with influential cases removed")
  }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  7. (Low) risk of bias                                                    #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  if (which.rob[1] == "general"){
    m.for.rob = mGeneral
    data.for.rob = mGeneral$data
  }

  if (which.outliers[1] == "combined"){
    m.for.rob = mComb
    data.for.rob = mComb$data
  }

  if (!(which.rob[1] %in% c("general", "combined"))){
    stop("'which.rob' must either be 'general' or 'combined'.")
  }

  if ("rob" %in% which.run){
    robVar = strsplit(low.rob.filter, " ")[[1]][1]
    m.for.rob$data[[robVar]] = as.numeric(m.for.rob$data[[robVar]])
    robFilter = paste0("data.for.rob$", low.rob.filter, " & !is.na(data.for.rob$", robVar, ")")
    robMask = eval(parse(text = robFilter))
    mRob = update.meta(m.for.rob, exclude = !robMask)
  } else {
    robMask = rep(TRUE, nrow(mGeneral$data))
    mRob = update.meta(mGeneral, exclude = !robMask)
  }

  mRobRes = with(mRob,{

    data.frame(
      k = k,
      g = ifelse(method.tau == "FE", TE.fixed, TE.random) %>% round(round.digits),
      g.ci = paste0("[", ifelse(method.tau == "FE",
                                lower.fixed, lower.random) %>% round(round.digits),"; ",
                    ifelse(method.tau == "FE",
                           upper.fixed, upper.random) %>% round(round.digits), "]"),
      p = ifelse(method.tau == "FE", pval.fixed, pval.random) %>% scales::pvalue(),
      i2 = round(I2*100, 2),
      i2.ci = paste0("[", round(lower.I2*100, 2), "; ", round(upper.I2*100, 2), "]"),
      prediction.ci = paste0("[", round(lower.predict, round.digits), "; ",
                             round(upper.predict, round.digits), "]"),
      nnt = metapsyNNT(ifelse(method.tau == "FE", TE.fixed, TE.random), nnt.cer) %>%
        round(round.digits) %>% abs()
    )
  })
  rownames(mRobRes) = paste("Only", low.rob.filter)
  mRobRes$excluded = paste(mRob$data$comparison[!robMask], collapse = ", ")

  if ("rob" %in% which.run){
    message("- [OK] Calculated effect size using only low RoB information")
  }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  8. Three-Level Model                                                     #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  if (method.tau != "REML"){
    message("- [OK] 3L-model tau(s) estimated using 'REML', since '",
            method.tau, "' is not applicable.")
  }

  data$es.id = 1:nrow(data)
  formula = paste0("~ 1 | ", colnames(data[study.var]), "/ es.id")

  mThreeLevel = metafor::rma.mv(yi = data[[es.var]],
                                V = data[[se.var]]*data[[se.var]],
                                slab = data[[study.var]],
                                data = data,
                                random = as.formula(formula),
                                test = ifelse(hakn == TRUE, "t", "z"),
                                method = "REML")

  # Calculate total I2
  W = diag(1/(data[[se.var]]^2))
  X = model.matrix(mThreeLevel)
  P = W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  mThreeLevel$I2 = 100 * sum(mThreeLevel$sigma2) /
    (sum(mThreeLevel$sigma2) + (mThreeLevel$k-mThreeLevel$p)/sum(diag(P)))

  # Calculate I2 per level
  with(mThreeLevel, {
    I2.between.studies = (100 * sigma2 / (sum(sigma2) + (k-p)/sum(diag(P))))[1]
  }) -> mThreeLevel$I2.between.studies

  with(mThreeLevel, {
    I2.between.studies = (100 * sigma2 / (sum(sigma2) + (k-p)/sum(diag(P))))[2]
  }) -> mThreeLevel$I2.within.studies

  # Get tau and I2
  data.frame(tau2 = c(mThreeLevel$sigma2, sum(mThreeLevel$sigma2)),
             i2 = c(mThreeLevel$I2.between.studies, mThreeLevel$I2.within.studies,
                    mThreeLevel$I2)) -> mThreeLevel$variance.components
  mThreeLevel$variance.components$tau2 = round(mThreeLevel$variance.components$tau2, 4)
  mThreeLevel$variance.components$i2 = round(mThreeLevel$variance.components$i2, 1)
  rownames(mThreeLevel$variance.components) = c("Between Studies",
                                                "Within Studies",
                                                "Total")

  mThreeLevelRes = with(mThreeLevel, {
    data.frame(k = k.all,
               g = as.numeric(b[,1]) %>%  round(round.digits),
               g.ci = paste0("[", ci.lb %>% round(round.digits), "; ",
                             ci.ub %>% round(round.digits), "]"),
               p = pval %>% scales::pvalue(),
               i2 = round(I2, 1),
               i2.ci = "-",
               prediction.ci = paste0("[", round(predict(mThreeLevel)$pi.lb, round.digits), "; ",
                                      round(predict(mThreeLevel)$pi.ub, round.digits), "]"),
               nnt = metapsyNNT(as.numeric(b[,1]), nnt.cer) %>%
                 round(round.digits) %>% abs())
  })
  rownames(mThreeLevelRes) = "Three-Level Model"
  mThreeLevelRes$excluded = paste("Number of clusters/studies:", mThreeLevel$s.nlevels[1])

  if ("threelevel" %in% which.run){
    message("- [OK] Calculated effect size using three-level MA model")
  }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #                                                                           #
  #  RETURN                                                                   #
  #                                                                           #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


  # Combine everything
  rbind(mGeneralRes, mLowestRes, mHighestRes, mOutliersRes,
        mInfluenceRes, mRobRes, mCombRes, mThreeLevelRes) -> summary

  # Return
  list(summary = summary,
       model.overall = mGeneral,
       model.combined = mComb,
       model.lowest = mLowest,
       model.highest = mHighest,
       model.outliers = mOutliers,
       model.influence = mInfluence,
       model.rob = mRob,
       model.threelevel = mThreeLevel,
       model.threelevel.var.comp = mThreeLevel$variance.components,
       influence.analysis = influenceRes,
       which.run = which.run,
       data = data.original,
       html = html,
       nnt.cer = nnt.cer) -> returnlist

  class(returnlist) = c("runMetaAnalysis", "list")

  if ("lowest.highest" %in% which.run){
    message("- [OK] Done!")
  }

  return(returnlist)

}



print.runMetaAnalysis = function(x, ...){

  models = list("overall" = 1, "lowest.highest" = c(2,3), "outliers" = 4,
                "influence" = 5, "rob" = 6, "combined" = 7, "threelevel" = 8)

  cat("Model results ")
  cat("---------------------- \n")
  dat = x$summary[unlist(models[x$which.run]),1:8]
  print(dplyr::as_tibble(cbind(model = rownames(dat), dat)))

  if ("threelevel" %in% x$which.run){
    cat("\n")
    cat("Variance components (three-level model) ")
    cat("---------------------- \n")
    print(x$model.threelevel.var.comp)
  }

  if (x$html == TRUE){

    # Add footnote labels
    fn.rows = x$summary$excluded != "none"
    rownames(x$summary)[fn.rows] = paste0(rownames(x$summary[fn.rows,]), "<sup>",
                                 letters[1:nrow(x$summary[fn.rows,])], "</sup>")

    x$summary %>%
      {.$Analysis = rownames(.); rownames(.) = NULL; .} %>%
      dplyr::select(Analysis, dplyr::everything(), -excluded) %>%
      magrittr::set_colnames(c(".", "<i>k</i>", "<i>g</i>", "CI", "<i>p</i>",
                               "<i>I</i><sup>2</sup>",
                               "CI", "PI", "NNT")) %>%
      knitr::kable(escape = FALSE) %>%
      kableExtra::kable_styling(font_size = 8, full_width = FALSE) %>%
      kableExtra::column_spec(1, bold = TRUE, width_min = "13em") %>%
      kableExtra::footnote(general = "Excluded effect sizes/studies:",
                           alphabet = x$summary$excluded[fn.rows]) %>%
      print()
  }
}


plot.runMetaAnalysis = function(x, which = NULL, ...){

  models = list("overall" = "model.overall",
                "lowest.highest" = c("model.lowest", "model.highest"),
                "outliers" = "model.outliers", "influence" = "model.influence",
                "rob" = "model.rob", "combined" = "model.combined",
                "threelevel" = "model.threelevel")

  leftCols = c("studlab", "comparison.only", "instrument", "TE", "seTE")
  leftLabs = c("Study", "Comparison", "Instrument", "g", "S.E.")

  # print forest plot by default
  if (is.null(which)){
    if (models[[x$which.run[1]]][1] != "model.threelevel"){
      message("- [OK] Generating forest plot ('", which.run[1], "' model)")
      meta::forest.meta(x[[models[[x$which.run[1]]][1]]], smlab = " ",
                        leftcols = leftCols, leftlabs = leftLabs,
                        text.random = "RE Model", ...)
    }
    if (models[[x$which.run[1]]][1] == "lowest.highest"){
      message("- [OK] Generating forest plot ('", "highest", "' model)")
      meta::forest.meta(x$model.highest, smlab = " ",
                        leftcols = leftCols, leftlabs = leftLabs,
                        text.random = "RE Model", ...)
    }
    if (models[[x$which.run[1]]][1] == "model.threelevel"){
      metafor::forest.rma(x$model.threelevel, ...)
    }
  } else {

    if (which[1] == "overall"){
      message("- [OK] Generating forest plot ('overall' model)")
      meta::forest.meta(x$model.overall, smlab = " ",
                        leftcols = leftCols, leftlabs = leftLabs,
                        text.random = "RE Model", ...)
    }

    if (which[1] == "lowest.highest"){
      message("- [OK] Generating forest plot ('lowest' model)")
      meta::forest.meta(x$model.lowest, smlab = " ",
                        leftcols = leftCols, leftlabs = leftLabs,
                        text.random = "RE Model", ...)
      message("- [OK] Generating forest plot ('highest' model)")
      meta::forest.meta(x$model.highest, smlab = " ",
                        leftcols = leftCols, leftlabs = leftLabs,
                        text.random = "RE Model", ...)
    }

    if (which[1] == "outliers"){
      message("- [OK] Generating forest plot ('outliers' model)")
      meta::forest.meta(x$model.outliers, smlab = " ",
                        leftcols = leftCols, leftlabs = leftLabs,
                        text.random = "RE Model", ...)
    }

    if (which[1] == "influence"){
      message("- [OK] Generating forest plot ('influence' model)")
      meta::forest.meta(x$model.influence, smlab = " ",
                        leftcols = leftCols, leftlabs = leftLabs,
                        text.random = "RE Model", ...)
    }

    if (which[1] == "rob"){
      message("- [OK] Generating forest plot ('rob' model)")
      meta::forest.meta(x$model.rob, smlab = " ",
                        leftcols = leftCols, leftlabs = leftLabs,
                        text.random = "RE Model", ...)
    }

    if (which[1] == "combined"){
      message("- [OK] Generating forest plot ('combined' model)")
      meta::forest.meta(x$model.combined, smlab = " ",
                        leftcols = leftCols, leftlabs = leftLabs,
                        text.random = "RE Model", ...)
    }

    if (which[1] == "threelevel"){
      message("- [OK] Generating forest plot ('threelevel' model)")
      metafor::forest.rma(x$model.threelevel, ...)
    }

    if (which[1] == "baujat"){
      message("- [OK] Generating baujat plot")
      plot(x$influence.analysis$BaujatPlot)
    }

    if (which[1] == "loo" | which[1] == "loo-es"){
      message("- [OK] Generating leave-one-out forest plot")
      suppressWarnings({plot(x$influence.analysis$ForestEffectSize)})
    }

    if (which[1] == "loo-i2"){
      message("- [OK] Generating leave-one-out forest plot (sorted by I2)")
      suppressWarnings({plot(x$influence.analysis$ForestI2)})
    }

    if (which[1] == "summary"){

      stringr::str_replace_all(x$summary$g.ci, ";|\\]|\\[", "") %>%
        strsplit(" ") %>% purrr::map(~as.numeric(.)) %>% do.call(rbind,.) %>%
        {colnames(.) = c("lower", "upper");.} %>%
        cbind(model = rownames(x$summary), g = x$summary$g,.) %>%
        data.frame() %>%
        dplyr::mutate(g = as.numeric(g),
                      lower = as.numeric(lower),
                      upper = as.numeric(upper)) %>%
        meta::metagen(TE = g, lower = lower, upper = upper, studlab = model,
                      data = .) %>%
        meta::forest.meta(col.square = "lightblue",
                          rightcols = FALSE,
                          overall.hetstat = FALSE,
                          test.overall = FALSE, overall = FALSE,
                          leftlabs = c("Model", "g", "SE"))
  }
  }
}



