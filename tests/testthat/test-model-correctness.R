# Correctness of model results, checked against independent calculations with
# {meta} and {metafor}, plus regression tests for previously detected errors.

get_data <- function(rom = FALSE) {
  d <- suppressMessages(calculateEffectSizes(depressionPsyCtr, calculate.rom = rom))
  d <- filterPoolingData(d, condition_arm2 %in% c("wl", "cau", "other ctr"))
  d[order(d$study), ]
}
quiet <- function(x) suppressWarnings(suppressMessages(x))
ci_num <- function(s) as.numeric(strsplit(gsub("\\[|\\]|;", "", s), " ")[[1]])
all_models <- c("overall", "combined", "lowest.highest", "outliers", "influence",
                "rob", "threelevel", "threelevel.che", "ccrem", "ccrem.che", "waap.wls")

test_that("pooled effects agree with independent calculations", {
  d <- get_data(); dd <- d[!is.na(d$.g), ]; dd$es.id <- seq_len(nrow(dd))
  m <- quiet(runMetaAnalysis(d, which.run = all_models, html = FALSE)); s <- m$summary
  
  f <- metafor::rma(yi = .g, sei = .g_se, data = dd, method = "REML", test = "knha")
  expect_equal(s["Overall", "g"], round(as.numeric(f$b), 2))
  expect_equal(ci_num(s["Overall", "g.ci"]), round(c(f$ci.lb, f$ci.ub), 2))
  expect_equal(s["Overall", "nnt"], 
               round(1/(pnorm(abs(as.numeric(f$b)) + qnorm(.2)) - .2), 2), tolerance = .01)
  
  dd$cl <- paste(dd$study, ifelse(is.na(dd$multi_arm1), "", dd$multi_arm1))
  ag <- metafor::aggregate.escalc(metafor::escalc(yi = .g, sei = .g_se, data = dd, measure = "GEN"), 
                                  cluster = dd$cl, rho = .6)
  fc <- metafor::rma(yi = yi, vi = vi, data = ag, method = "REML", test = "knha")
  expect_equal(s["Combined", "g"], round(as.numeric(fc$b), 2))
  expect_equal(s["Combined", "k"], fc$k)
  
  f3 <- metafor::rma.mv(yi = .g, V = .g_se^2, random = ~ 1 | study/es.id, data = dd, test = "t")
  r3 <- clubSandwich::conf_int(f3, "CR2")
  expect_equal(s["Three-Level Model", "g"], round(as.numeric(f3$b), 2))
  expect_equal(ci_num(s["Three-Level Model", "g.ci"]), round(c(r3$CI_L, r3$CI_U), 2))
  
  vi <- ave(dd$.g_se^2, dd$study)
  V <- metafor::bldiag(lapply(split(vi, dd$study), function(v) {
    x <- matrix(.6*v[1], length(v), length(v)); diag(x) <- v[1]; x }))
  fc3 <- metafor::rma.mv(yi = .g, V = V, random = ~ 1 | study/es.id, data = dd, test = "t")
  rc3 <- clubSandwich::conf_int(fc3, "CR2")
  expect_equal(s["Three-Level Model (CHE)", "g"], round(as.numeric(fc3$b), 2))
  expect_equal(ci_num(s["Three-Level Model (CHE)", "g.ci"]), round(c(rc3$CI_L, rc3$CI_U), 2))
  
  z <- dd$.g/dd$.g_se; prec <- 1/dd$.g_se
  pw <- abs(coef(lm(z ~ 0 + prec)))/(qnorm(.975) + qnorm(.8)) >= dd$.g_se
  expect_equal(s["WAAP-WLS Model", "k"], sum(pw))
  expect_equal(s["WAAP-WLS Model", "g"], round(as.numeric(coef(lm(z[pw] ~ 0 + prec[pw]))), 2))
})

test_that("risk ratios and NNTs agree with independent calculations", {
  d <- get_data(); dr <- d[!is.na(d$.event_arm1), ]
  raw <- quiet(runMetaAnalysis(d, es.measure = "RR", es.type = "raw", which.run = "overall", html = FALSE))
  mh <- meta::metabin(.event_arm1, .totaln_arm1, .event_arm2, .totaln_arm2, data = dr, sm = "RR",
                      method.tau = "REML", method.random.ci = "HK", common = FALSE, random = TRUE)
  expect_equal(raw$summary["Overall", "rr"], round(exp(mh$TE.random), 2))
  rd <- getFromNamespace("update.meta", "meta")(mh, sm = "RD")
  expect_equal(raw$summary["Overall", "nnt"], round(1/abs(rd$TE.random), 2))
  
  e <- metafor::escalc("RR", ai = event_arm1, n1i = totaln_arm1, ci = event_arm2, n2i = totaln_arm2, 
                       data = dr)
  expect_equal(dr$.log_rr, as.numeric(e$yi), tolerance = 1e-6)
  expect_equal(dr$.log_rr_se, sqrt(as.numeric(e$vi)), tolerance = 1e-6)
  
  expect_equal(as.numeric(metapsyNNT(0.5, 0.2)), 1/(pnorm(.5 + qnorm(.2)) - .2))
  expect_equal(as.numeric(metapsyNNT(RR = 1.5, CER = 0.2)), 10)
  # counts are often stored as integers
  expect_equal(as.numeric(metapsyNNT(event.e = 30L, n.e = 100L, event.c = 20L, n.c = 100L)), 10)
})

test_that("effect size signs are as documented", {
  m <- data.frame(mean_arm1 = 10, mean_arm2 = 15, sd_arm1 = 5, sd_arm2 = 5, n_arm1 = 50, n_arm2 = 50)
  expect_equal(g.m.sd(m)$es, -1*(1 - 3/(4*98 - 1)))
  b <- data.frame(event_arm1 = 30, event_arm2 = 10, totaln_arm1 = 50, totaln_arm2 = 50)
  expect_lt(g.binary(b)$es, 0)
  expect_gt(rr.binary(b)$es, 0)
  ch <- data.frame(mean_change_arm1 = .5, mean_change_arm2 = 1.5, sd_change_arm1 = 2.42, 
                   sd_change_arm2 = 2.53, n_change_arm1 = 10, n_change_arm2 = 8)
  expect_equal(g.change.m.sd(ch)$es, -0.3858, tolerance = 1e-3)
})

test_that("event rates use a correct continuity correction", {
  logit <- eventRateLogit(c(0, 5, 10), c(10, 10, 10))
  expect_equal(logit$es, c(log(.5/10.5), 0, log(10.5/.5)))
  expect_equal(logit$se[2], sqrt(1/5 + 1/5))
  
  d <- get_data(); rows <- which(!is.na(d$.event_arm1)); d$.event_arm1[rows[1]] <- 0
  m <- quiet(runMetaAnalysis(d, es.measure = "EER", which.run = "overall", html = FALSE))
  dz <- d[rows, ]
  mp <- meta::metaprop(.event_arm1, .totaln_arm1, data = dz, method = "Inverse", method.tau = "REML",
                       method.random.ci = "HK", sm = "PLOGIT", common = FALSE, incr = .5)
  expect_equal(m$summary["Overall", "eer"], round(plogis(mp$TE.random), 2))
})

test_that("event rate models work without robust variance estimation", {
  d <- get_data()
  m <- quiet(runMetaAnalysis(d, es.measure = "EER", use.rve = FALSE, 
                             which.run = c("overall", "threelevel", "ccrem"), html = FALSE))
  expect_false(any(m$summary$eer.ci == "[-; -]"))
  expect_gt(m$summary["Three-Level Model", "eer"], 0)
  expect_lt(m$summary["Three-Level Model", "eer"], 1)
})

test_that("ratio of means is not back-transformed from rounded values", {
  d <- get_data(rom = TRUE); dr <- d[!is.na(d$.log_rom), ]
  f <- metafor::rma(yi = .log_rom, sei = .log_rom_se, data = dr, test = "knha")
  m <- quiet(runMetaAnalysis(d, es.measure = "ROM", which.run = "overall", html = FALSE, round.digits = 3))
  expect_equal(m$summary["Overall", "rom"], round(exp(as.numeric(f$b)), 3))
  expect_equal(ci_num(m$summary["Overall", "rom.ci"]), round(exp(c(f$ci.lb, f$ci.ub)), 3))
  m1 <- quiet(runMetaAnalysis(d, es.measure = "ROM", which.run = "overall", html = FALSE, round.digits = 1))
  expect_equal(m1$summary["Overall", "rom"], round(exp(as.numeric(f$b)), 1))
})

test_that("ratio of means is undefined for means with opposite signs", {
  x <- data.frame(mean_change_arm1 = -2, mean_change_arm2 = 3, sd_change_arm1 = 2, sd_change_arm2 = 2,
                  n_change_arm1 = 10, n_change_arm2 = 10)
  expect_true(is.na(rom.change.m.sd(x)$es))
  expect_true(is.na(rom.change.m.sd(x)$se))
  y <- data.frame(mean_arm1 = -2, mean_arm2 = -3, sd_arm1 = 2, sd_arm2 = 2, n_arm1 = 10, n_arm2 = 10)
  expect_equal(rom.m.sd(y)$es, log(2/3))
})

test_that("prediction intervals of multilevel models are available with z-tests", {
  d <- get_data()
  m <- quiet(runMetaAnalysis(d, hakn = FALSE, html = FALSE,
                             which.run = c("threelevel", "threelevel.che", "ccrem", "ccrem.che")))
  rows <- c("Three-Level Model", "Three-Level Model (CHE)", "Three-Level CCREM", "Three-Level CCREM (CHE)")
  expect_false(any(m$summary[rows, "prediction.ci"] == "[-; -]"))
})

test_that("a custom study variable can be used in all models", {
  d <- get_data(); names(d)[names(d) == "study"] <- "author"
  ref <- quiet(runMetaAnalysis(get_data(), which.run = c("combined", "threelevel", "threelevel.che", "ccrem"), 
                               html = FALSE))
  m <- quiet(runMetaAnalysis(d, study.var = "author", html = FALSE,
                             which.run = c("combined", "threelevel", "threelevel.che", "ccrem")))
  expect_equal(m$summary[, c("k", "g", "g.ci")], ref$summary[, c("k", "g", "g.ci")])
  sg <- quiet(subgroupAnalysis(m, country, .which.run = "threelevel", .html = FALSE))
  sg.ref <- quiet(subgroupAnalysis(ref, country, .which.run = "threelevel", .html = FALSE))
  expect_equal(sg$summary, sg.ref$summary)
  mc <- quiet(runMetaAnalysis(d, study.var = "author", vcov = "complex", html = FALSE,
                              which.run = c("combined", "threelevel.che")))
  expect_false("combined" %in% unlist(mc$error.model.list))
})

test_that("complex vcov does not need the which.combine.var column", {
  d <- get_data(); d$multi_arm1 <- NULL
  m <- quiet(runMetaAnalysis(d, vcov = "complex", which.run = c("combined", "threelevel.che"), html = FALSE))
  expect_true(all(c("Combined", "Three-Level Model (CHE)") %in% rownames(m$summary)))
})

test_that("subgroup analyses of multilevel models test subgroup differences", {
  d <- get_data(); dd <- d[!is.na(d$.g), ]; dd$es.id <- seq_len(nrow(dd)); dd$cty <- factor(dd$country)
  m <- quiet(runMetaAnalysis(d, which.run = c("threelevel", "threelevel.che", "outliers"), html = FALSE))
  sg <- quiet(subgroupAnalysis(m, country, .which.run = "threelevel", .html = FALSE))
  f <- metafor::rma.mv(yi = .g, V = .g_se^2, random = ~ 1 | study/es.id, mods = ~ cty, data = dd, test = "t")
  expect_equal(sg$summary$p[1], scales::pvalue(f$QMp))
  
  # The correlated sampling errors of the CHE model must be kept
  sg.che <- quiet(subgroupAnalysis(m, country, .which.run = "threelevel.che", .html = FALSE))
  fit.che <- sg.che$subgroup.analysis.list$country
  expect_equal(as.matrix(fit.che$V), as.matrix(m$model.threelevel.che$V), ignore_attr = TRUE)
  
  # Effects excluded in the selected model stay excluded
  sg.o <- quiet(subgroupAnalysis(m, country, .which.run = "outliers", .html = FALSE))
  expect_equal(sum(sg.o$summary$n.comp), m$model.outliers$k)
})

test_that("subgroupAnalysis gives clear errors for unavailable models", {
  d <- get_data()
  m <- quiet(runMetaAnalysis(d, which.run = c("lowest.highest"), html = FALSE))
  expect_error(quiet(subgroupAnalysis(m, country)), "'lowest' or 'highest'")
  expect_no_error(quiet(subgroupAnalysis(m, country, .which.run = "lowest", .html = FALSE)))
  m2 <- quiet(runMetaAnalysis(d, which.run = "waap.wls", html = FALSE))
  expect_error(quiet(subgroupAnalysis(m2, country)), "not available")
})

test_that("low risk of bias filters are evaluated correctly", {
  d <- get_data()
  ref <- quiet(runMetaAnalysis(d, which.run = "rob", low.rob.filter = "rob > 2", html = FALSE))
  k <- ref$model.rob$k
  m <- quiet(runMetaAnalysis(d, which.run = "rob", low.rob.filter = "rob>2", html = FALSE))
  expect_equal(m$model.rob$k, k)
  m2 <- quiet(runMetaAnalysis(d, which.run = "rob", low.rob.filter = "rob %in% c(3,4)", html = FALSE))
  expect_equal(m2$model.rob$k, k)
  expect_error(quiet(runMetaAnalysis(d, which.run = "rob", low.rob.filter = "nonexistent > 2", html = FALSE)),
               "not found in data")
  mask <- evalRobFilter("rob > 2 & !is.na(year)", d)
  expect_type(mask, "logical")
  expect_false(anyNA(mask))
})

test_that("input columns are not overwritten in the returned data", {
  d <- get_data(rom = TRUE)
  m <- quiet(runMetaAnalysis(d, es.measure = "ROM", which.run = "overall", html = FALSE))
  keep <- !is.na(d$.log_rom)
  expect_equal(m$data$.g, d$.g[keep])
  e <- quiet(runMetaAnalysis(d, es.measure = "EER", which.run = "overall", html = FALSE))
  expect_equal(e$data$.g, d$.g[!is.na(d$.event_arm1)])
})

test_that("plot and summary methods work for every single model", {
  d <- get_data(); pdf(NULL); on.exit(dev.off())
  for (w in all_models){
    m <- quiet(runMetaAnalysis(d, which.run = w, html = FALSE))
    expect_no_error(quiet(plot(m)))
    expect_no_error(quiet(capture.output(summary(m))))
  }
  m <- quiet(runMetaAnalysis(d, which.run = c("combined", "waap.wls", "ccrem", "ccrem.che"), html = FALSE))
  expect_identical(m$which.waap.wls, "combined")
  for (w in c("waap.wls", "ccrem", "ccrem.che")) expect_no_error(quiet(plot(m, w)))
  m.none <- quiet(runMetaAnalysis(d, low.rob.filter = "rob > 100", which.run = c("overall", "rob"), html = FALSE))
  expect_error(quiet(plot(m.none, "rob")), "was not fitted")
})

test_that("empirical Bayes estimates work with raw event data", {
  d <- get_data()
  m <- quiet(runMetaAnalysis(d, es.measure = "RR", es.type = "raw", which.run = "overall", html = FALSE))
  expect_s3_class(quiet(eb(m)), "list.rma")
})

test_that("summary states the assumed correlation of the combined model", {
  d <- get_data()
  m <- quiet(runMetaAnalysis(d, which.run = "combined", rho.within.study = .7, html = FALSE))
  out <- paste(capture.output(quiet(summary(m, forest = FALSE))), collapse = " ")
  expect_match(out, "rho = 0.7")
  expect_match(out, "arms")
})

test_that("printed tables are complete without a color console", {
  d <- get_data()
  m <- quiet(runMetaAnalysis(d, which.run = c("overall", "threelevel"), html = FALSE))
  out <- capture.output(print(m))
  expect_true(any(grepl("^Overall", out)))
  expect_true(any(grepl("^Three-Level Model", out)))
  expect_true(any(grepl(m$summary["Overall", "g.ci"], out, fixed = TRUE)))
  expect_true(any(grepl("^Between Studies", out)))
})

test_that("correctPublicationBias agrees with independent calculations", {
  d <- get_data(); dd <- d[!is.na(d$.g), ]
  m <- quiet(runMetaAnalysis(d, which.run = "overall", html = FALSE))
  p <- quiet(correctPublicationBias(m))$correctPublicationBias$summary
  f <- metafor::rma(yi = .g, sei = .g_se, data = dd, method = "REML")
  sm <- metafor::selmodel(f, type = "stepfun", alternative = "less", steps = .05)
  expect_equal(p["Selection model", "g"], round(as.numeric(sm$b), 2))
  lim <- metasens::limitmeta(m$model.overall)
  expect_equal(p["Limit meta-analysis", "g"], round(lim$TE.adjust, 2))
  tf <- meta::trimfill(m$model.overall)
  expect_equal(p["Trim-and-fill method", "g"], round(tf$TE.random, 2))
})

test_that("selection model direction follows the coding of the effect sizes", {
  d <- get_data()
  m <- quiet(runMetaAnalysis(d, which.run = "overall", html = FALSE))
  neg <- quiet(correctPublicationBias(m))$correctPublicationBias$summary
  pos.d <- d; pos.d$.g <- -pos.d$.g
  m.pos <- quiet(runMetaAnalysis(pos.d, which.run = "overall", html = FALSE))
  pos <- quiet(correctPublicationBias(m.pos))$correctPublicationBias$summary
  # Mirrored effect sizes lead to a mirrored estimate
  expect_equal(pos["Selection model", "g"], -neg["Selection model", "g"])
  # Explicitly contradicting the sign of the pooled effect gives a warning
  expect_warning(suppressMessages(correctPublicationBias(m.pos, lower.is.better = TRUE)),
                 "wrong direction")
  expect_no_warning(suppressMessages(correctPublicationBias(m.pos, lower.is.better = FALSE)))
})

test_that("correctPublicationBias respects the type of model", {
  d <- get_data(); dd <- d[!is.na(d$.g), ]
  m <- quiet(runMetaAnalysis(d, method.tau = "FE", which.run = "overall", html = FALSE))
  p <- quiet(correctPublicationBias(m))
  s <- p$correctPublicationBias$summary
  expect_equal(p$correctPublicationBias$model.selection$method, "FE")
  f <- metafor::selmodel(metafor::rma(yi = .g, sei = .g_se, data = dd, method = "FE"), 
                         type = "stepfun", alternative = "less", steps = .05)
  expect_equal(s["Selection model", "g"], round(as.numeric(f$b), 2))
  tf <- meta::trimfill(m$model.overall, common = TRUE, random = FALSE)
  expect_equal(s["Trim-and-fill method", "g"], round(tf$TE.common, 2))
})

test_that("'lowest' and 'highest' can be used in other functions", {
  d <- get_data()
  rd <- list(domains = c("sg", "ac"), domain.names = c("Sequence", "Allocation"),
             categories = c("0", "1", "sr"), symbols = c("-", "+", "s"), 
             colors = c("red", "green", "yellow"))
  m <- quiet(runMetaAnalysis(d, which.run = "lowest.highest", rob.data = rd, html = FALSE))
  expect_s3_class(m$model.lowest$rob, "rob")
  p <- quiet(correctPublicationBias(m, which.run = "highest"))
  expect_equal(p$correctPublicationBias$which.run, "highest")
  expect_error(quiet(correctPublicationBias(m)), "'lowest' or 'highest'")
  pdf(NULL); on.exit(dev.off())
  expect_no_error(quiet(createRobSummary(m, name.low = "1", name.high = "0", name.unclear = "sr",
                                         which.run = "lowest")))
})

test_that("limit meta-analysis of event rates shows the adjusted estimate", {
  d <- get_data()
  m <- quiet(runMetaAnalysis(d, es.measure = "EER", which.run = "overall", html = FALSE))
  p <- quiet(correctPublicationBias(m))
  s <- p$correctPublicationBias$summary
  expect_equal(as.numeric(s["Limit meta-analysis", "eer"]), 
               round(plogis(p$correctPublicationBias$model.limitmeta$TE.adjust), 2))
  expect_match(s["Limit meta-analysis", "excluded"], "G-squared")
  expect_equal(s["Selection model", "eer"], "-")
  pdf(NULL); on.exit(dev.off())
  expect_error(quiet(plot(p, "selection")), "could not be calculated")
})

test_that("NNTs of pre-calculated risk ratios are available after correction", {
  d <- get_data()
  m <- quiet(runMetaAnalysis(d, es.measure = "RR", which.run = "overall", html = FALSE))
  p <- quiet(correctPublicationBias(m))$correctPublicationBias$summary
  expect_false(p["Selection model", "nnt"] == "-")
  expect_equal(as.numeric(p["Selection model", "nnt"]),
               round(as.numeric(metapsyNNT(RR = as.numeric(p["Selection model", "rr"]), CER = .2)), 2),
               tolerance = .1)
})

test_that("bootstrap of the G2 confidence interval is not silently empty", {
  d <- get_data()
  m <- quiet(runMetaAnalysis(d, which.run = "overall", i2.ci.boot = TRUE, nsim.boot = 30, html = FALSE))
  set.seed(1)
  invisible(capture.output(p <- quiet(correctPublicationBias(m))))
  ci <- p$correctPublicationBias$model.limitmeta$G.squared.ci
  expect_length(ci, 2)
  expect_false(anyNA(ci))
  expect_false(p$correctPublicationBias$summary["Limit meta-analysis", "i2.ci"] %in% c("-", "[NA; NA]"))
})

test_that("imputeResponse counts responders when higher values are better", {
  r <- quiet(imputeResponse(m.trt.pre = 10, m.trt.post = 18, sd.trt.post = 4, n.trt = 100,
                            m.ctr.pre = 10, m.ctr.post = 12, sd.ctr.post = 4, n.ctr = 100,
                            lower.is.better = FALSE))
  # Responders: score above a doubling of the baseline value
  expect_equal(r$trtResponder, round(100*(1 - pnorm((20 - 18)/4))))
  expect_equal(r$ctrResponder, round(100*(1 - pnorm((20 - 12)/4))))
  expect_gt(r$logRR, 0)
  expect_error(quiet(imputeResponse(10, 8, 4, 100, sd.trt.pre = 4, rho = 2, method = "rci")), "between 0 and 1")
  rci <- quiet(imputeResponse(m.trt.pre = 20, m.trt.post = 10, sd.trt.post = 4, n.trt = 100,
                              sd.trt.pre = 4, rho = .8, method = "rci"))
  expect_equal(names(rci), c("trtResponder", "nTrt"))
})

test_that("calculateEffectSizes handles signs, switched arms and missing values", {
  d <- depressionPsyCtr[1:6, ]
  d$flip <- c(TRUE, FALSE, NA, TRUE, FALSE, NA)
  a <- quiet(calculateEffectSizes(d))
  b <- quiet(calculateEffectSizes(d, change.sign = "flip"))
  expect_equal(b$.g, a$.g * c(-1, 1, 1, -1, 1, 1))
  s <- quiet(calculateEffectSizes(depressionPsyCtr[1:6, ], calculate.rom = TRUE, include.switched.arms = TRUE))
  orig <- s[!grepl("switched", s$.id), ]; sw <- s[grepl("switched", s$.id), ]
  expect_equal(sw$.log_rom, -orig$.log_rom)
  expect_equal(sw$.g, -orig$.g)
  changed <- a; changed$.g[1] <- 99
  expect_message(calculateEffectSizes(changed), "replaced")
  x <- d; x$mean_arm1[1] <- NA
  expect_message(calculateEffectSizes(x), "No effect size could be calculated")
  y <- d; y$mean_change_arm1 <- NULL; y$change_m_arm1 <- NA
  expect_message(calculateEffectSizes(y), "change_m_arm1")
})

test_that("filterPoolingData can drop rows without effect size by default", {
  d <- suppressMessages(calculateEffectSizes(depressionPsyCtr[1:6, ])); d$.g[1] <- NA
  expect_equal(nrow(filterPoolingData(d, .filter.missing.rows = TRUE)), 5)
})

test_that("proportionMID is correct for different models", {
  d <- get_data(rom = TRUE)
  m <- quiet(runMetaAnalysis(d, which.run = c("overall", "threelevel", "waap.wls", "ccrem"), html = FALSE))
  pm <- quiet(proportionMID(m, mid = -.24))
  q <- pm$qhat
  expect_true(all(c("Overall", "Three-Level Model", "Three-Level CCREM") %in% rownames(q)))
  o <- m$model.overall
  expect_equal(as.numeric(q["Overall", 3]), round(100*pnorm((-.24 - o$TE.random)/o$tau), 2), tolerance = .01)
  t3 <- m$model.threelevel
  expect_equal(as.numeric(q["Three-Level Model", 3]),
               round(100*pnorm((-.24 - as.numeric(t3$b))/sqrt(sum(t3$sigma2))), 2), tolerance = .01)
  
  fe <- quiet(runMetaAnalysis(d, method.tau = "FE", which.run = "overall", html = FALSE))
  expect_equal(quiet(proportionMID(fe, mid = -.24))$qhat$g, round(fe$model.overall$TE.common, 2))
  
  rom <- quiet(runMetaAnalysis(d, es.measure = "ROM", which.run = "overall", html = FALSE))
  pr <- quiet(proportionMID(rom, mid = .8))
  or <- rom$model.overall
  expect_equal(as.numeric(pr$qhat["Overall", 3]), 
               round(100*pnorm((log(.8) - or$TE.random)/or$tau), 2), tolerance = .01)
  eer <- quiet(runMetaAnalysis(d, es.measure = "EER", which.run = "overall", html = FALSE))
  expect_error(proportionMID(eer, mid = .5), "not available")
})

test_that("average number of treatments accounts for response rate decay", {
  s <- simulateTreatmentCycles(.5, .2)
  resp <- diff(c(0, s$data$cum.response))
  expect_equal(s$avg.no.treatments, sum(s$data$cycle*resp)/sum(resp))
  s0 <- simulateTreatmentCycles(.5, 0)
  # Close to the mean of a geometric distribution with p = 0.5
  expect_equal(s0$avg.no.treatments, 2, tolerance = .05)
  s2 <- simulateTreatmentCycles(c(.5, .3), .1)
  expect_equal(s2$data$treated.per.100, round(s2$data$prop.treated*100))
})

test_that("risk of bias summary works with one domain", {
  d <- get_data(); pdf(NULL); on.exit(dev.off())
  rd <- list(domains = "sg", domain.names = "Sequence", categories = c("0", "1", "sr"), 
             symbols = c("-", "+", "s"), colors = c("red", "green", "yellow"))
  m <- quiet(runMetaAnalysis(d, which.run = "overall", rob.data = rd, html = FALSE))
  expect_no_error(quiet(createRobSummary(m, name.low = "1", name.high = "0", name.unclear = "sr")))
  expect_error(quiet(createRobSummary(m, name.low = "x", name.high = "y", name.unclear = "z")), 
               "specified ratings")
})
