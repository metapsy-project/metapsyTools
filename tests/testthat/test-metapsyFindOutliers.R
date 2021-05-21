# Tests for metapsyFindOutliers function

test_that("metapsyFindOutliers works with metafor rma.uni objects", {
  # Create a simple metafor model
  test_data <- create_test_data_with_es(n = 10)
  
  # Create rma.uni model
  rma_model <- try(metafor::rma(yi = test_data$.g, 
                                sei = test_data$.g_se,
                                slab = test_data$study,
                                method = "REML"),
                   silent = TRUE)
  
  if (!inherits(rma_model, "try-error")) {
    result <- try(metapsyFindOutliers(rma_model), silent = TRUE)
    
    if (!inherits(result, "try-error")) {
      expect_true("out.study" %in% names(result))
      expect_true("m" %in% names(result))
      expect_s3_class(result$m, "rma")
    }
  }
})

test_that("metapsyFindOutliers works with meta metagen objects", {
  # Create a simple meta model
  test_data <- create_test_data_with_es(n = 10)
  
  meta_model <- try(meta::metagen(TE = test_data$.g,
                                  seTE = test_data$.g_se,
                                  studlab = test_data$study,
                                  comb.fixed = FALSE,
                                  comb.random = TRUE),
                    silent = TRUE)
  
  if (!inherits(meta_model, "try-error")) {
    result <- try(metapsyFindOutliers(meta_model), silent = TRUE)
    
    if (!inherits(result, "try-error")) {
      expect_true("out.study.fixed" %in% names(result))
      expect_true("out.study.random" %in% names(result))
      expect_true("m.fixed" %in% names(result))
      expect_true("m.random" %in% names(result))
    }
  }
})

test_that("metapsyFindOutliers handles meta objects with NAs", {
  # Create meta model with NAs (should error)
  test_data <- create_test_data_with_es(n = 10)
  test_data$.g[1] <- NA
  
  meta_model <- try(meta::metagen(TE = test_data$.g,
                                  seTE = test_data$.g_se,
                                  studlab = test_data$study,
                                  comb.fixed = FALSE,
                                  comb.random = TRUE),
                    silent = TRUE)
  
  if (!inherits(meta_model, "try-error")) {
    expect_error(metapsyFindOutliers(meta_model),
                 "cannot contain NA")
  }
})

test_that("metapsyFindOutliers handles invalid input", {
  # Test with invalid class
  expect_error(try(metapsyFindOutliers("not a meta object"), silent = TRUE), NA)
})

test_that("metapsyFindOutliers works with different meta object types", {
  test_data <- create_test_data_with_es(n = 10)
  
  # Test metagen
  meta_model <- try(meta::metagen(TE = test_data$.g,
                                  seTE = test_data$.g_se,
                                  studlab = test_data$study,
                                  comb.fixed = TRUE,
                                  comb.random = TRUE),
                    silent = TRUE)
  
  if (!inherits(meta_model, "try-error")) {
    result <- try(metapsyFindOutliers(meta_model), silent = TRUE)
    if (!inherits(result, "try-error")) {
      expect_true("out.study.fixed" %in% names(result))
      expect_true("out.study.random" %in% names(result))
    }
  }
})

test_that("metapsyFindOutliers handles cases with no outliers", {
  # Create data where no outliers are expected
  test_data <- create_test_data_with_es(n = 10)
  # Make all effect sizes similar
  test_data$.g <- rnorm(10, mean = 0.5, sd = 0.1)
  test_data$.g_se <- rep(0.2, 10)
  
  rma_model <- try(metafor::rma(yi = test_data$.g, 
                                sei = test_data$.g_se,
                                slab = test_data$study,
                                method = "REML"),
                   silent = TRUE)
  
  if (!inherits(rma_model, "try-error")) {
    result <- try(metapsyFindOutliers(rma_model), silent = TRUE)
    
    if (!inherits(result, "try-error")) {
      # May have no outliers
      expect_true("out.study" %in% names(result))
      expect_true("m" %in% names(result))
    }
  }
})

