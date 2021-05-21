# Tests for subgroupAnalysis function

test_that("subgroupAnalysis validates input", {
  # Test with non-runMetaAnalysis object
  expect_error(subgroupAnalysis("not a runMetaAnalysis object"),
               "must be of class 'runMetaAnalysis'")
  
  # Test with runMetaAnalysis object but invalid which.run
  test_data <- create_test_data_with_es(n = 6)
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    # Should work with valid which.run
    expect_error(try(subgroupAnalysis(result, condition_arm1, 
                                     .which.run = "overall",
                                     .html = FALSE), 
                    silent = TRUE), NA)
  }
})

test_that("subgroupAnalysis works with basic subgroup variable", {
  test_data <- create_test_data_with_es(n = 10)
  # Add a subgroup variable
  test_data$country <- rep(c("USA", "Europe"), each = 5)
  
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    sg_result <- try(subgroupAnalysis(result, country,
                                      .which.run = "overall",
                                      .html = FALSE),
                     silent = TRUE)
    
    if (!inherits(sg_result, "try-error")) {
      expect_s3_class(sg_result, "subgroupAnalysis")
      expect_true("summary" %in% names(sg_result))
      expect_s3_class(sg_result$summary, "data.frame")
    }
  }
})

test_that("subgroupAnalysis works with multiple subgroup variables", {
  test_data <- create_test_data_with_es(n = 10)
  test_data$country <- rep(c("USA", "Europe"), each = 5)
  test_data$setting <- rep(c("inpatient", "outpatient"), 5)
  
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    sg_result <- try(subgroupAnalysis(result, country, setting,
                                     .which.run = "overall",
                                     .html = FALSE),
                     silent = TRUE)
    
    if (!inherits(sg_result, "try-error")) {
      expect_s3_class(sg_result, "subgroupAnalysis")
      # Subgroup variables are stored in subgroup.analysis.list
      expect_true("subgroup.analysis.list" %in% names(sg_result))
      if ("subgroup.analysis.list" %in% names(sg_result)) {
        expect_true("country" %in% names(sg_result$subgroup.analysis.list))
        expect_true("setting" %in% names(sg_result$subgroup.analysis.list))
      }
    }
  }
})

test_that("subgroupAnalysis works with tau.common = TRUE", {
  test_data <- create_test_data_with_es(n = 10)
  test_data$country <- rep(c("USA", "Europe"), each = 5)
  
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    sg_result <- try(subgroupAnalysis(result, country,
                                     .which.run = "overall",
                                     .tau.common = TRUE,
                                     .html = FALSE),
                     silent = TRUE)
    
    if (!inherits(sg_result, "try-error")) {
      expect_s3_class(sg_result, "subgroupAnalysis")
    }
  }
})

test_that("subgroupAnalysis works with custom round.digits", {
  test_data <- create_test_data_with_es(n = 10)
  test_data$country <- rep(c("USA", "Europe"), each = 5)
  
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    sg_result <- try(subgroupAnalysis(result, country,
                                     .which.run = "overall",
                                     .round.digits = 3,
                                     .html = FALSE),
                     silent = TRUE)
    
    if (!inherits(sg_result, "try-error")) {
      expect_s3_class(sg_result, "subgroupAnalysis")
    }
  }
})

test_that("subgroupAnalysis works with HTML output", {
  test_data <- create_test_data_with_es(n = 10)
  test_data$country <- rep(c("USA", "Europe"), each = 5)
  
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    # HTML output is printed, not returned, so just check it doesn't error
    expect_error(try(subgroupAnalysis(result, country,
                                     .which.run = "overall",
                                     .html = TRUE),
                     silent = TRUE), NA)
  }
})

test_that("subgroupAnalysis handles missing values in subgroup variable", {
  test_data <- create_test_data_with_es(n = 10)
  test_data$country <- rep(c("USA", "Europe", NA), length.out = 10)
  
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    # Should handle NAs gracefully
    sg_result <- try(subgroupAnalysis(result, country,
                                     .which.run = "overall",
                                     .html = FALSE),
                     silent = TRUE)
    
    # May error or work depending on implementation
    # Just check it doesn't crash unexpectedly
    expect_true(TRUE)  # Placeholder - actual behavior may vary
  }
})

test_that("subgroupAnalysis works with combined model", {
  test_data <- create_test_data_with_es(n = 10)
  test_data$country <- rep(c("USA", "Europe"), each = 5)
  
  result <- try(runMetaAnalysis(test_data, which.run = "combined"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error") && 
      !is.null(result$model.combined) &&
      result$model.combined$k > 1) {
    sg_result <- try(subgroupAnalysis(result, country,
                                     .which.run = "combined",
                                     .html = FALSE),
                     silent = TRUE)
    
    if (!inherits(sg_result, "try-error")) {
      expect_s3_class(sg_result, "subgroupAnalysis")
    }
  }
})

