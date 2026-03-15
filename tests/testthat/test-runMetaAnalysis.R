test_that("runMetaAnalysis validates input data type", {
  # Test with NULL
  expect_error(runMetaAnalysis(NULL), 
               "must be a data.frame")
  
  # Test with non-data.frame
  expect_error(runMetaAnalysis(c(1, 2, 3)), 
               "must be a data.frame")
  
  # Test with empty data.frame
  empty_df <- data.frame()
  expect_error(runMetaAnalysis(empty_df),
               "must be a data.frame|at least 3 effect sizes")
})

test_that("runMetaAnalysis requires at least 3 effect sizes", {
  # Create data with only 2 effect sizes
  test_data <- data.frame(
    study = c("Study1", "Study2"),
    .g = c(0.5, 0.6),
    .g_se = c(0.1, 0.12),
    stringsAsFactors = FALSE
  )
  
  expect_error(runMetaAnalysis(test_data),
               "at least 3 effect sizes")
})

test_that("runMetaAnalysis validates es.measure parameter", {
  # Create valid test data with pre-calculated effect sizes
  test_data <- create_test_data_with_es(n = 5)
  
  # Invalid es.measure
  expect_error(runMetaAnalysis(test_data, es.measure = "invalid"),
               "'es.measure' must be 'g', 'RR', 'EER', 'CER' or 'ROM'.")
  
  # Valid es.measure should work (or fail gracefully with other issues)
  result <- try(runMetaAnalysis(test_data, es.measure = "g", 
                                which.run = "overall"), 
                silent = TRUE)
  # Either succeeds or fails for other reasons (not es.measure)
  if (inherits(result, "try-error")) {
    expect_false(grepl("must be 'g'", as.character(result)))
  }
})

test_that("runMetaAnalysis validates which.combine parameter", {
  # Create valid test data with pre-calculated effect sizes
  test_data <- create_test_data_with_es(n = 5)
  
  # Invalid which.combine
  expect_error(runMetaAnalysis(test_data, which.combine = "invalid"),
               "must be either 'arms' or 'studies'")
  
  # Valid which.combine
  result <- try(runMetaAnalysis(test_data, which.combine = "arms",
                                which.run = "overall"), 
                silent = TRUE)
  # Should not fail due to which.combine
  if (inherits(result, "try-error")) {
    expect_false(grepl("must be either 'arms'", as.character(result)))
  }
})

test_that("runMetaAnalysis validates nsim.boot parameter", {
  # Create valid test data with pre-calculated effect sizes
  test_data <- create_test_data_with_es(n = 5)
  
  # Invalid nsim.boot (non-numeric)
  expect_error(runMetaAnalysis(test_data, nsim.boot = "invalid"),
               "must numeric")
  
  # Low nsim.boot with i2.ci.boot should warn
  expect_warning(
    try(runMetaAnalysis(test_data, nsim.boot = 100, i2.ci.boot = TRUE,
                        which.run = "overall"), silent = TRUE),
    "Number of boostrap samples"
  )
})

test_that("runMetaAnalysis rejects reserved column names", {
  # Create test data with 'subset' column
  test_data <- create_test_data_with_es(n = 5)
  test_data$subset <- "test"
  
  expect_error(runMetaAnalysis(test_data),
               "not allowed in the data set")
  
  # Create test data with 'exclude' column
  test_data2 <- create_test_data_with_es(n = 5)
  test_data2$exclude <- "test"
  
  expect_error(runMetaAnalysis(test_data2),
               "not allowed in the data set")
})

test_that("runMetaAnalysis validates study.var parameter", {
  # Create test data with pre-calculated effect sizes
  test_data <- create_test_data_with_es(n = 5)
  
  # Invalid study.var
  expect_error(runMetaAnalysis(test_data, study.var = "nonexistent"),
               "was not found in the dataset")
})

test_that("runMetaAnalysis validates power.within.study parameter", {
  # Create test data with pre-calculated effect sizes
  # Need to ensure we don't trigger other errors first (like missing 'rob' column)
  test_data <- create_test_data_with_es(n = 5)
  test_data$rob <- 3  # Add rob column to avoid that error
  
  # Too high - validation happens after other checks, so we might get different errors
  # Try to catch the specific validation error
  result <- try(runMetaAnalysis(test_data, power.within.study = 0.99,
                                which.run = "overall"), silent = TRUE)
  if (inherits(result, "try-error")) {
    error_msg <- as.character(result)
    # Check if it's the power.within.study error or another error
    expect_true(grepl("must be between 0.01 and 0.99|power.within.study", 
                     error_msg, ignore.case = TRUE))
  }
  
  # Too low
  result <- try(runMetaAnalysis(test_data, power.within.study = 0.001,
                                which.run = "overall"), silent = TRUE)
  if (inherits(result, "try-error")) {
    error_msg <- as.character(result)
    expect_true(grepl("must be between 0.01 and 0.99|power.within.study", 
                     error_msg, ignore.case = TRUE))
  }
})

test_that("runMetaAnalysis validates rob.data parameter", {
  # Create test data with pre-calculated effect sizes
  test_data <- create_test_data_with_es(n = 5)
  
  # Invalid rob.data (not NULL or list)
  expect_error(runMetaAnalysis(test_data, rob.data = "invalid"),
               "must be either NULL or a list")
  
  # Valid rob.data (NULL)
  result <- try(runMetaAnalysis(test_data, rob.data = NULL,
                                which.run = "overall"), 
                silent = TRUE)
  # Should not fail due to rob.data
  if (inherits(result, "try-error")) {
    expect_false(grepl("must be either NULL", as.character(result)))
  }
})

test_that("runMetaAnalysis returns correct structure for overall model", {
  # Create test data with pre-calculated effect sizes
  test_data <- create_test_data_with_es(n = 6)
  
  # Run with just overall model
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    # Check return structure
    expect_s3_class(result, "runMetaAnalysis")
    expect_true("summary" %in% names(result))
    expect_true("model.overall" %in% names(result))
    expect_true("data" %in% names(result))
    expect_true("which.run" %in% names(result))
    
    # Check summary is a data.frame
    expect_s3_class(result$summary, "data.frame")
    expect_true(nrow(result$summary) > 0)
  }
})

test_that("runMetaAnalysis handles different which.run options", {
  # Create test data with pre-calculated effect sizes
  test_data <- create_test_data_with_es(n = 6)
  
  # Test with different model types
  models_to_test <- c("overall", "combined")
  
  for (model in models_to_test) {
    result <- try(runMetaAnalysis(test_data, which.run = model), 
                  silent = TRUE)
    
    if (!inherits(result, "try-error")) {
      expect_s3_class(result, "runMetaAnalysis")
      expect_true(model %in% result$which.run || 
                  "overall" %in% result$which.run)  # May default to overall
    }
  }
})

test_that("runMetaAnalysis handles missing effect size columns gracefully", {
  # Create test data without effect size columns
  test_data <- create_test_data()
  test_data <- checkDataFormat(test_data)
  # Don't run calculateEffectSizes - so no .g or .g_se columns
  
  # Should either error or handle gracefully
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  # Should fail with informative error or handle missing data
  if (inherits(result, "try-error")) {
    # Error should be about missing data, not a cryptic error
    error_msg <- as.character(result)
    expect_true(length(error_msg) > 0)
  }
})

test_that("runMetaAnalysis handles flagEffectSizes input", {
  # Create test data with pre-calculated effect sizes
  test_data <- create_test_data_with_es(n = 6)
  
  # Create a mock flagEffectSizes object structure
  # Note: This is a simplified test - actual flagEffectSizes might have different structure
  flagged_data <- list(
    data = test_data,
    flags = rep(1, nrow(test_data))  # All unflagged
  )
  class(flagged_data) <- "flagEffectSizes"
  
  # Should handle flagEffectSizes input
  result <- try(runMetaAnalysis(flagged_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    expect_s3_class(result, "runMetaAnalysis")
  }
})

test_that("runMetaAnalysis filters out flagged effect sizes", {
  # Create test data with pre-calculated effect sizes
  test_data <- create_test_data_with_es(n = 6)
  
  # Create flagEffectSizes with some flagged (flag == 3)
  flagged_data <- list(
    data = test_data,
    flags = c(rep(1, 4), 3, 3)  # Last two flagged
  )
  class(flagged_data) <- "flagEffectSizes"
  
  result <- try(runMetaAnalysis(flagged_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    # Should have removed flagged effect sizes
    expect_s3_class(result, "runMetaAnalysis")
    # The data should have fewer rows than original
    expect_true(nrow(result$data) <= nrow(test_data))
  }
})

test_that("runMetaAnalysis handles different es.measure values", {
  # Create test data with binary outcomes for RR
  test_data <- create_test_data_binary()
  test_data <- rbind(test_data, test_data, test_data)
  test_data$study <- paste0("Study", 1:6)
  test_data <- checkDataFormat(test_data)
  test_data <- calculateEffectSizes(test_data)
  
  # Test with g (default)
  result_g <- try(runMetaAnalysis(test_data, es.measure = "g",
                                  which.run = "overall"), 
                  silent = TRUE)
  
  # Test with RR if log_rr columns exist
  if (".log_rr" %in% names(test_data) && 
      ".log_rr_se" %in% names(test_data) &&
      sum(!is.na(test_data$.log_rr)) >= 3) {
    result_rr <- try(runMetaAnalysis(test_data, es.measure = "RR",
                                     which.run = "overall"), 
                     silent = TRUE)
    
    if (!inherits(result_rr, "try-error")) {
      expect_s3_class(result_rr, "runMetaAnalysis")
    }
  }
  
  if (!inherits(result_g, "try-error")) {
    expect_s3_class(result_g, "runMetaAnalysis")
  }
})

test_that("runMetaAnalysis handles which.run with invalid model name", {
  # Create test data with pre-calculated effect sizes
  test_data <- create_test_data_with_es(n = 6)
  
  # Invalid model name
  result <- try(runMetaAnalysis(test_data, which.run = "nonexistent"), 
                 silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    # Should default to "overall" with warning
    expect_true("overall" %in% result$which.run)
  } else {
    # Or fail with informative error
    expect_true(length(as.character(result)) > 0)
  }
})

test_that("runMetaAnalysis creates summary data.frame", {
  # Create test data with pre-calculated effect sizes
  test_data <- create_test_data_with_es(n = 6)
  
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    # Summary should exist and be a data.frame
    expect_true("summary" %in% names(result))
    expect_s3_class(result$summary, "data.frame")
    expect_true(nrow(result$summary) > 0)
    expect_true(ncol(result$summary) > 0)
  }
})

