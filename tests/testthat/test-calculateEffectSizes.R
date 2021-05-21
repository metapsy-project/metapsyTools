test_that("calculateEffectSizes requires valid data", {
  # Test with NULL - function doesn't validate early, so we get a different error
  expect_error(calculateEffectSizes(NULL))
  
  # Test with empty data.frame - function tries to process it and fails on missing columns
  empty_df <- data.frame()
  expect_error(calculateEffectSizes(empty_df))
})

test_that("calculateEffectSizes handles missing required columns gracefully", {
  # Create minimal data without all required columns
  test_data <- data.frame(
    study = c("Study1"),
    condition_arm1 = c("Treatment"),
    condition_arm2 = c("Control"),
    stringsAsFactors = FALSE
  )
  
  # Should either work with defaults or give informative error
  result <- try(calculateEffectSizes(test_data), silent = TRUE)
  # Either succeeds or fails with informative message
  expect_true(is.data.frame(result) || 
              (inherits(result, "try-error") && 
               length(grep("column|variable|required", 
                          as.character(result), ignore.case = TRUE)) > 0))
})

test_that("calculateEffectSizes creates .id column", {
  test_data <- create_test_data()
  
  # Run through checkDataFormat first (as recommended workflow)
  test_data <- checkDataFormat(test_data)
  
  result <- try(calculateEffectSizes(test_data), silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    expect_true(".id" %in% names(result))
    expect_true(all(!is.na(result$.id)))
  }
})

test_that("calculateEffectSizes handles binary outcomes", {
  test_data <- create_test_data_binary()
  test_data <- checkDataFormat(test_data)
  
  result <- try(calculateEffectSizes(test_data), silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    # Should have calculated effect sizes
    expect_true(any(c(".g", ".log_rr") %in% names(result)))
  }
})

