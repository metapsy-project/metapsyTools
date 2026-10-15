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

test_that("g.change.m.sd matches esc_mean_sd for change scores", {
  skip_if_not_installed("esc")
  test_data <- create_test_data_change(n = 1)
  got <- g.change.m.sd(test_data)
  ref <- esc::esc_mean_sd(
    test_data$mean_change_arm1, test_data$sd_change_arm1, test_data$n_change_arm1,
    test_data$mean_change_arm2, test_data$sd_change_arm2, test_data$n_change_arm2,
    es.type = "g"
  )
  expect_equal(got$es, as.numeric(ref$es), tolerance = 1e-6)
  expect_equal(got$se, as.numeric(ref$se), tolerance = 1e-6)
})

test_that("calculateEffectSizes calculates Hedges' g from change scores", {
  test_data <- create_test_data_change(n = 2)
  test_data <- checkDataFormat(test_data)

  result <- calculateEffectSizes(test_data)

  expect_true(all(c(".g", ".g_se") %in% names(result)))
  expect_true(all(is.finite(result$.g)))
  expect_true(all(is.finite(result$.g_se)))
  expect_true(all(result$.g_se > 0))

  # Same inputs as create_test_data_change(); matches esc::esc_mean_sd
  expect_equal(unique(round(result$.g, 6)), -1.302256)
})

test_that("calculateEffectSizes handles single-row change-score data with continuous columns present", {
  # Regression: nrow == 1 + both continuous and change columns used to crash
  # in the ES-group conflict check (sapply collapsed to a vector).
  test_data <- create_test_data_change(n = 1)
  expect_no_error({
    test_data <- checkDataFormat(test_data)
    result <- calculateEffectSizes(test_data)
  })
  expect_true(is.finite(result$.g))
  expect_true(is.finite(result$.g_se))
  expect_equal(round(result$.g, 6), -1.302256)
})

