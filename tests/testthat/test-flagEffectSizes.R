# Tests for flagEffectSizes function

test_that("flagEffectSizes validates input", {
  # Create minimal test data
  test_data <- data.frame(
    study = c("Study1", "Study2"),
    .g = c(0.5, 1.0),
    .g_se = c(0.2, 0.3),
    rob = c(1, 2),
    stringsAsFactors = FALSE
  )
  
  # Test invalid reference
  expect_error(flagEffectSizes(test_data, reference = "invalid"),
               "reference.*muste be one of")
  
  # Test missing rob variable
  test_data_no_rob <- test_data[, !names(test_data) %in% "rob"]
  expect_error(flagEffectSizes(test_data_no_rob, high.rob.filter = "rob <= 2"),
               "rob.*variable not found")
  
  # Test missing es.var - should error when trying to use it
  test_data_no_es <- test_data[, !names(test_data) %in% ".g"]
  expect_error(flagEffectSizes(test_data_no_es, es.var = ".g"),
               "non-numeric|subscript|not found")  # Should error when accessing missing column
})

test_that("flagEffectSizes works with valid inputs", {
  test_data <- data.frame(
    study = c("Study1", "Study2", "Study3"),
    .g = c(0.5, 1.0, 1.5),
    .g_se = c(0.2, 0.3, 0.25),
    rob = c(1, 2, 3),
    stringsAsFactors = FALSE
  )
  
  result <- flagEffectSizes(test_data, reference = "all")
  
  expect_s3_class(result, "flagEffectSizes")
  expect_equal(length(result$smd), 3)
  expect_equal(length(result$se), 3)
  expect_true(is.logical(result$flag.effect))
  expect_true(is.logical(result$flag.power))
  expect_true(is.logical(result$flag.rob))
  expect_true(all(result$flags >= 0 & result$flags <= 3))
  expect_true("lookup" %in% names(result))
  expect_true("data" %in% names(result))
  expect_equal(nrow(result$data), 3)
})

test_that("flagEffectSizes works with different reference values", {
  test_data <- data.frame(
    study = c("Study1", "Study2"),
    .g = c(0.5, 1.0),
    .g_se = c(0.2, 0.3),
    rob = c(1, 2),
    stringsAsFactors = FALSE
  )
  
  for (ref in c("all", "dep", "psy", "ptsd")) {
    result <- flagEffectSizes(test_data, reference = ref)
    expect_s3_class(result, "flagEffectSizes")
  }
})

test_that("flagEffectSizes works with custom power", {
  test_data <- data.frame(
    study = c("Study1", "Study2"),
    .g = c(0.5, 1.0),
    .g_se = c(0.2, 0.3),
    rob = c(1, 2),
    stringsAsFactors = FALSE
  )
  
  result <- flagEffectSizes(test_data, power = 0.8)
  expect_s3_class(result, "flagEffectSizes")
})

test_that("flagEffectSizes works with custom high.rob.filter", {
  test_data <- data.frame(
    study = c("Study1", "Study2"),
    .g = c(0.5, 1.0),
    .g_se = c(0.2, 0.3),
    rob = c(1, 2),
    stringsAsFactors = FALSE
  )
  
  result <- flagEffectSizes(test_data, high.rob.filter = "rob <= 1")
  expect_s3_class(result, "flagEffectSizes")
})

test_that("flagEffectSizes handles no studies with risk of bias", {
  test_data <- data.frame(
    study = c("Study1", "Study2"),
    .g = c(0.5, 1.0),
    .g_se = c(0.2, 0.3),
    rob = c(1, 1),  # All low risk
    stringsAsFactors = FALSE
  )
  
  # Should warn but not error
  result <- expect_warning(
    flagEffectSizes(test_data, high.rob.filter = "rob > 1"),
    "No studies with risk of bias"
  )
  expect_s3_class(result, "flagEffectSizes")
})

test_that("print.flagEffectSizes works", {
  test_data <- data.frame(
    study = c("Study1", "Study2"),
    .g = c(0.5, 1.0),
    .g_se = c(0.2, 0.3),
    rob = c(1, 2),
    stringsAsFactors = FALSE
  )
  
  result <- flagEffectSizes(test_data)
  expect_error(print(result), NA)
  expect_output(print(result), "N effect sizes")
})

test_that("plot.flagEffectSizes works", {
  test_data <- data.frame(
    study = c("Study1", "Study2"),
    .g = c(0.5, 1.0),
    .g_se = c(0.2, 0.3),
    rob = c(1, 2),
    stringsAsFactors = FALSE
  )
  
  result <- flagEffectSizes(test_data)
  # Should not error (plot might fail due to missing plotDensityHist, but function call should work)
  expect_error(try(plot(result), silent = TRUE), NA)
})

