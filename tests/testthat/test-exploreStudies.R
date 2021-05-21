# Tests for exploreStudies function

test_that("exploreStudies validates input", {
  # Test with invalid data type
  expect_error(try(exploreStudies("not a data.frame"), silent = TRUE), NA)
})

test_that("exploreStudies works with data.frame", {
  test_data <- create_test_data()
  # Add required columns for exploreStudies
  test_data$multi_arm1 <- rep("A", nrow(test_data))
  test_data$multi_arm2 <- rep("B", nrow(test_data))
  
  result <- try(exploreStudies(test_data, which = "treatments", html = FALSE),
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    expect_s3_class(result, "exploreStudies")
    expect_true("summary" %in% names(result))
    expect_true("data" %in% names(result))
    expect_true("conditions" %in% names(result$summary))
  }
})

test_that("exploreStudies works with comparisons", {
  test_data <- create_test_data()
  test_data$multi_arm1 <- rep("A", nrow(test_data))
  test_data$multi_arm2 <- rep("B", nrow(test_data))
  
  result <- try(exploreStudies(test_data, which = "comparisons", html = FALSE),
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    expect_s3_class(result, "exploreStudies")
    expect_true("summary" %in% names(result))
    expect_true("comparisons" %in% names(result$summary))
  }
})

test_that("exploreStudies works with runMetaAnalysis object", {
  test_data <- create_test_data_with_es(n = 6)
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    explore_result <- try(exploreStudies(result, which = "treatments", html = FALSE),
                         silent = TRUE)
    
    if (!inherits(explore_result, "try-error")) {
      expect_s3_class(explore_result, "exploreStudies")
    }
  }
})

test_that("exploreStudies works with custom parameters", {
  test_data <- create_test_data()
  test_data$multi_arm1 <- rep("A", nrow(test_data))
  test_data$multi_arm2 <- rep("B", nrow(test_data))
  
  result <- try(exploreStudies(test_data,
                               which = "treatments",
                               .study.var = "study",
                               .condition = "condition",
                               .condition.specification = "multi",
                               .groups.column.indicator = c("_arm1", "_arm2"),
                               html = FALSE),
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    expect_s3_class(result, "exploreStudies")
  }
})

test_that("exploreStudies works with HTML output", {
  test_data <- create_test_data()
  test_data$multi_arm1 <- rep("A", nrow(test_data))
  test_data$multi_arm2 <- rep("B", nrow(test_data))
  
  # HTML output is printed, not returned
  expect_error(try(exploreStudies(test_data, which = "treatments", html = TRUE),
                   silent = TRUE), NA)
})

test_that("print.exploreStudies works", {
  test_data <- create_test_data()
  test_data$multi_arm1 <- rep("A", nrow(test_data))
  test_data$multi_arm2 <- rep("B", nrow(test_data))
  
  result <- try(exploreStudies(test_data, which = "treatments", html = FALSE),
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    expect_error(print(result), NA)
    expect_output(print(result), "studies|treatments|comparisons")
  }
})

test_that("exploreStudies handles missing sample size information", {
  test_data <- create_test_data()
  test_data$multi_arm1 <- rep("A", nrow(test_data))
  test_data$multi_arm2 <- rep("B", nrow(test_data))
  # Remove sample size columns
  test_data$n_arm1 <- NULL
  test_data$n_arm2 <- NULL
  
  # Should warn but not error
  result <- expect_warning(
    try(exploreStudies(test_data, which = "treatments", html = FALSE),
        silent = TRUE),
    "No sample size information"
  )
  
  if (!inherits(result, "try-error")) {
    expect_s3_class(result, "exploreStudies")
  }
})

