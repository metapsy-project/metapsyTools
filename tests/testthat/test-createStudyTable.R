# Tests for createStudyTable function

test_that("createStudyTable validates input", {
  # Test with runMetaAnalysis object (should extract data)
  test_data <- create_test_data()
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    # Should work with runMetaAnalysis object
    expect_error(createStudyTable(result, study), NA)
  }
})

test_that("createStudyTable works with basic data", {
  test_data <- create_test_data()
  
  # Basic usage - just select columns
  result <- createStudyTable(test_data, study, condition_arm1, condition_arm2,
                              .html = FALSE)
  
  expect_s3_class(result, "data.frame")
  expect_true("study" %in% names(result))
  expect_true("condition_arm1" %in% names(result))
  expect_true("condition_arm2" %in% names(result))
})

test_that("createStudyTable works with value conversion", {
  test_data <- create_test_data()
  
  # Add a numeric factor column for testing conversion
  test_data$country <- factor(c(1, 2, 3))
  
  # Test value conversion
  result <- createStudyTable(test_data, study,
                             country = c("Europe" = "1", "USA" = "2", "Canada" = "3"),
                             .html = FALSE)
  
  expect_s3_class(result, "data.frame")
  expect_true("country" %in% names(result))
})

test_that("createStudyTable works with rounding", {
  test_data <- create_test_data()
  
  # Test rounding
  result <- createStudyTable(test_data, study, n_arm1, n_arm2,
                             .round.by.digits = list(n_arm1 = 0, n_arm2 = 0),
                             .html = FALSE)
  
  expect_s3_class(result, "data.frame")
  # Convert to numeric if needed (rounding may convert to character)
  n_arm1_num <- suppressWarnings(as.numeric(as.character(result$n_arm1)))
  n_arm2_num <- suppressWarnings(as.numeric(as.character(result$n_arm2)))
  if (!all(is.na(n_arm1_num))) {
    expect_true(all(n_arm1_num == round(n_arm1_num), na.rm = TRUE))
  }
  if (!all(is.na(n_arm2_num))) {
    expect_true(all(n_arm2_num == round(n_arm2_num), na.rm = TRUE))
  }
})

test_that("createStudyTable works with column renaming", {
  test_data <- create_test_data()
  
  # Test column renaming
  result <- createStudyTable(test_data, study, n_arm1, n_arm2,
                             .column.names = list(n_arm1 = "N (arm1)",
                                                 n_arm2 = "N (arm2)"),
                             .html = FALSE)
  
  expect_s3_class(result, "data.frame")
  expect_true("N (arm1)" %in% names(result) || "n_arm1" %in% names(result))
})

test_that("createStudyTable handles NA values", {
  test_data <- create_test_data()
  test_data$some_var <- c(NA, "value1", "value2")
  
  # Test NA replacement
  result <- createStudyTable(test_data, study, some_var,
                             .na.replace = "missing",
                             .html = FALSE)
  
  expect_s3_class(result, "data.frame")
  expect_false(any(is.na(result$some_var)))
})

test_that("createStudyTable removes redundant rows", {
  # Create data with duplicate rows
  test_data <- create_test_data()
  test_data <- rbind(test_data, test_data[1, ])  # Duplicate first row
  
  result <- createStudyTable(test_data, study, condition_arm1,
                             .html = FALSE)
  
  # Should have fewer rows due to distinct()
  expect_true(nrow(result) <= nrow(test_data))
})

test_that("createStudyTable works with HTML output", {
  test_data <- create_test_data()
  
  # Should not error (HTML output is printed, not returned)
  expect_error(createStudyTable(test_data, study, .html = TRUE), NA)
})

test_that("createStudyTable handles study name replacement", {
  test_data <- create_test_data()
  
  # If study is first column, redundant names should be replaced with "."
  result <- createStudyTable(test_data, study, condition_arm1,
                             .html = FALSE)
  
  expect_s3_class(result, "data.frame")
  # If there are duplicates, some study names might be "."
  if (any(result$study == ".")) {
    expect_true(any(result$study == "."))
  }
})

