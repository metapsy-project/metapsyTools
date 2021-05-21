test_that("checkDataFormat validates required columns", {
  # Create minimal test data
  test_data <- create_test_data()
  
  # Should pass with all required columns
  result <- checkDataFormat(test_data)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), nrow(test_data))
})

test_that("checkDataFormat detects missing columns", {
  # Create data with missing columns
  test_data <- create_test_data()
  test_data$study <- NULL
  
  # Should detect missing column
  expect_message(
    checkDataFormat(test_data),
    "does not contain variable"
  )
})

test_that("checkDataFormat converts variable classes", {
  # Create data with wrong class
  test_data <- create_test_data()
  test_data$time_weeks <- as.character(test_data$time_weeks)
  
  result <- checkDataFormat(test_data)
  expect_type(result$time_weeks, "double")
})

test_that("checkDataFormat renames reserved column names", {
  test_data <- create_test_data()
  test_data$subset <- "test"
  test_data$exclude <- "test"
  
  result <- checkDataFormat(test_data)
  expect_false("subset" %in% names(result))
  expect_false("exclude" %in% names(result))
  expect_true("subset.1" %in% names(result) || "exclude.1" %in% names(result))
})

test_that("checkDataFormat handles custom must.contain", {
  test_data <- data.frame(
    study = c("A", "B"),
    condition = c("T", "C"),
    primary = c(TRUE, FALSE),
    year = c(2020, 2021),
    stringsAsFactors = FALSE
  )
  
  result <- checkDataFormat(
    test_data,
    must.contain = c("study", "condition", "primary", "year"),
    variable.class = list(
      study = "character",
      condition = "character",
      primary = "logical",
      year = "numeric"
    )
  )
  
  expect_s3_class(result, "data.frame")
})

test_that("checkDataFormat handles empty must.contain", {
  test_data <- create_test_data()
  
  result <- checkDataFormat(test_data, must.contain = character(0))
  expect_s3_class(result, "data.frame")
})

test_that("checkDataFormat handles empty variable.class", {
  test_data <- create_test_data()
  
  result <- checkDataFormat(test_data, variable.class = list())
  expect_s3_class(result, "data.frame")
})

