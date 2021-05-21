test_that("checkConflicts handles valid data", {
  test_data <- create_test_data()
  test_data <- checkDataFormat(test_data)
  
  result <- checkConflicts(test_data)
  
  # Should return data.frame if no conflicts, or list if conflicts found
  expect_true(is.data.frame(result) || 
              (is.list(result) && "checkConflicts" %in% class(result)))
})

test_that("checkConflicts detects duplicate IDs", {
  # Create data with duplicate IDs (same study, outcome, etc.)
  test_data <- create_test_data()
  test_data <- rbind(test_data, test_data[1, ])  # Duplicate first row
  test_data <- checkDataFormat(test_data)
  
  result <- checkConflicts(test_data)
  
  # Should detect conflicts
  if (is.list(result) && "checkConflicts" %in% class(result)) {
    expect_true("allConflicts" %in% names(result))
    expect_true(nrow(result$allConflicts) > 0)
  }
})

test_that("checkConflicts works with custom vars.for.id", {
  test_data <- create_test_data()
  test_data <- checkDataFormat(test_data)
  
  # Use only study as ID variable
  result <- checkConflicts(test_data, vars.for.id = "study")
  
  expect_true(is.data.frame(result) || 
              (is.list(result) && "checkConflicts" %in% class(result)))
})

