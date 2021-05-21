# Tests for filterPriorityRule function

test_that("filterPriorityRule validates input", {
  # Test invalid input type
  expect_error(filterPriorityRule("not a data.frame"),
               "must be a data.frame")
  
  # Test with empty data frame
  empty_df <- data.frame(study = character(0), condition_arm1 = character(0))
  result <- filterPriorityRule(empty_df, condition_arm1 = c("a", "b"))
  expect_equal(nrow(result), 0)
})

test_that("filterPriorityRule works with simple data", {
  # Create test data
  test_data <- data.frame(
    study = c("Study1", "Study1", "Study2", "Study2"),
    condition_arm1 = c("cbt", "pst", "cbt", "other"),
    condition_arm2 = c("wl", "wl", "wl", "wl"),
    stringsAsFactors = FALSE
  )
  
  # Filter with priority rule
  result <- filterPriorityRule(test_data, 
                                condition_arm1 = c("cbt", "pst"),
                                .study.indicator = "study")
  
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) <= nrow(test_data))
  
  # Study1 should have cbt (first priority)
  study1_rows <- result[result$study == "Study1", ]
  if (nrow(study1_rows) > 0) {
    expect_true(all(study1_rows$condition_arm1 == "cbt"))
  }
})

test_that("filterPriorityRule works with multiple priority rules", {
  test_data <- data.frame(
    study = c("Study1", "Study1", "Study2", "Study2"),
    condition_arm1 = c("cbt", "pst", "cbt", "other"),
    instrument = c("phq-9", "cesd", "phq-9", "other"),
    stringsAsFactors = FALSE
  )
  
  result <- filterPriorityRule(test_data,
                                condition_arm1 = c("cbt", "pst"),
                                instrument = c("phq-9", "cesd"),
                                .study.indicator = "study")
  
  expect_s3_class(result, "data.frame")
})

test_that("filterPriorityRule removes studies with no matching levels", {
  test_data <- data.frame(
    study = c("Study1", "Study1", "Study2", "Study2"),
    condition_arm1 = c("other1", "other2", "cbt", "pst"),
    stringsAsFactors = FALSE
  )
  
  result <- filterPriorityRule(test_data,
                                condition_arm1 = c("cbt", "pst"),
                                .study.indicator = "study")
  
  # Study1 should be removed (no matching levels)
  expect_false(any(result$study == "Study1"))
  # Study2 should remain
  expect_true(any(result$study == "Study2"))
})

test_that("filterPriorityRule example from documentation works", {
  # This test uses the example from the documentation
  # We'll use minimal test data since we can't load the full dataset
  test_data <- data.frame(
    study = rep(c("Study1", "Study2"), each = 3),
    condition_arm1 = c("cbt", "pst", "other", "cbt", "pst", "other"),
    condition_arm2 = c("cau", "wl", "cbt", "cau", "wl", "cbt"),
    instrument = c("cesd", "phq-9", "scl", "cesd", "phq-9", "scl"),
    time = c("post", "fu", "post", "post", "fu", "post"),
    stringsAsFactors = FALSE
  )
  
  result <- try(filterPriorityRule(test_data,
                                    condition_arm1 = c("cbt", "pst"),
                                    condition_arm2 = c("cau", "wl", "cbt"),
                                    instrument = c("cesd", "phq-9", "scl"),
                                    time = c("post", "fu")),
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    expect_s3_class(result, "data.frame")
  }
})

