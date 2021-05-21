# Tests for expandMultiarmTrials function (deprecated but should still work)

test_that("expandMultiarmTrials validates input", {
  # Test with checkConflicts object
  conflicts_obj <- structure(list(), class = "checkConflicts")
  expect_error(expandMultiarmTrials(conflicts_obj),
               "Data format conflicts were detected")
})

test_that("expandMultiarmTrials works with long format data", {
  # Create simple test data in long format
  test_data <- data.frame(
    study = c("Study1", "Study1", "Study2", "Study2"),
    condition = c("ig", "cg", "ig", "cg"),
    primary = rep(1, 4),
    Outc_measure = rep("phq-9", 4),
    Time = rep("post", 4),
    Time_weeks = rep(8, 4),
    sr_clinician = rep(0, 4),  # Add required column
    is.multiarm = rep(0, 4),
    no.arms = rep(2, 4),
    Cond_spec = c("cbt", "wl", "cbt", "wl"),
    stringsAsFactors = FALSE
  )
  
  expect_warning(
    result <- expandMultiarmTrials(test_data, data.format = "long"),
    "deprecated"
  )
  
  expect_s3_class(result, "data.frame")
  # Should have id and study.id columns
  expect_true("id" %in% names(result) || "study.id" %in% names(result))
})

test_that("expandMultiarmTrials works with custom parameters", {
  test_data <- data.frame(
    study = c("Study1", "Study1", "Study2", "Study2"),
    condition = c("treatment", "control", "treatment", "control"),
    primary = rep(1, 4),
    Outc_measure = rep("phq-9", 4),
    Time = rep("post", 4),
    Time_weeks = rep(8, 4),
    sr_clinician = rep(0, 4),  # Add required column
    is.multiarm = rep(0, 4),
    no.arms = rep(2, 4),
    Cond_spec = c("cbt", "wl", "cbt", "wl"),
    stringsAsFactors = FALSE
  )
  
  expect_warning(
    result <- expandMultiarmTrials(test_data,
                                   data.format = "long",
                                   group.indicator = "condition",
                                   group.names = list("ig" = "treatment",
                                                      "cg" = "control")),
    "deprecated"
  )
  
  expect_s3_class(result, "data.frame")
})

test_that("expandMultiarmTrials works with wide format", {
  # Create wide format data
  test_data <- data.frame(
    study = c("Study1", "Study2"),
    condition_arm1 = c("cbt", "cbt"),
    condition_arm2 = c("wl", "wl"),
    primary = rep(1, 2),
    Outc_measure = rep("phq-9", 2),
    Time = rep("post", 2),
    Time_weeks = rep(8, 2),
    sr_clinician = rep(0, 2),
    stringsAsFactors = FALSE
  )
  
  expect_warning(
    result <- expandMultiarmTrials(test_data, data.format = "wide"),
    "deprecated"
  )
  
  expect_s3_class(result, "data.frame")
  expect_true("id" %in% names(result))
  expect_true("study.id" %in% names(result))
})

test_that("expandMultiarmTrials handles missing variables", {
  test_data <- data.frame(
    study = c("Study1", "Study2"),
    condition = c("ig", "cg"),
    stringsAsFactors = FALSE
  )
  
  # Should error about missing variables
  expect_error(expandMultiarmTrials(test_data, data.format = "long"),
               "variable.*not found")
})

