# Tests for addTrialArmInfo function (deprecated but should still work)

test_that("addTrialArmInfo validates input", {
  # Test invalid input type
  expect_error(addTrialArmInfo("not a data.frame"),
               "must be of class 'data.frame'")
})

test_that("addTrialArmInfo works with valid expanded multiarm data", {
  # Create test data that mimics expanded multiarm trial structure
  test_data <- data.frame(
    study = rep(c("Study1", "Study2"), each = 2),
    condition = rep(c("ig", "cg"), 2),
    primary = rep(1, 4),
    Outc_measure = rep("phq-9", 4),
    Time = rep("post", 4),
    Time_weeks = rep(8, 4),
    Post_N = c(50, 48, 45, 47),  # Trial arm variable
    stringsAsFactors = FALSE
  )
  
  # Should warn about deprecation
  expect_warning(
    result <- addTrialArmInfo(test_data, Post_N),
    "deprecated"
  )
  
  expect_s3_class(result, "data.frame")
  # Should have added columns for Post_N
  expect_true(any(grepl("Post_N", names(result))))
})

test_that("addTrialArmInfo works with multiple variables", {
  test_data <- data.frame(
    study = rep(c("Study1", "Study2"), each = 2),
    condition = rep(c("ig", "cg"), 2),
    primary = rep(1, 4),
    Outc_measure = rep("phq-9", 4),
    Time = rep("post", 4),
    Time_weeks = rep(8, 4),
    Post_N = c(50, 48, 45, 47),
    Rand_N = c(52, 50, 47, 49),
    stringsAsFactors = FALSE
  )
  
  expect_warning(
    result <- addTrialArmInfo(test_data, Post_N, Rand_N),
    "deprecated"
  )
  
  expect_s3_class(result, "data.frame")
})

test_that("addTrialArmInfo works with custom group names", {
  test_data <- data.frame(
    study = rep(c("Study1", "Study2"), each = 2),
    condition = rep(c("treatment", "control"), 2),
    primary = rep(1, 4),
    Outc_measure = rep("phq-9", 4),
    Time = rep("post", 4),
    Time_weeks = rep(8, 4),
    Post_N = c(50, 48, 45, 47),
    stringsAsFactors = FALSE
  )
  
  expect_warning(
    result <- addTrialArmInfo(test_data, Post_N,
                              .group.indicator = "condition",
                              .name.intervention.group = "treatment",
                              .name.control.group = "control"),
    "deprecated"
  )
  
  expect_s3_class(result, "data.frame")
})

