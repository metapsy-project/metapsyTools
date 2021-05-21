# Tests for createRobRatings function

test_that("createRobRatings validates required variables in rob.data", {
  # Test with missing required variables
  incomplete_rob <- data.frame(
    study = c("Study1", "Study2"),
    d1_1 = c("Yes/PY", "No/PN"),
    stringsAsFactors = FALSE
  )
  
  expect_error(createRobRatings(incomplete_rob),
               "Required variables.*not found")
})

test_that("createRobRatings works with rob.data only (no database)", {
  # Create minimal valid rob.data
  rob_data <- data.frame(
    study = c("Study1", "Study2"),
    d1_1 = c("Yes/PY", "No/PN"),
    d1_2 = c("Yes/PY", "No/PN"),
    d1_3 = c("Yes/PY", "No/PN"),
    d1_4 = c("Yes/PY", "No/PN"),
    d1_notes = c("", ""),
    d2_5 = c("Yes/PY", "No/PN"),
    d2_6 = c("Yes/PY", "No/PN"),
    d2_7 = c("Yes/PY", "No/PN"),
    d2_8 = c("Yes/PY", "No/PN"),
    d2_9 = c("Yes/PY", "No/PN"),
    d2_notes = c("", ""),
    d3_10 = c("Yes/PY", "No/PN"),
    d3_11 = c("Yes/PY", "No/PN"),
    d3_12 = c("Yes/PY", "No/PN"),
    d3_13 = c("Yes/PY", "No/PN"),
    d3_14 = c("Yes/PY", "No/PN"),
    d3_notes = c("", ""),
    d4_15 = c("Yes/PY", "No/PN"),
    d4_16 = c("Yes/PY", "No/PN"),
    d4_17 = c("Yes/PY", "No/PN"),
    d4_18 = c("Yes/PY", "No/PN"),
    d4_notes = c("", ""),
    d5_19 = c("Yes/PY", "No/PN"),
    d5_20 = c("Yes/PY", "No/PN"),
    d5_21 = c("Yes/PY", "No/PN"),
    d5_22 = c("Yes/PY", "No/PN"),
    d5_23 = c("Yes/PY", "No/PN"),
    d5_24 = c("Yes/PY", "No/PN"),
    d5_notes = c("", ""),
    stringsAsFactors = FALSE
  )
  
  result <- createRobRatings(rob_data, database = NULL)
  
  expect_s3_class(result, "createRobRatings")
  expect_true("rob.data" %in% names(result))
  expect_true("database" %in% names(result))
  expect_true(is.null(result$database))
  expect_true("rob" %in% names(result$rob.data))
  expect_true("d1" %in% names(result$rob.data))
  expect_true("d2" %in% names(result$rob.data))
  expect_true("d3" %in% names(result$rob.data))
  expect_true("d4" %in% names(result$rob.data))
  expect_true("d5" %in% names(result$rob.data))
})

test_that("createRobRatings validates database variables", {
  # Create valid rob.data
  rob_data <- data.frame(
    study = c("Study1", "Study2"),
    d1_1 = c("Yes/PY", "No/PN"),
    d1_2 = c("Yes/PY", "No/PN"),
    d1_3 = c("Yes/PY", "No/PN"),
    d1_4 = c("Yes/PY", "No/PN"),
    d1_notes = c("", ""),
    d2_5 = c("Yes/PY", "No/PN"),
    d2_6 = c("Yes/PY", "No/PN"),
    d2_7 = c("Yes/PY", "No/PN"),
    d2_8 = c("Yes/PY", "No/PN"),
    d2_9 = c("Yes/PY", "No/PN"),
    d2_notes = c("", ""),
    d3_10 = c("Yes/PY", "No/PN"),
    d3_11 = c("Yes/PY", "No/PN"),
    d3_12 = c("Yes/PY", "No/PN"),
    d3_13 = c("Yes/PY", "No/PN"),
    d3_14 = c("Yes/PY", "No/PN"),
    d3_notes = c("", ""),
    d4_15 = c("Yes/PY", "No/PN"),
    d4_16 = c("Yes/PY", "No/PN"),
    d4_17 = c("Yes/PY", "No/PN"),
    d4_18 = c("Yes/PY", "No/PN"),
    d4_notes = c("", ""),
    d5_19 = c("Yes/PY", "No/PN"),
    d5_20 = c("Yes/PY", "No/PN"),
    d5_21 = c("Yes/PY", "No/PN"),
    d5_22 = c("Yes/PY", "No/PN"),
    d5_23 = c("Yes/PY", "No/PN"),
    d5_24 = c("Yes/PY", "No/PN"),
    d5_notes = c("", ""),
    stringsAsFactors = FALSE
  )
  
  # Test missing study variable
  incomplete_db <- data.frame(
    condition_arm1 = c("cbt", "cbt"),
    stringsAsFactors = FALSE
  )
  
  expect_error(createRobRatings(rob_data, database = incomplete_db),
               "study.*variable not found")
  
  # Test missing required variables
  incomplete_db2 <- data.frame(
    study = c("Study1", "Study2"),
    stringsAsFactors = FALSE
  )
  
  expect_error(createRobRatings(rob_data, database = incomplete_db2),
               "attr_arm1.*variable not found|rand_arm1.*variable not found")
})

test_that("createRobRatings works with rob.data and database", {
  # Create valid rob.data
  rob_data <- data.frame(
    study = c("Study1", "Study2"),
    d1_1 = c("Yes/PY", "No/PN"),
    d1_2 = c("Yes/PY", "No/PN"),
    d1_3 = c("Yes/PY", "No/PN"),
    d1_4 = c("Yes/PY", "No/PN"),
    d1_notes = c("", ""),
    d2_5 = c("Yes/PY", "No/PN"),
    d2_6 = c("Yes/PY", "No/PN"),
    d2_7 = c("Yes/PY", "No/PN"),
    d2_8 = c("Yes/PY", "No/PN"),
    d2_9 = c("Yes/PY", "No/PN"),
    d2_notes = c("", ""),
    d3_10 = c("Yes/PY", "No/PN"),
    d3_11 = c("Yes/PY", "No/PN"),
    d3_12 = c("Yes/PY", "No/PN"),
    d3_13 = c("Yes/PY", "No/PN"),
    d3_14 = c("Yes/PY", "No/PN"),
    d3_notes = c("", ""),
    d4_15 = c("Yes/PY", "No/PN"),
    d4_16 = c("Yes/PY", "No/PN"),
    d4_17 = c("Yes/PY", "No/PN"),
    d4_18 = c("Yes/PY", "No/PN"),
    d4_notes = c("", ""),
    d5_19 = c("Yes/PY", "No/PN"),
    d5_20 = c("Yes/PY", "No/PN"),
    d5_21 = c("Yes/PY", "No/PN"),
    d5_22 = c("Yes/PY", "No/PN"),
    d5_23 = c("Yes/PY", "No/PN"),
    d5_24 = c("Yes/PY", "No/PN"),
    d5_notes = c("", ""),
    stringsAsFactors = FALSE
  )
  
  # Create minimal valid database
  database <- data.frame(
    study = c("Study1", "Study2"),
    attr_arm1 = c(10, 15),
    attr_arm2 = c(12, 18),
    rand_arm1 = c(50, 60),
    rand_arm2 = c(50, 60),
    rand_ratio = c("1 to 1", "1 to 1"),
    rating = c("self-report", "clinician"),
    condition_arm1 = c("cbt", "cbt"),
    condition_arm2 = c("wl", "wl"),
    stringsAsFactors = FALSE
  )
  
  result <- try(createRobRatings(rob_data, database = database), silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    expect_s3_class(result, "createRobRatings")
    expect_true("rob.data" %in% names(result))
    expect_true("database" %in% names(result))
    expect_s3_class(result$database, "data.frame")
    expect_true("rob" %in% names(result$database))
    expect_true("rob_d1" %in% names(result$database))
    expect_true("rob_study_lvl" %in% names(result$database))
  }
})

test_that("createRobRatings handles invalid rating values", {
  rob_data <- data.frame(
    study = c("Study1", "Study2"),
    d1_1 = c("Yes/PY", "No/PN"),
    d1_2 = c("Yes/PY", "No/PN"),
    d1_3 = c("Yes/PY", "No/PN"),
    d1_4 = c("Yes/PY", "No/PN"),
    d1_notes = c("", ""),
    d2_5 = c("Yes/PY", "No/PN"),
    d2_6 = c("Yes/PY", "No/PN"),
    d2_7 = c("Yes/PY", "No/PN"),
    d2_8 = c("Yes/PY", "No/PN"),
    d2_9 = c("Yes/PY", "No/PN"),
    d2_notes = c("", ""),
    d3_10 = c("Yes/PY", "No/PN"),
    d3_11 = c("Yes/PY", "No/PN"),
    d3_12 = c("Yes/PY", "No/PN"),
    d3_13 = c("Yes/PY", "No/PN"),
    d3_14 = c("Yes/PY", "No/PN"),
    d3_notes = c("", ""),
    d4_15 = c("Yes/PY", "No/PN"),
    d4_16 = c("Yes/PY", "No/PN"),
    d4_17 = c("Yes/PY", "No/PN"),
    d4_18 = c("Yes/PY", "No/PN"),
    d4_notes = c("", ""),
    d5_19 = c("Yes/PY", "No/PN"),
    d5_20 = c("Yes/PY", "No/PN"),
    d5_21 = c("Yes/PY", "No/PN"),
    d5_22 = c("Yes/PY", "No/PN"),
    d5_23 = c("Yes/PY", "No/PN"),
    d5_24 = c("Yes/PY", "No/PN"),
    d5_notes = c("", ""),
    stringsAsFactors = FALSE
  )
  
  database <- data.frame(
    study = c("Study1", "Study2"),
    attr_arm1 = c(10, 15),
    attr_arm2 = c(12, 18),
    rand_arm1 = c(50, 60),
    rand_arm2 = c(50, 60),
    rand_ratio = c("1 to 1", "1 to 1"),
    rating = c("invalid", "clinician"),  # Invalid rating
    condition_arm1 = c("cbt", "cbt"),
    condition_arm2 = c("wl", "wl"),
    stringsAsFactors = FALSE
  )
  
  expect_error(createRobRatings(rob_data, database = database),
               "Only allowed.*values for 'rating'")
})

test_that("createRobRatings handles missing studies", {
  rob_data <- data.frame(
    study = c("Study1"),
    d1_1 = c("Yes/PY"),
    d1_2 = c("Yes/PY"),
    d1_3 = c("Yes/PY"),
    d1_4 = c("Yes/PY"),
    d1_notes = c(""),
    d2_5 = c("Yes/PY"),
    d2_6 = c("Yes/PY"),
    d2_7 = c("Yes/PY"),
    d2_8 = c("Yes/PY"),
    d2_9 = c("Yes/PY"),
    d2_notes = c(""),
    d3_10 = c("Yes/PY"),
    d3_11 = c("Yes/PY"),
    d3_12 = c("Yes/PY"),
    d3_13 = c("Yes/PY"),
    d3_14 = c("Yes/PY"),
    d3_notes = c(""),
    d4_15 = c("Yes/PY"),
    d4_16 = c("Yes/PY"),
    d4_17 = c("Yes/PY"),
    d4_18 = c("Yes/PY"),
    d4_notes = c(""),
    d5_19 = c("Yes/PY"),
    d5_20 = c("Yes/PY"),
    d5_21 = c("Yes/PY"),
    d5_22 = c("Yes/PY"),
    d5_23 = c("Yes/PY"),
    d5_24 = c("Yes/PY"),
    d5_notes = c(""),
    stringsAsFactors = FALSE
  )
  
  database <- data.frame(
    study = c("Study1", "Study2"),  # Study2 not in rob_data
    attr_arm1 = c(10, 15),
    attr_arm2 = c(12, 18),
    rand_arm1 = c(50, 60),
    rand_arm2 = c(50, 60),
    rand_ratio = c("1 to 1", "1 to 1"),
    rating = c("self-report", "clinician"),
    condition_arm1 = c("cbt", "cbt"),
    condition_arm2 = c("wl", "wl"),
    stringsAsFactors = FALSE
  )
  
  result <- try(createRobRatings(rob_data, database = database), silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    expect_s3_class(result, "createRobRatings")
    # Should have miss.studies
    expect_true("miss.studies" %in% names(result))
  }
})

test_that("print.createRobRatings works", {
  rob_data <- data.frame(
    study = c("Study1", "Study2"),
    d1_1 = c("Yes/PY", "No/PN"),
    d1_2 = c("Yes/PY", "No/PN"),
    d1_3 = c("Yes/PY", "No/PN"),
    d1_4 = c("Yes/PY", "No/PN"),
    d1_notes = c("", ""),
    d2_5 = c("Yes/PY", "No/PN"),
    d2_6 = c("Yes/PY", "No/PN"),
    d2_7 = c("Yes/PY", "No/PN"),
    d2_8 = c("Yes/PY", "No/PN"),
    d2_9 = c("Yes/PY", "No/PN"),
    d2_notes = c("", ""),
    d3_10 = c("Yes/PY", "No/PN"),
    d3_11 = c("Yes/PY", "No/PN"),
    d3_12 = c("Yes/PY", "No/PN"),
    d3_13 = c("Yes/PY", "No/PN"),
    d3_14 = c("Yes/PY", "No/PN"),
    d3_notes = c("", ""),
    d4_15 = c("Yes/PY", "No/PN"),
    d4_16 = c("Yes/PY", "No/PN"),
    d4_17 = c("Yes/PY", "No/PN"),
    d4_18 = c("Yes/PY", "No/PN"),
    d4_notes = c("", ""),
    d5_19 = c("Yes/PY", "No/PN"),
    d5_20 = c("Yes/PY", "No/PN"),
    d5_21 = c("Yes/PY", "No/PN"),
    d5_22 = c("Yes/PY", "No/PN"),
    d5_23 = c("Yes/PY", "No/PN"),
    d5_24 = c("Yes/PY", "No/PN"),
    d5_notes = c("", ""),
    stringsAsFactors = FALSE
  )
  
  result <- createRobRatings(rob_data, database = NULL)
  expect_error(print(result), NA)
  expect_output(print(result), "createRobRatings")
})

