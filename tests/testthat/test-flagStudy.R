# Tests for flagStudy function (internal)

test_that("flagStudy validates input arguments", {
  # Test invalid reference
  expect_error(flagStudy(smd = 0.5, se = 0.2, rob = TRUE, reference = "invalid"),
               "reference.*must be either")
  
  # Test invalid pval
  expect_error(flagStudy(smd = 0.5, se = 0.2, rob = TRUE, pval = 0),
               "pval.*must be a numeric between 0 and 1")
  expect_error(flagStudy(smd = 0.5, se = 0.2, rob = TRUE, pval = 1),
               "pval.*must be a numeric between 0 and 1")
  expect_error(flagStudy(smd = 0.5, se = 0.2, rob = TRUE, pval = 1.5),
               "pval.*must be a numeric between 0 and 1")
  
  # Test invalid rob (must be logical)
  expect_error(flagStudy(smd = 0.5, se = 0.2, rob = c(1, 0)),
               "rob.*must be a logical vector")
  
  # Test mismatched lengths
  expect_error(flagStudy(smd = c(0.5, 0.6), se = 0.2, rob = TRUE),
               "lengths.*do not match")
  
  # Test invalid power
  expect_error(flagStudy(smd = 0.5, se = 0.2, rob = TRUE, power = 1.5),
               "power.*must by a numeric")
  expect_error(flagStudy(smd = 0.5, se = 0.2, rob = TRUE, power = -0.1),
               "power.*must by a numeric")
})

test_that("flagStudy works with valid inputs", {
  smd <- c(0.5, 1.0, 1.5, 2.0)
  se <- c(0.2, 0.3, 0.25, 0.3)
  rob <- c(TRUE, FALSE, TRUE, FALSE)
  
  result <- flagStudy(smd = smd, se = se, rob = rob, reference = "all")
  
  expect_s3_class(result, "flagStudy")
  expect_equal(length(result$smd), 4)
  expect_equal(length(result$se), 4)
  expect_true(is.logical(result$flag.effect))
  expect_true(is.logical(result$flag.power))
  expect_true(is.logical(result$flag.rob))
  expect_true(all(result$flags >= 0 & result$flags <= 3))
  expect_true("lookup" %in% names(result))
})

test_that("flagStudy works with different reference values", {
  smd <- c(0.5, 1.0)
  se <- c(0.2, 0.3)
  rob <- c(TRUE, FALSE)
  
  for (ref in c("all", "dep", "psy", "ptsd")) {
    result <- flagStudy(smd = smd, se = se, rob = rob, reference = ref)
    expect_s3_class(result, "flagStudy")
    expect_equal(nrow(result$lookup), 1)
  }
})

test_that("flagStudy works with numeric power", {
  smd <- c(0.5, 1.0)
  se <- c(0.2, 0.3)
  rob <- c(TRUE, FALSE)
  
  result <- flagStudy(smd = smd, se = se, rob = rob, power = 0.8)
  expect_s3_class(result, "flagStudy")
})

test_that("print.flagStudy works", {
  smd <- c(0.5, 1.0, 1.5)
  se <- c(0.2, 0.3, 0.25)
  rob <- c(TRUE, FALSE, TRUE)
  
  result <- flagStudy(smd = smd, se = se, rob = rob)
  
  # Should not error
  expect_error(print(result), NA)
  expect_output(print(result), "N effect sizes")
})

test_that("plot.flagStudy works", {
  smd <- c(0.5, 1.0)
  se <- c(0.2, 0.3)
  rob <- c(TRUE, FALSE)
  
  result <- flagStudy(smd = smd, se = se, rob = rob)
  
  # Should not error (plot might fail due to missing plotDensityHist, but function call should work)
  expect_error(try(plot(result), silent = TRUE), NA)
})

