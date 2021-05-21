# Tests for simulateTreatmentCycles function

test_that("simulateTreatmentCycles validates input", {
  # Test missing response.rate
  expect_error(simulateTreatmentCycles(decay = 0.1),
               "response.rate.*must be provided")
  
  # Test missing decay
  expect_error(simulateTreatmentCycles(response.rate = 0.5),
               "decay.*must be provided")
  
  # Test invalid response.rate (too low)
  expect_error(simulateTreatmentCycles(response.rate = 0.005, decay = 0.1),
               "response.rate.*must be numeric and between 0.01 and 0.99")
  
  # Test invalid response.rate (too high)
  expect_error(simulateTreatmentCycles(response.rate = 0.995, decay = 0.1),
               "response.rate.*must be numeric and between 0.01 and 0.99")
  
  # Test invalid decay (too low)
  expect_error(simulateTreatmentCycles(response.rate = 0.5, decay = -1.0),
               "decay.*must be numeric and between -0.99 and 0.99")
  
  # Test invalid decay (too high)
  expect_error(simulateTreatmentCycles(response.rate = 0.5, decay = 1.0),
               "decay.*must be numeric and between -0.99 and 0.99")
  
  # Test invalid max.cycles
  expect_error(simulateTreatmentCycles(response.rate = 0.5, decay = 0.1, max.cycles = 0),
               "max.cycles.*must be numeric and >= 1")
})

test_that("simulateTreatmentCycles works with simple scenario", {
  # Simple scenario: 50% response rate, no decay
  result <- simulateTreatmentCycles(response.rate = 0.5, decay = 0)
  
  expect_s3_class(result, "simulateTreatmentCycles")
  expect_true("data" %in% names(result))
  expect_true("total.cycles" %in% names(result))
  expect_true("total.treatments" %in% names(result))
  expect_true("excess.treatments" %in% names(result))
  expect_true("avg.no.treatments" %in% names(result))
  expect_s3_class(result$data, "data.frame")
  expect_true(nrow(result$data) > 0)
  expect_true(result$total.cycles > 0)
  expect_true(result$total.treatments > 0)
})

test_that("simulateTreatmentCycles works with decay", {
  # Test with positive decay (decreasing response)
  result <- simulateTreatmentCycles(response.rate = 0.5, decay = 0.1)
  
  expect_s3_class(result, "simulateTreatmentCycles")
  expect_true("data" %in% names(result))
  expect_true(result$total.cycles > 0)
})

test_that("simulateTreatmentCycles works with negative decay", {
  # Test with negative decay (increasing response)
  result <- simulateTreatmentCycles(response.rate = 0.5, decay = -0.1)
  
  expect_s3_class(result, "simulateTreatmentCycles")
  expect_true("data" %in% names(result))
})

test_that("simulateTreatmentCycles works with vector response rates", {
  # Manual response rates
  result <- simulateTreatmentCycles(response.rate = c(0.57, 0.4), decay = 0)
  
  expect_s3_class(result, "simulateTreatmentCycles")
  expect_true("data" %in% names(result))
  expect_true(nrow(result$data) > 0)
})

test_that("simulateTreatmentCycles works with vector response rates and decay", {
  # Manual response rates with decay
  result <- simulateTreatmentCycles(response.rate = c(0.56, 0.3, 0.28, 0.4), decay = 0.1)
  
  expect_s3_class(result, "simulateTreatmentCycles")
  expect_true("data" %in% names(result))
})

test_that("simulateTreatmentCycles works with max.cycles", {
  # Test with capped cycles
  result <- simulateTreatmentCycles(response.rate = 0.5, decay = 0.1, max.cycles = 10)
  
  expect_s3_class(result, "simulateTreatmentCycles")
  expect_true("capped" %in% names(result))
  expect_true(result$max.cycles == 10)
  expect_true(nrow(result$data) <= 10)
})

test_that("print.simulateTreatmentCycles works", {
  result <- simulateTreatmentCycles(response.rate = 0.5, decay = 0)
  
  expect_error(print(result), NA)
  expect_output(print(result), "Total number of treatment cycles")
})

test_that("plot.simulateTreatmentCycles works", {
  result <- simulateTreatmentCycles(response.rate = 0.5, decay = 0)
  
  # Should not error (plot is created)
  expect_error(plot(result), NA)
})

test_that("simulateTreatmentCycles handles edge cases", {
  # Very high response rate
  result <- simulateTreatmentCycles(response.rate = 0.99, decay = 0)
  expect_s3_class(result, "simulateTreatmentCycles")
  
  # Very low response rate
  result <- simulateTreatmentCycles(response.rate = 0.01, decay = 0, max.cycles = 50)
  expect_s3_class(result, "simulateTreatmentCycles")
})

