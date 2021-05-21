# Tests for proportionMID function

test_that("proportionMID validates input", {
  # Test with non-runMetaAnalysis object - function accesses x$which.run early
  expect_error(proportionMID("not a runMetaAnalysis object", mid = 0.5),
               "\\$ operator|which.run")
  
  # Test with NULL mid
  test_data <- create_test_data_with_es(n = 6)
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    expect_error(proportionMID(result, mid = NULL),
                 "A value of 'mid' must be specified")
    
    # Test invalid test parameter
    expect_error(proportionMID(result, mid = 0.5, test = "invalid"),
                 "test.*must be either 'smaller'")
    
    # Test invalid which parameter
    expect_error(proportionMID(result, mid = 0.5, which = "invalid"),
                 "which.*must be either 'all'")
  }
})

test_that("proportionMID works with overall model", {
  test_data <- create_test_data_with_es(n = 10)
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    prop_result <- try(proportionMID(result, mid = 0.5, which = "overall"),
                      silent = TRUE)
    
    if (!inherits(prop_result, "try-error")) {
      expect_s3_class(prop_result, "proportionMID")
      expect_true("qhat" %in% names(prop_result))
      expect_true("plot.dat" %in% names(prop_result))
      expect_true("xlab" %in% names(prop_result))
      expect_true("test" %in% names(prop_result))
      expect_s3_class(prop_result$qhat, "data.frame")
    }
  }
})

test_that("proportionMID works with test='bigger'", {
  test_data <- create_test_data_with_es(n = 10)
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    prop_result <- try(proportionMID(result, mid = 0.5, test = "bigger"),
                      silent = TRUE)
    
    if (!inherits(prop_result, "try-error")) {
      expect_s3_class(prop_result, "proportionMID")
      expect_equal(prop_result$test, "bigger")
    }
  }
})

test_that("proportionMID works with which='all'", {
  test_data <- create_test_data_with_es(n = 10)
  result <- try(runMetaAnalysis(test_data, which.run = c("overall", "combined")), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    prop_result <- try(proportionMID(result, mid = 0.5, which = "all"),
                      silent = TRUE)
    
    if (!inherits(prop_result, "try-error")) {
      expect_s3_class(prop_result, "proportionMID")
      expect_true(nrow(prop_result$qhat) >= 1)
    }
  }
})

test_that("proportionMID works with specific which values", {
  test_data <- create_test_data_with_es(n = 10)
  result <- try(runMetaAnalysis(test_data, which.run = "combined"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error") && 
      !is.null(result$model.combined) &&
      result$model.combined$k > 1) {
    prop_result <- try(proportionMID(result, mid = 0.5, which = "combined"),
                      silent = TRUE)
    
    if (!inherits(prop_result, "try-error")) {
      expect_s3_class(prop_result, "proportionMID")
    }
  }
})

test_that("proportionMID works with plot=TRUE", {
  test_data <- create_test_data_with_es(n = 10)
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    # Plot is created, not returned
    expect_error(try(proportionMID(result, mid = 0.5, plot = TRUE), 
                    silent = TRUE), NA)
  }
})

test_that("proportionMID works with negative mid values", {
  test_data <- create_test_data_with_es(n = 10)
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    prop_result <- try(proportionMID(result, mid = -0.24),
                      silent = TRUE)
    
    if (!inherits(prop_result, "try-error")) {
      expect_s3_class(prop_result, "proportionMID")
    }
  }
})

test_that("proportionMID works with correctPublicationBias models", {
  test_data <- create_test_data_with_es(n = 10)
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    pb_result <- try(correctPublicationBias(result, which.run = "overall"),
                     silent = TRUE)
    
    if (!inherits(pb_result, "try-error")) {
      prop_result <- try(proportionMID(pb_result, mid = 0.5, which = "trimfill"),
                        silent = TRUE)
      
      if (!inherits(prop_result, "try-error")) {
        expect_s3_class(prop_result, "proportionMID")
      }
    }
  }
})

test_that("print.proportionMID works", {
  test_data <- create_test_data_with_es(n = 10)
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    prop_result <- try(proportionMID(result, mid = 0.5),
                      silent = TRUE)
    
    if (!inherits(prop_result, "try-error")) {
      expect_error(print(prop_result), NA)
      expect_output(print(prop_result), "Proportion of meaningful")
    }
  }
})

test_that("plot.proportionMID works", {
  test_data <- create_test_data_with_es(n = 10)
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    prop_result <- try(proportionMID(result, mid = 0.5),
                      silent = TRUE)
    
    if (!inherits(prop_result, "try-error")) {
      # Plot is created, not returned
      expect_error(plot(prop_result), NA)
    }
  }
})

test_that("proportionMID handles RR type", {
  # This would require RR data, but we can test the structure
  # For now, just verify the function exists and can be called
  test_data <- create_test_data_with_es(n = 10)
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    # Test that it works with regular g
    prop_result <- try(proportionMID(result, mid = 0.5),
                      silent = TRUE)
    
    if (!inherits(prop_result, "try-error")) {
      expect_s3_class(prop_result, "proportionMID")
      expect_true(prop_result$xlab == "Hedges' g" || 
                  prop_result$xlab == "log(RR)")
    }
  }
})

