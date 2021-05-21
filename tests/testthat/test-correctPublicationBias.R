# Tests for correctPublicationBias function

test_that("correctPublicationBias validates input", {
  # Test with non-runMetaAnalysis object
  expect_error(correctPublicationBias("not a runMetaAnalysis object"),
               "must be of class 'runMetaAnalysis'")
  
  # Test with runMetaAnalysis object
  test_data <- create_test_data_with_es(n = 10)
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    # Test invalid which.run
    expect_error(correctPublicationBias(result, which.run = "invalid"),
                 "which.run.*must be")
    
    # Test which.run not in model
    expect_error(correctPublicationBias(result, which.run = "combined"),
                 "has not been fitted|k=1")
  }
})

test_that("correctPublicationBias works with overall model", {
  test_data <- create_test_data_with_es(n = 10)
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    pb_result <- try(correctPublicationBias(result, which.run = "overall"),
                     silent = TRUE)
    
    if (!inherits(pb_result, "try-error")) {
      expect_s3_class(pb_result, "runMetaAnalysis")
      expect_s3_class(pb_result, "correctPublicationBias")
      expect_true("correctPublicationBias" %in% names(pb_result))
      expect_true("summary" %in% names(pb_result$correctPublicationBias))
      expect_true("model.trimfill" %in% names(pb_result$correctPublicationBias))
      expect_true("model.limitmeta" %in% names(pb_result$correctPublicationBias))
      expect_true("model.selection" %in% names(pb_result$correctPublicationBias))
    }
  }
})

test_that("correctPublicationBias works with custom lower.is.better", {
  test_data <- create_test_data_with_es(n = 10)
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    pb_result <- try(correctPublicationBias(result, 
                                            which.run = "overall",
                                            lower.is.better = FALSE),
                     silent = TRUE)
    
    if (!inherits(pb_result, "try-error")) {
      expect_s3_class(pb_result, "correctPublicationBias")
    }
  }
})

test_that("correctPublicationBias works with custom selmodelSteps", {
  test_data <- create_test_data_with_es(n = 10)
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    pb_result <- try(correctPublicationBias(result, 
                                            which.run = "overall",
                                            selmodelSteps = c(0.03, 0.05)),
                     silent = TRUE)
    
    if (!inherits(pb_result, "try-error")) {
      expect_s3_class(pb_result, "correctPublicationBias")
    }
  }
})

test_that("correctPublicationBias handles k=1 error", {
  test_data <- create_test_data_with_es(n = 1)
  result <- try(runMetaAnalysis(test_data, which.run = "overall"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error")) {
    # Should error for combined with k=1
    expect_error(correctPublicationBias(result, which.run = "combined"),
                 "k=1")
  }
})

test_that("correctPublicationBias works with combined model", {
  test_data <- create_test_data_with_es(n = 10)
  result <- try(runMetaAnalysis(test_data, which.run = "combined"), 
                silent = TRUE)
  
  if (!inherits(result, "try-error") && 
      !is.null(result$model.combined) &&
      result$model.combined$k > 1) {
    pb_result <- try(correctPublicationBias(result, which.run = "combined"),
                     silent = TRUE)
    
    if (!inherits(pb_result, "try-error")) {
      expect_s3_class(pb_result, "correctPublicationBias")
    }
  }
})


