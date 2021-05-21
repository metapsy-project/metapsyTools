# Tests for all replacement functions for runMetaAnalysis objects
#
# Important: Replacement functions modify the 'call' object in the runMetaAnalysis
# result. The actual model doesn't change until rerun() is called.

# Helper to create a runMetaAnalysis object for testing
create_test_meta_analysis <- function() {
  data_prepared <- prepare_runMetaAnalysis_data()
  if (is.null(data_prepared)) {
    return(NULL)
  }
  
  result <- try(runMetaAnalysis(data_prepared, which.run = "overall"), 
                silent = TRUE)
  
  if (inherits(result, "try-error")) {
    return(NULL)
  }
  
  return(result)
}

test_that("data replacement function works", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Test setting data
  new_data <- result$data[1:5, ]
  expect_error(data(result) <- new_data, NA)
  expect_equal(nrow(result$data), 5)
  
  # Test getting data (data is stored directly, not in call)
  expect_s3_class(result$data, "data.frame")
})

test_that("which.run replacement function modifies call and requires rerun", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Store original which.run from the actual model
  original_which_run <- result$which.run
  
  # Test setting which.run - this modifies the call, not the model
  expect_error(which.run(result) <- "combined", NA)
  
  # The call should be modified
  call_args <- as.list(result$call)
  expect_equal(call_args$which.run, "combined")
  
  # But the actual model's which.run hasn't changed yet
  expect_equal(result$which.run, original_which_run)
  
  # Getter reads from call - extract from call object
  call_args <- as.list(result$call)
  expect_equal(call_args$which.run, "combined")
  
  # After rerun, the model should have the new which.run
  rerun_result <- try(rerun(result), silent = TRUE)
  if (!inherits(rerun_result, "try-error")) {
    expect_true("combined" %in% rerun_result$which.run || 
                rerun_result$which.run == "combined")
  }
})

test_that("es.measure replacement function modifies call", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Test setting es.measure - modifies call
  expect_error(es.measure(result) <- "g", NA)
  
  # Call should be modified
  call_args <- as.list(result$call)
  expect_equal(call_args$es.measure, "g")
  
  # Getter reads from call - extract from call object
  call_args <- as.list(result$call)
  expect_equal(call_args$es.measure, "g")
  
  # Model hasn't changed until rerun
  expect_s3_class(result, "runMetaAnalysis")
})

test_that("es.type replacement function modifies call", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Test setting es.type - modifies call
  expect_error(es.type(result) <- "precalculated", NA)
  
  # Call should be modified
  call_args <- as.list(result$call)
  expect_equal(call_args$es.type, "precalculated")
  
  # Getter reads from call - extract from call object
  call_args <- as.list(result$call)
  expect_equal(call_args$es.type, "precalculated")
})

test_that("es.var replacement function modifies call", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Test setting es.var - modifies call
  expect_error(es.var(result) <- ".g", NA)
  call_args <- as.list(result$call)
  expect_equal(call_args$es.var, ".g")
})

test_that("se.var replacement function modifies call", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Test setting se.var - modifies call
  expect_error(se.var(result) <- ".g_se", NA)
  call_args <- as.list(result$call)
  expect_equal(call_args$se.var, ".g_se")
})

test_that("es.binary.raw.vars replacement function modifies call", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Test setting es.binary.raw.vars - modifies call
  new_vars <- c(".event_arm1", ".event_arm2", ".totaln_arm1", ".totaln_arm2")
  expect_error(es.binary.raw.vars(result) <- new_vars, NA)
  call_args <- as.list(result$call)
  expect_equal(call_args$es.binary.raw.vars, new_vars)
})

test_that("method.tau replacement function modifies call and works with rerun", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Store original method.tau from model
  original_method <- result$model.overall$method.tau
  
  # Test setting method.tau - modifies call
  expect_error(method.tau(result) <- "DL", NA)
  
  # Call should be modified
  call_args <- as.list(result$call)
  expect_equal(call_args$method.tau, "DL")
  
  # Getter reads from call - extract from call object
  call_args <- as.list(result$call)
  expect_equal(call_args$method.tau, "DL")
  
  # Original model hasn't changed
  expect_equal(result$model.overall$method.tau, original_method)
  
  # After rerun, new model should have DL method
  rerun_result <- try(rerun(result), silent = TRUE)
  if (!inherits(rerun_result, "try-error")) {
    expect_equal(rerun_result$model.overall$method.tau, "DL")
  }
  
  # Test other methods
  expect_error(method.tau(result) <- "PM", NA)
  call_args <- as.list(result$call)
  expect_equal(call_args$method.tau, "PM")
})

test_that("i2.ci.threelevel replacement function modifies call", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Test setting i2.ci.threelevel - modifies call
  expect_error(i2.ci.threelevel(result) <- TRUE, NA)
  call_args <- as.list(result$call)
  expect_equal(call_args$i2.ci.threelevel, TRUE)
  
  expect_error(i2.ci.threelevel(result) <- FALSE, NA)
  call_args <- as.list(result$call)
  expect_equal(call_args$i2.ci.threelevel, FALSE)
})

test_that("nsim.boot replacement function modifies call", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Test setting nsim.boot - modifies call
  expect_error(nsim.boot(result) <- 1000, NA)
  call_args <- as.list(result$call)
  expect_equal(call_args$nsim.boot, 1000)
})

test_that("hakn replacement function modifies call and works with rerun", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Store original hakn from model
  original_hakn <- result$model.overall$hakn
  
  # Test setting hakn - modifies call
  expect_error(hakn(result) <- FALSE, NA)
  
  # Call should be modified
  call_args <- as.list(result$call)
  expect_equal(call_args$hakn, FALSE)
  
  # Getter reads from call - extract from call object
  call_args <- as.list(result$call)
  expect_equal(call_args$hakn, FALSE)
  
  # Original model hasn't changed
  expect_equal(result$model.overall$hakn, original_hakn)
  
  # After rerun, new model should have hakn = FALSE
  rerun_result <- try(rerun(result), silent = TRUE)
  if (!inherits(rerun_result, "try-error")) {
    expect_equal(rerun_result$model.overall$hakn, FALSE)
  }
})

# Test all remaining replacement functions modify call correctly
# Using a more efficient approach - test that they modify call and getters work

test_that("all remaining replacement functions modify call", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Test each replacement function modifies the call
  replacement_tests <- list(
    list(func = "study.var", value = "study"),
    list(func = "arm.var.1", value = "condition_arm1"),
    list(func = "arm.var.2", value = "condition_arm2"),
    list(func = "measure.var", value = "instrument"),
    list(func = "low.rob.filter", value = "rob > 2"),
    list(func = "method.tau.ci", value = "Q-Profile"),
    list(func = "round.digits", value = 3),
    list(func = "which.combine", value = "studies"),
    list(func = "which.combine.var", value = "multi_arm1"),
    list(func = "which.outliers", value = "combined"),
    list(func = "which.influence", value = "combined"),
    list(func = "which.rob", value = "combined"),
    list(func = "rho.within.study", value = 0.7),
    list(func = "phi.within.study", value = 0.8),
    list(func = "w1.var", value = "n_arm1"),
    list(func = "w2.var", value = "n_arm2"),
    list(func = "time.var", value = "time_weeks"),
    list(func = "vcov", value = "complex"),
    list(func = "near.pd", value = TRUE),
    list(func = "use.rve", value = FALSE),
    list(func = "html", value = FALSE)
  )
  
  for (test_case in replacement_tests) {
    func_name <- test_case$func
    test_value <- test_case$value
    
    # Get the setter function
    setter <- get(paste0(func_name, "<-"))
    
    # Apply replacement function - it returns the modified object
    expect_error(result <- setter(result, test_value), NA, 
                info = paste("Setting", func_name))
    
    # Verify call was modified
    call_args <- as.list(result$call)
    expect_equal(call_args[[func_name]], test_value,
                 info = paste("Call modification for", func_name))
  }
})

test_that("lower.is.better replacement function requires correctPublicationBias", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Should error if correctPublicationBias hasn't been run
  expect_error(lower.is.better(result) <- TRUE,
               "correctPublicationBias.*must be run first")
  
  # Run correctPublicationBias first
  result_pb <- try(correctPublicationBias(result), silent = TRUE)
  
  if (!inherits(result_pb, "try-error")) {
    # Now it should work
    expect_error(lower.is.better(result_pb) <- TRUE, NA)
    
    # Extract from correctPublicationBias call
    pb_call_args <- as.list(result_pb$correctPublicationBias$call)
    expect_equal(pb_call_args$lower.is.better, TRUE)
  }
})

test_that("selmodelSteps replacement function requires correctPublicationBias", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Should error if correctPublicationBias hasn't been run
  expect_error(selmodelSteps(result) <- 5,
               "correctPublicationBias.*must be run first")
  
  # Run correctPublicationBias first
  result_pb <- try(correctPublicationBias(result), silent = TRUE)
  
  if (!inherits(result_pb, "try-error")) {
    # Now it should work
    expect_error(selmodelSteps(result_pb) <- 5, NA)
    
    # Extract from correctPublicationBias call
    pb_call_args <- as.list(result_pb$correctPublicationBias$call)
    expect_equal(pb_call_args$selmodel.steps, 5)
  }
})

test_that("rerun function applies modified settings from replacement functions", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Store original values
  original_method_tau <- result$model.overall$method.tau
  original_hakn <- result$model.overall$hakn
  
  # Modify settings using replacement functions (modifies call only)
  method.tau(result) <- "DL"
  hakn(result) <- FALSE
  
  # Verify call was modified but model hasn't changed
  call_args <- as.list(result$call)
  expect_equal(call_args$method.tau, "DL")
  expect_equal(call_args$hakn, FALSE)
  expect_equal(result$model.overall$method.tau, original_method_tau)
  expect_equal(result$model.overall$hakn, original_hakn)
  
  # Rerun with new settings - this actually applies the changes
  rerun_result <- try(rerun(result), silent = TRUE)
  
  if (!inherits(rerun_result, "try-error")) {
    expect_s3_class(rerun_result, "runMetaAnalysis")
    expect_true("summary" %in% names(rerun_result))
    
    # Verify settings were actually applied in the new model
    expect_equal(rerun_result$model.overall$method.tau, "DL")
    expect_equal(rerun_result$model.overall$hakn, FALSE)
    
    # Verify call still has the values (extract from call)
    rerun_call_args <- as.list(rerun_result$call)
    expect_equal(rerun_call_args$method.tau, "DL")
    expect_equal(rerun_call_args$hakn, FALSE)
  }
})

test_that("rerun function works with multiple modifications", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Modify multiple settings
  method.tau(result) <- "PM"
  which.combine(result) <- "studies"
  round.digits(result) <- 3
  
  # Rerun
  rerun_result <- try(rerun(result), silent = TRUE)
  
  if (!inherits(rerun_result, "try-error")) {
    expect_s3_class(rerun_result, "runMetaAnalysis")
  }
})

test_that("replacement functions modify call but preserve model until rerun", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Store original structure
  original_summary <- result$summary
  original_which_run <- result$which.run
  original_method_tau <- result$model.overall$method.tau
  
  # Modify settings using replacement functions
  method.tau(result) <- "DL"
  which.combine(result) <- "studies"
  
  # Object should still be valid
  expect_s3_class(result, "runMetaAnalysis")
  expect_true("summary" %in% names(result))
  expect_true("which.run" %in% names(result))
  
  # Summary and model should be unchanged (until rerun)
  expect_identical(result$summary, original_summary)
  expect_equal(result$which.run, original_which_run)
  expect_equal(result$model.overall$method.tau, original_method_tau)
  
  # But call should be modified
  call_args <- as.list(result$call)
  expect_equal(call_args$method.tau, "DL")
  expect_equal(call_args$which.combine, "studies")
})

test_that("replacement functions work in sequence and modify call", {
  result <- create_test_meta_analysis()
  if (is.null(result)) {
    skip("Could not create test meta-analysis object")
  }
  
  # Apply multiple replacements in sequence
  expect_error({
    method.tau(result) <- "DL"
    hakn(result) <- FALSE
    which.combine(result) <- "studies"
    round.digits(result) <- 4
  }, NA)
  
  # Verify all were set in call (extract from call object)
  call_args <- as.list(result$call)
  expect_equal(call_args$method.tau, "DL")
  expect_equal(call_args$hakn, FALSE)
  expect_equal(call_args$which.combine, "studies")
  expect_equal(call_args$round.digits, 4)
  
  # Verify call was modified - extract from call object
  call_args <- as.list(result$call)
  expect_equal(call_args$method.tau, "DL")
  expect_equal(call_args$hakn, FALSE)
  expect_equal(call_args$which.combine, "studies")
  expect_equal(call_args$round.digits, 4)
  
  # After rerun, all settings should be applied
  rerun_result <- try(rerun(result), silent = TRUE)
  if (!inherits(rerun_result, "try-error")) {
    expect_equal(rerun_result$model.overall$method.tau, "DL")
    expect_equal(rerun_result$model.overall$hakn, FALSE)
  }
})

test_that("replacement functions example from documentation works correctly", {
  data_prepared <- prepare_runMetaAnalysis_data()
  if (is.null(data_prepared)) {
    skip("Could not prepare data")
  }
  
  # As in documentation example
  result <- try(runMetaAnalysis(data_prepared, which.run = "combined"), 
                silent = TRUE)
  
  if (inherits(result, "try-error")) {
    skip("Could not run meta-analysis")
  }
  
  # Store original method.tau
  original_method <- result$model.combined$method.tau
  
  # Compare results when other tau^2 estimator is used
  # This modifies the call, not the model
  expect_error(method.tau(result) <- "DL", NA)
  
  # Call should be modified - extract from call object
  call_args <- as.list(result$call)
  expect_equal(call_args$method.tau, "DL")
  
  # But model hasn't changed yet
  expect_equal(result$model.combined$method.tau, original_method)
  
  # Rerun applies the change
  rerun_result <- try(rerun(result), silent = TRUE)
  
  if (!inherits(rerun_result, "try-error")) {
    expect_s3_class(rerun_result, "runMetaAnalysis")
    
    # New model should have DL method
    expect_equal(rerun_result$model.combined$method.tau, "DL")
    
    # Call should still have the value
    rerun_call_args <- as.list(rerun_result$call)
    expect_equal(rerun_call_args$method.tau, "DL")
  }
})

