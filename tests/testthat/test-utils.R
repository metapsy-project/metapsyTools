test_that("isNAorNaN works correctly", {
  expect_true(isNAorNaN(NA))
  expect_true(isNAorNaN(NaN))
  expect_false(isNAorNaN(1))
  expect_false(isNAorNaN("a"))
  # NULL returns logical(0) from is.na(), so the result is logical(0), not FALSE
  # This is technically correct behavior - NULL is not NA or NaN
  # Check that any() of the result is FALSE (meaning no TRUE values)
  result_null <- isNAorNaN(NULL)
  expect_false(any(result_null))
  
  # Vectorized
  x <- c(1, NA, NaN, 2)
  result <- isNAorNaN(x)
  expect_equal(result, c(FALSE, TRUE, TRUE, FALSE))
})

test_that("tryCatch2 captures errors and warnings", {
  # Successful execution
  result <- tryCatch2(1 + 1)
  expect_equal(result$value, 2)
  expect_null(result$warning)
  expect_null(result$error)
  expect_false(result$has.error)
  
  # Error case
  result <- tryCatch2(stop("test error"))
  expect_null(result$value)
  expect_null(result$warning)
  expect_true(!is.null(result$error))
  expect_true(result$has.error)
  
  # Warning case - warning() returns NULL but the expression is evaluated
  # The value is the result of the expression, not NULL
  result <- tryCatch2({warning("test warning"); "test_value"})
  expect_equal(result$value, "test_value")
  expect_true(!is.null(result$warning))
  expect_null(result$error)
  expect_false(result$has.error)
})

test_that("selectArguments filters function arguments correctly", {
  test_fun <- function(a, b, c) { a + b + c }
  dots <- list(a = 1, b = 2, d = 3, e = 4)
  
  result <- selectArguments(test_fun, dots)
  expect_equal(result, list(a = 1, b = 2))
  expect_false("d" %in% names(result))
  expect_false("e" %in% names(result))
})

