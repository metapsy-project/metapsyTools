test_that("validate_data_frame works correctly", {
  # Valid data.frame
  df <- data.frame(x = 1:3)
  expect_invisible(validate_data_frame(df))
  
  # Not a data.frame
  expect_error(validate_data_frame(c(1, 2, 3)), 
               "must be a data.frame")
  
  # Empty data.frame
  empty_df <- data.frame()
  expect_error(validate_data_frame(empty_df), 
               "must have at least one row")
})

test_that("validate_required_columns works correctly", {
  df <- data.frame(a = 1, b = 2, c = 3)
  
  # All columns present
  expect_equal(validate_required_columns(df, c("a", "b")), character(0))
  
  # Some missing
  missing <- validate_required_columns(df, c("a", "b", "d", "e"))
  expect_equal(sort(missing), c("d", "e"))
  
  # Empty required
  expect_equal(validate_required_columns(df, character(0)), character(0))
})

test_that("safe_convert_class works correctly", {
  # Already correct class
  x <- c(1, 2, 3)
  result <- safe_convert_class(x, "numeric", "test_var")
  expect_type(result, "double")
  
  # Convertible
  x_char <- c("1", "2", "3")
  result <- safe_convert_class(x_char, "numeric", "test_var")
  expect_type(result, "double")
  
  # Not convertible (should return original with warning)
  # R's as.numeric() gives "NAs introduced by coercion" warning
  x_invalid <- c("a", "b", "c")
  expect_warning(
    result <- safe_convert_class(x_invalid, "numeric", "test_var"),
    "NAs introduced|Could not convert"
  )
  # Result should have NAs
  expect_true(any(is.na(result)))
})

test_that("is_reserved_name works correctly", {
  expect_true(is_reserved_name("subset"))
  expect_true(is_reserved_name("exclude"))
  expect_false(is_reserved_name("study"))
  expect_false(is_reserved_name("other"))
})

test_that("rename_reserved_columns works correctly", {
  # No reserved names
  df <- data.frame(study = c("A", "B"), value = 1:2)
  result <- rename_reserved_columns(df)
  expect_equal(result$data, df)
  expect_false(result$renamed)
  
  # With reserved names
  df_reserved <- data.frame(study = c("A", "B"), subset = c(1, 2))
  result <- rename_reserved_columns(df_reserved)
  expect_false("subset" %in% names(result$data))
  expect_true("subset.1" %in% names(result$data))
  expect_true(result$renamed)
  
  # Multiple reserved names
  df_multi <- data.frame(study = c("A", "B"), 
                         subset = c(1, 2), 
                         exclude = c(3, 4))
  result <- rename_reserved_columns(df_multi)
  expect_false("subset" %in% names(result$data))
  expect_false("exclude" %in% names(result$data))
  expect_true("subset.1" %in% names(result$data))
  expect_true("exclude.1" %in% names(result$data))
})

