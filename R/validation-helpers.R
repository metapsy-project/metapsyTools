# Validate that data is a data.frame (internal helper)
validate_data_frame <- function(data, arg_name = "data") {
  if (!is.data.frame(data)) {
    stop(sprintf("'%s' must be a data.frame, but got class: %s", 
                 arg_name, paste(class(data), collapse = ", ")),
         call. = FALSE)
  }
  if (nrow(data) == 0) {
    stop(sprintf("'%s' must have at least one row", arg_name),
         call. = FALSE)
  }
  invisible(TRUE)
}

# Validate that required columns exist (internal helper)
validate_required_columns <- function(data, required_cols, arg_name = "data") {
  if (length(required_cols) == 0) {
    return(character(0))
  }
  
  missing <- required_cols[!required_cols %in% colnames(data)]
  return(missing)
}

# Safely convert variable to target class (internal helper)
safe_convert_class <- function(x, target_class, var_name) {
  current_class <- class(x)[1]
  
  if (current_class == target_class) {
    return(x)
  }
  
  # Try conversion
  result <- tryCatch({
    methods::as(x, target_class)
  }, error = function(e) {
    warning(sprintf("Could not convert '%s' from %s to %s: %s",
                    var_name, current_class, target_class, e$message),
            call. = FALSE)
    return(x)
  })
  
  return(result)
}

# Check if column name is reserved (internal helper)
is_reserved_name <- function(col_name) {
  reserved <- c("subset", "exclude")
  col_name %in% reserved
}

# Rename reserved column names (internal helper)
rename_reserved_columns <- function(data) {
  reserved <- c("subset", "exclude")
  renamed <- FALSE
  
  for (reserved_name in reserved) {
    if (reserved_name %in% colnames(data)) {
      new_name <- paste0(reserved_name, ".1")
      data[[new_name]] <- data[[reserved_name]]
      data[[reserved_name]] <- NULL
      renamed <- TRUE
    }
  }
  
  return(list(data = data, renamed = renamed))
}

