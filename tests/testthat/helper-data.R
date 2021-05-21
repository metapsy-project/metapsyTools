# Helper functions and test data for metapsyTools tests

# Create minimal valid test dataset
create_test_data <- function() {
  data.frame(
    study = c("Study1", "Study2", "Study3"),
    condition_arm1 = c("Treatment", "Treatment", "Treatment"),
    condition_arm2 = c("Control", "Control", "Control"),
    multi_arm1 = c("A", "A", "A"),
    multi_arm2 = c("B", "B", "B"),
    outcome_type = c("continuous", "continuous", "continuous"),
    instrument = c("PHQ-9", "PHQ-9", "PHQ-9"),
    time = c("post", "post", "post"),
    time_weeks = c(8, 8, 8),
    rating = c("self", "self", "self"),
    mean_arm1 = c(10.5, 11.2, 9.8),
    mean_arm2 = c(15.3, 16.1, 14.9),
    sd_arm1 = c(3.2, 3.5, 2.9),
    sd_arm2 = c(4.1, 4.3, 3.8),
    n_arm1 = c(50, 45, 55),
    n_arm2 = c(48, 47, 52),
    event_arm1 = c(NA, NA, NA),
    event_arm2 = c(NA, NA, NA),
    totaln_arm1 = c(NA, NA, NA),
    totaln_arm2 = c(NA, NA, NA),
    stringsAsFactors = FALSE
  )
}

# Create test data with binary outcomes
create_test_data_binary <- function() {
  data.frame(
    study = c("Study1", "Study2"),
    condition_arm1 = c("Treatment", "Treatment"),
    condition_arm2 = c("Control", "Control"),
    multi_arm1 = c("A", "A"),
    multi_arm2 = c("B", "B"),
    outcome_type = c("dichotomous", "dichotomous"),
    instrument = c("Response", "Response"),
    time = c("post", "post"),
    time_weeks = c(8, 8),
    rating = c("self", "self"),
    mean_arm1 = c(NA, NA),
    mean_arm2 = c(NA, NA),
    sd_arm1 = c(NA, NA),
    sd_arm2 = c(NA, NA),
    n_arm1 = c(NA, NA),
    n_arm2 = c(NA, NA),
    event_arm1 = c(25, 30),
    event_arm2 = c(15, 18),
    totaln_arm1 = c(50, 60),
    totaln_arm2 = c(50, 60),
    stringsAsFactors = FALSE
  )
}

# Create test data with pre-calculated effect sizes (for runMetaAnalysis tests)
create_test_data_with_es <- function(n = 5) {
  data.frame(
    study = paste0("Study", 1:n),
    condition_arm1 = rep("Treatment", n),
    condition_arm2 = rep("Control", n),
    multi_arm1 = rep("A", n),
    multi_arm2 = rep("B", n),
    outcome_type = rep("continuous", n),
    instrument = rep("PHQ-9", n),
    time = rep("post", n),
    time_weeks = rep(8, n),
    rating = rep("self", n),
    .g = rnorm(n, mean = 0.5, sd = 0.2),  # Pre-calculated effect sizes
    .g_se = runif(n, 0.1, 0.3),          # Pre-calculated standard errors
    n_arm1 = rep(50, n),
    n_arm2 = rep(48, n),
    stringsAsFactors = FALSE
  )
}

# Helper function to prepare data for runMetaAnalysis examples
prepare_runMetaAnalysis_data <- function() {
  data("depressionPsyCtr", package = "metapsyTools", envir = environment())
  
  data_formatted <- checkDataFormat(depressionPsyCtr)
  data_checked <- checkConflicts(data_formatted)
  
  if (inherits(data_checked, "checkConflicts")) {
    return(NULL)
  }
  
  data_es <- try(calculateEffectSizes(data_checked), silent = TRUE)
  if (inherits(data_es, "try-error")) {
    return(NULL)
  }
  
  data_filtered <- try(filterPoolingData(data_es, 
                                         condition_arm2 %in% c("wl", "other ctr")),
                       silent = TRUE)
  
  if (inherits(data_filtered, "try-error")) {
    return(NULL)
  }
  
  return(data_filtered)
}

