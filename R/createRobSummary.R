#' Create a summary risk of bias plot
#'
#' If the `rob.data` argument has been specified, this function allows to create
#' a summary risk of bias plot for results of the [runMetaAnalysis()] function.
#'
#' @usage createRobSummary(model, 
#'                  name.low, 
#'                  name.high, 
#'                  name.unclear, 
#'                  which.run = model$which.run[1])
#'
#' @param model An object of class \code{runMetaAnalysis}, created by the [runMetaAnalysis()] function.
#' @param name.low A `character` vector, specifying which code(s) have been used in the original data for 
#' studies with a low risk of bias.
#' @param name.high A `character` vector, specifying which code(s) have been used in the original data for 
#' studies with a high risk of bias.
#' @param name.unclear A `character` vector, specifying which code(s) have been used in the original data for 
#' studies with unclear risk of bias.
#' @param which.run The model in \code{model} that should be used for the summary risk of bias plot. 
#' Uses the default analysis in \code{model} if no value is specified by the user. Possible values are
#' \code{"overall"}, \code{"combined"}, \code{"lowest"}, \code{"highest"}, \code{"outliers"},
#' \code{"influence"} and \code{"rob"}.
#'
#' @return Creates a RevMan-type risk of bias summary plot.
#'
#'
#' @examples
#' \dontrun{
#' 
#' 
#' # Define ROB data to be added to the models
#' robData = list(
#'   # Names of ROB variables included in 'data'
#'   domains = c("sg", "ac", "ba", "itt"),
#'   # Long-format labels for each ROB domain
#'   domain.names = c("Sequence Generation", 
#'                    "Allocation Concealment", 
#'                    "Blinding of Assessors", 
#'                    "ITT Analyses"),
#'   # Codes used to rate the risk of bias (sr=self-report)
#'   categories = c("0", "1", "sr"),
#'   # Symbols that should be used for these codes in forest plots
#'   symbols = c("-", "+", "s"),
#'   # Colors to be used in forest plots for each of these codes
#'   colors = c("red", "green", "yellow"))
#' 
#' # Run meta-analyses with ROB data
#' res <- depressionPsyCtr %>% 
#'   filterPoolingData(condition_arm1 %in% c("cbt", "pst", "3rd")) %>% 
#'   runMetaAnalysis(rob.data = robData)
#' 
#' # Create a summary plot
#' createRobSummary(res, 
#'                  name.low = "1", 
#'                  name.high = "0", 
#'                  name.unclear = "sr")
#' 
#' # Create a summary plot for the "combined" model
#' # - Recode 'sr' (self-report) as low risk of bias
#' createRobSummary(res, 
#'                  name.low = c("1", "sr"), 
#'                  name.high = "0", 
#'                  name.unclear = NULL,
#'                  which.run = "combined")
#'                  
#'                  
#' }
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso [dmetar::rob.summary], [robvis::rob_summary], [meta::rob]
#' 
#' @importFrom dplyr recode
#' @importFrom stringr str_remove_all
#' @export createRobSummary

createRobSummary = function(model, name.low, name.high, name.unclear, 
                            which.run = model$which.run[1]) {
  
  # Check class
  if (class(model)[1] != "runMetaAnalysis"){
    stop("Input must be of class 'runMetaAnalysis'. Did you apply 'runMetaAnalysis' first?")
  }
  
  # Throw an error if which.run does not contain required model
  if (!(which.run[1] %in% c("overall", "combined",
                            "outliers", "rob", "highest",
                            "lowest", "influence"))){
    stop("'which.run' must be 'overall', 'combined', 'outliers', 'rob',",
         " 'highest', 'lowest' or 'influence'.")
  }
  if (!(which.run[1] %in% model$which.run)){
    stop("'", which.run[1], "' model has not been fitted using runMetaAnalysis.")
  }
  if (which.run[1] %in% unlist(model$error.model.list)){
    stop("There was an error fitting the '", which.run[1], "' model using runMetaAnalysis;",
         " the correctPublicationBias function cannot be applied.")
  }
  
  # Get model type
  model.type = paste0("model.", which.run)
  M = model[[model.type]]
  
  # Create recoder for ROB plot
  unlist(list("High" = name.high, "Low" = name.low, 
              "Unclear" = name.unclear)) -> recoder
  stringr::str_remove_all(names(recoder), "[0-9]") %>% 
    {names(.) = recoder;.} -> recoder
  
  # Extract rob element
  if (!class(M$rob)[1] == "rob") {
    stop("Model does not contain ROB data.",
         " Specify 'rob.data' argument when running 'runMetaAnalysis'.")
  }
  message("- ", crayon::green("[OK] "), "'", model.type, 
          "' used to generate summary ROB plot.")
  M$rob -> rob.elem
  rob.df = data.frame(rob.elem)[-c(1, ncol(rob.elem))]
  rob.domains = attr(rob.elem, "domains")
  colnames(rob.df) = rob.domains
  
  # Recode the judgments
  lapply(as.list(rob.df), function(x) { dplyr::recode(x, !!!recoder) }) %>% 
    as.data.frame() -> rob.df
  
  # Create plot
  p = try(robSummary(rob.df), silent=TRUE)
  invisible(p)
  
}






