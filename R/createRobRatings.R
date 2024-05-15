#' Generate RoB rating for a Metapsy database
#'
#' Based on a prepared extraction sheet, this function allows you 
#' to automatically generate risk of bias judgments for an entire
#' Metapsy database, using the RoB-2 Metapsy extension.
#'
#' @usage createRobRatings(database, rob.data)
#'
#' @param database A meta-analytic database conforming with the Metapsy data standard. In particular, the
#' variables `rand_ratio`, `rand_arm1`/`rand_arm2`, `attr_arm1`/`attr_arm2` and `rating` must be populated
#' according to the Metapsy [data standard](https://docs.metapsy.org/data-preparation/format/). 
#' @param rob.data An RoB extraction sheet. Columns of this file must have the same name as in the
#' [RoB extraction sheet template](https://www.metapsy.org/assets/files/rob-template.xlsx) 
#' provided by the Metapsy initative. If the Metapsy template has been
#' used, make sure to delete the top rows before importing, so that only the metapsyTools variables
#' remain as the column names. Required columns are: `study`, `d1_1`, `d1_2`, `d1_3`, `d1_4`, `d1_notes`, 
#' `d2_5`, `d2_6`, `d2_7`, `d2_8`, `d2_9`, `d2_notes`, `d3_10`, `d3_11`, `d3_12`, `d3_13`, 
#' `d3_14`, `d3_notes`, `d4_15`, `d4_16`, `d4_17`, `d4_18`, `d4_notes`, `d5_19`, 
#' `d5_20`, `d5_21`, `d5_22`, `d5_23`, `d5_24`, `d5_notes`. 
#'
#' @return Returns an object of class \code{"createRobRatings"}. This includes: 
#' 
#' - `database`. The database with the RoB ratings appended on a domain level: 
#' `rob_d1`, `rob_d2`, etc. This dataset also includes the comparison-level (`rob`) and
#' study-level (`rob_study_lvl`) overall RoB judgments.
#' - `rob.data`. The original extraction sheet, with signalling question responses 
#' based on the database information appended. If a study has more than one comparison 
#' in the database, ratings are displayed for each comparison in a study.
#' - `miss.studies`. A vector of studies that appear in the database, but were not
#' found in the extraction sheet. 
#'
#' @examples
#' \dontrun{
#' library(readxl)
#' 
#' # Get example database from metapsy.org/assets/files/data.xlsx
#' data <- read_excel("data.xlsx")
#' 
#' # Get example extraction sheet from metapsy.org/assets/files/rob_data.xlsx
#' rob_data <- read_excel("rob_data.xlsx")
#' 
#' # Create ratings
#' tmp <- metapsyTools:::createRobRatings(database = data, 
#'                                        rob.data = rob_data)
#' tmp
#' 
#' # Show database
#' tmp$database
#' 
#' # Show extraction sheet
#' tmp$rob.data
#' }
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Clara Miguel Sanz \email{clara.miguelsanz@@vu.nl}, 
#' Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @import dplyr
#' @importFrom crayon green yellow cyan bold
#' @importFrom stringr str_replace_all
#' @keywords internal

createRobRatings = function(database, rob.data) {
  
  
  # Define required variables
  data = database
  req.variables = c("study", "d1_1", "d1_2", "d1_3", "d1_4", "d1_notes", "d2_5", 
                    "d2_6", "d2_7", "d2_8", "d2_9", "d2_notes", "d3_10", "d3_11", 
                    "d3_12", "d3_13", "d3_14", "d3_notes", "d4_15", "d4_16", 
                    "d4_17", "d4_18", "d4_notes", "d5_19", "d5_20", "d5_21", 
                    "d5_22", "d5_23", "d5_24", "d5_notes")
  
  # Check if required variables are included
  if (sum(!req.variables %in% colnames(rob.data)) > 0) {
    stop("Required variables ", 
         paste(req.variables[!req.variables %in% colnames(rob.data)], 
               collapse = ", "), 
         " not found in rob.data.",
         call. = FALSE)
  }
  
  # Select required variables
  rob.data[req.variables] -> d.rob
  d.rob[,-1] %>% as.list() %>% 
    lapply(function(x){x[!x %in% c("Yes", "No", "NI")] = NA; x}) %>% 
    as.data.frame() %>% cbind(d.rob[,1],.) -> d.rob
  
  # Check if required variables are in data
  if (!"study" %in% colnames(data)) {
    stop("'study' variable not found in data.", call. = FALSE)
  }
  if (!"attr_arm1" %in% colnames(data)) {
    stop("'attr_arm1' variable not found in data.", call. = FALSE)
  }
  if (!"attr_arm2" %in% colnames(data)) {
    stop("'attr_arm2' variable not found in data.", call. = FALSE)
  }
  if (!"rand_arm1" %in% colnames(data)) {
    stop("'rand_arm1' variable not found in data.", call. = FALSE)
  }
  if (!"rand_arm2" %in% colnames(data)) {
    stop("'rand_arm2' variable not found in data.", call. = FALSE)
  }
  if (!"rating" %in% colnames(data)) {
    stop("'rating' variable not found in data.", call. = FALSE)
  }
  
  # Sort data
  data = data[order(data$study),]
  
  # Get studies not found in data
  if (nrow(data[data$study %in% d.rob$study,])==0) {
    stop("None of the studies in 'data.rob' were found", 
         " in the database. Check the 'study' variable.", call. = FALSE)
  }
  
  has.studies = unique(data[data$study %in% d.rob$study,][["study"]])
  miss.studies = unique(data[!data$study %in% d.rob$study,][["study"]])
  matrix(rep(NA, length(miss.studies) * ncol(d.rob)), ncol = ncol(d.rob)) %>% 
    {colnames(.) = colnames(d.rob); .[,"study"] = miss.studies;.} %>% 
    as.data.frame() %>% rbind(d.rob, .) -> d.rob
  
  # Expand rob data
  data %>% split(.$study) %>% lapply(nrow) -> recoder
  d.rob$.rep = recode(d.rob$study, !!!recoder)
  d.rob %>% split(.$study) %>% lapply(function(x) x[rep(1,x$.rep),]) %>% 
    do.call(rbind,.) -> d.rob
  
  # Check if join is possible; merge
  if (sum(!data$study==d.rob$study)==0) {
    message("- ", crayon::green("[OK] "), "Datasets joined.")}
  cbind(data, d.rob[,-c(1,ncol(d.rob))]) %>% {rownames(.)=NULL;.} -> d.merge
  
  # Prepare required MP-standard variables
  if (!is.numeric(d.merge$rand_arm1)) 
    as.numeric(d.merge$rand_arm1) -> d.merge$rand_arm1
  if (!is.numeric(d.merge$rand_arm2)) 
    as.numeric(d.merge$rand_arm2) -> d.merge$rand_arm2
  if (!is.numeric(d.merge$attr_arm1)) 
    as.numeric(d.merge$attr_arm1) -> d.merge$attr_arm1
  if (!is.numeric(d.merge$attr_arm2)) 
    as.numeric(d.merge$attr_arm2) -> d.merge$attr_arm2
  if (!is.numeric(d.merge$rand_arm1)) 
    as.numeric(d.merge$rand_arm1) -> d.merge$rand_arm1
  if (!is.numeric(d.merge$rand_arm2)) 
    as.numeric(d.merge$rand_arm2) -> d.merge$rand_arm2
  if (!is.numeric(d.merge$attr_arm2)) 
    as.numeric(d.merge$attr_arm2) -> d.merge$attr_arm2
  if (!is.numeric(d.merge$attr_arm2)) 
    as.numeric(d.merge$attr_arm2) -> d.merge$attr_arm2
  
  # Check if attrition is higher than randomized N
  within(data, {
    .attr_arm1 = attr_arm1; .attr_arm2 = attr_arm2
    attr_arm1 = attr_arm1/rand_arm1; attr_arm2 = attr_arm2/rand_arm2;
  }) -> data
  if (sum(data$attr_arm1>=1, na.rm=TRUE)>0)
    stop("'attr_arm1' contains sample sizes larger than 'rand_arm1'", call. = FALSE)
  if (sum(data$attr_arm2>=1, na.rm=TRUE)>0)
    stop("'attr_arm2' contains sample sizes larger than 'rand_arm2'", call. = FALSE)
  
  # Prepare ratio variables
  strsplit(d.merge$rand_ratio, "to") %>% 
    lapply(function(x) paste(x, collapse = ":")) %>% 
    stringr::str_replace_all(" ", "") %>% 
    {ifelse(grepl(":", .), ., NA)} -> d.merge$rand_ratio
  
  # Adapt D1_3 (Randomized Sample Size Balance)
  apply(d.merge, 1, function(x){
    try(testRandomizedProportion(
      as.numeric(x[["n_arm1"]]), as.numeric(x[["n_arm2"]]), 
      x[["study"]], ratio = x[["rand_ratio"]]), 
    silent = TRUE) -> tmp
    if (class(tmp)[1]=="try-error") {"NI"}
    else if (is.na(tmp$p)) {"NI"}
    else if (tmp$p < 0.01) {"No"}
    else {"Yes"}
  }) %>% as.vector() -> d.merge$d1_3
  
  message("- ", crayon::green("[OK] "), 
          "Rated item 3 (Domain 1) from database.")
  
  # Adapt D1_4 (Baseline Imbalance)
  apply(d.merge, 1, function(x){
    try(testBaselineImbalance(
      as.numeric(x[["baseline_m_arm1"]]), as.numeric(x[["baseline_m_arm2"]]), 
      as.numeric(x[["baseline_sd_arm1"]]), as.numeric(x[["baseline_sd_arm2"]]), 
      as.numeric(x[["baseline_n_arm1"]]), as.numeric(x[["baseline_n_arm2"]]), 
      study = x[["study"]]), silent = TRUE) -> tmp
    if (class(tmp)[1]=="try-error") {"NI"}
    else if (is.na(tmp$p)) {"NI"}
    else if (tmp$p < 0.01) {"No"}
    else {"Yes"}
  }) %>% as.vector() -> d.merge$d1_4
  
  message("- ", crayon::green("[OK] "), 
          "Rated item 4 (Domain 1) from database.")
  
  # Adapt D3_10 & D3_13 (missingness)
  ifelse(ifelse(d.merge$attr_arm1<.05, "Yes", "No") == "Yes" &
           ifelse(d.merge$attr_arm2<.05, "Yes", "No") == "Yes",
         "Yes", "No") %>% {.[is.na(d.merge$attr_arm1)]="NI";.} -> d.merge$d3_10
  ifelse(ifelse(d.merge$attr_arm1<=.3, "Yes", "No") == "Yes" &
           ifelse(d.merge$attr_arm2<=.3, "Yes", "No") == "Yes",
         "Yes", "No") %>% {.[is.na(d.merge$attr_arm1)]="NI";.} -> d.merge$d3_13
  ifelse(d.merge$rating=="self-report", "Yes", "No") %>% 
    {.[is.na(d.merge$rating)]="NI";.} -> d.merge$d4_16
  ifelse(d.merge$rating=="self-report", "No", "Yes") %>% 
    {.[is.na(d.merge$rating)]="NI";.} -> d.merge$d4_17
  
  # Adapt D3_12 (larger trial)
  ifelse(ifelse(d.merge$rand_arm1>=40, "Yes", "No") == "Yes" &
           ifelse(d.merge$rand_arm1>=40, "Yes", "No") == "Yes",
         "Yes", "No") %>% {.[is.na(d.merge$rand_arm1)]="NI";.} -> d.merge$d3_12
  
  message("- ", crayon::green("[OK] "), "Rated item 10 (Domain 3) from database.")
  message("- ", crayon::green("[OK] "), "Rated item 12 (Domain 3) from database.")
  message("- ", crayon::green("[OK] "), "Rated item 13 (Domain 3) from database.")
  message("- ", crayon::green("[OK] "), "Rated item 16 (Domain 4) from database.")
  message("- ", crayon::green("[OK] "), "Rated item 17 (Domain 4) from database.")
  
  # Loop through RoB rating algorithms
  apply(d.merge, 1, function(x){
    data.frame(
      d1 = domain1Algorithm(x[["d1_1"]], x[["d1_2"]], x[["d1_3"]], x[["d1_4"]]),
      d2 = domain2Algorithm(x[["d2_5"]], x[["d2_6"]], x[["d2_7"]], x[["d2_8"]], x[["d2_9"]]),
      d3 = domain3Algorithm(x[["d3_10"]], x[["d3_11"]], x[["d3_12"]], x[["d3_13"]], x[["d3_14"]]),
      d4 = domain4Algorithm(x[["d4_15"]], x[["d4_16"]], x[["d4_17"]], x[["d4_18"]]),
      d5 = domain5Algorithm(x[["d5_19"]], x[["d5_20"]], x[["d5_21"]], x[["d5_22"]], x[["d5_23"]], x[["d5_24"]])) %>% 
      {.$.rob = overallAlgorithm(.$d1, .$d2, .$d3, .$d4, .$d5);.}
  }) %>% do.call(rbind, .) %>% cbind(d.merge, .) %>% 
    within({rob = .rob; .rob = NULL}) -> d.merge
  
  # Re-attach old 'attr' variables
  d.merge$attr_arm1 = data$.attr_arm1; data$.attr_arm1 = NULL
  d.merge$attr_arm2 = data$.attr_arm2; data$.attr_arm2 = NULL
  
  # Create rob_study_lvl variable
  d.merge %>% split(.$study) %>% lapply(function(x){
    prio = c("high risk" = 3, "some concerns" = 2, "low risk" = 1)
    x$rob -> r.rating
    suppressWarnings(names(
      which(prio==max(prio[unique(r.rating)], na.rm = T)))) -> r.new
    if (length(r.new)>0)
      x$rob_study_lvl = rep(r.new, nrow(x))
    else 
      x$rob_study_lvl = rep(NA, nrow(x))
    x
  }) %>% do.call(rbind, .) %>% 
    {rownames(.) = NULL;.} -> d.merge
  
  d.merge[colnames(data)] %>% 
    {.$rob = NULL;.} %>% 
    {cbind(., d.merge[c("d1", "d2", "d3", "d4", "d5", 
                        "rob", "rob_study_lvl")] %>% 
             rename("rob_d1" = "d1", "rob_d2" = "d2",
                    "rob_d3" = "d3", "rob_d3" = "d3",
                    "rob_d4" = "d4", "rob_d5" = "d5"))} -> database
  
  d.merge[
    d.merge$study %in% has.studies,
    c("study", "d1_1", "d1_2", "d1_3", "d1_4", "d1_notes", "d2_5", "d2_6", "d2_7", "d2_8", 
      "d2_9", "d2_notes", "d3_10", "d3_11", "d3_12", "d3_13", "d3_14", "d3_notes", 
      "d4_15", "d4_16", "d4_17", "d4_18", "d4_notes", "d5_19", "d5_20", "d5_21", 
      "d5_22", "d5_23", "d5_24", "d5_notes", "attr_arm1", "d1", "d2", "d3", 
      "d4", "d5", "rob", "rob_study_lvl")] -> rob.data
  try({d.merge[d.merge$study %in% has.studies, ".id"] -> rob.data$.id}, 
      silent = TRUE)
  
  message("- ", crayon::green("[OK] "), "Done!")
  
  ret.obj = list(database = database,
                 rob.data = rob.data,
                 miss.studies = miss.studies)
  
  class(ret.obj) = c("createRobRatings", "list")
  return(ret.obj)

}


#' Print method for objects of class 'createRobRatings'
#'
#' Print S3 method for objects of class \code{createRobRatings}.
#'
#' @param x An object of class \code{createRobRatings}.
#' @param ... Additional arguments.
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @importFrom crayon green blue magenta bold
#'
#' @export
#' @method print createRobRatings
#' @keywords internal 

print.createRobRatings = function(x, ...){
  cat(crayon::bold("'createRobRatings' object \n"))
  cat("- ", crayon::bold(crayon::green("database")), ": Metapsy database (", 
      nrow(x$database), " rows, ", length(unique(x$database[["study"]])), 
      " studies) \n", sep = "")
  cat("- ", crayon::bold(crayon::green("rob.data")), ": RoB sheet (", 
      nrow(x$rob.data), " rows, ", length(unique(x$rob.data[["study"]])), 
      " studies)\n", sep = "")
  cat("- ", crayon::bold(crayon::green("miss.studies")), 
      ": studies not found in extraction sheet (", 
      length(x$miss.studies), ")\n", sep = "")
}




