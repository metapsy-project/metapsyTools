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
#' library(xlsx)
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
#' 
#' # Save both files into a new MS Excel file
#' xlsx::write.xlsx(tmp$database, "data_rated.xlsx", 
#'                  sheetName = "database", showNA = FALSE)
#' xlsx::write.xlsx(tmp$rob.data, "data_rated.xlsx", 
#'                  sheetName = "rob", showNA = FALSE,
#'                  append = TRUE)
#' }
#' 
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Clara Miguel Sanz \email{clara.miguelsanz@@vu.nl}, 
#' Pim Cuijpers \email{p.cuijpers@@vu.nl}
#' 
#' @seealso \code{\link{checkRobDiscrepancies}}
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
  if (!"rand_ratio" %in% colnames(data)) {
    stop("'rand_ratio' variable not found in data.",
         " This variable should provide the allocation ratio ('1 to 1', ",
         "'2 to 1', ...) for each effect.", 
         call. = FALSE)
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
  
  # Check for RoB studies not in database
  miss.studies.rob = unique(d.rob[!d.rob$study %in% data$study,][["study"]])
  d.rob = d.rob[!d.rob$study %in% miss.studies.rob,]
  
  # Expand rob data
  data %>% split(.$study) %>% lapply(nrow) -> recoder
  d.rob$.rep = dplyr::recode(d.rob$study, !!!recoder)
  d.rob %>% split(.$study) %>% lapply(function(x) x[rep(1,x$.rep),]) %>% 
    do.call(rbind,.) -> d.rob
  
  # Check if join is possible; merge
  if (sum(!data$study==d.rob$study)==0) {
    message("- ", crayon::green("[OK] "), "Datasets joined.")} else {
      stop("Datasets could not be merged.", call. = FALSE)
    }; cbind(data, d.rob[,-c(1,ncol(d.rob))]) %>% {rownames(.)=NULL;.} -> d.merge
  
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
  if (sum(data$attr_arm1>=1, na.rm=TRUE)>0) {
    stdy = data[data$attr_arm1>=1 & !is.na(data$attr_arm1),][["study"]] %>% 
      unique() %>% paste(collapse = "; ")
    stop("'attr_arm1' contains sample sizes larger than 'rand_arm1' in ", 
         stdy, ".", call. = FALSE)
  }
  if (sum(data$attr_arm2>=1, na.rm=TRUE)>0) {
    stdy = data[data$attr_arm2>=1 & !is.na(data$attr_arm2),][["study"]] %>% 
      unique() %>% paste(collapse = "; ")
    stop("'attr_arm2' contains sample sizes larger than 'rand_arm2' in ", 
         stdy, ".", call. = FALSE)
  }
  
  # Prepare ratio variables
  strsplit(d.merge$rand_ratio, "to") %>% 
    lapply(function(x) paste(x, collapse = ":")) %>% 
    stringr::str_replace_all(" ", "") %>% 
    {ifelse(grepl(":", .), ., NA)} -> d.merge$rand_ratio
  
  # Adapt D1_3 (Randomized Sample Size Balance)
  apply(d.merge, 1, function(x){
    try(testRandomizedProportion(
      as.numeric(x[["rand_arm1"]]), as.numeric(x[["rand_arm2"]]), 
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
  d.merge$attr_arm1/d.merge$rand_arm1 -> p.attr.arm1
  d.merge$attr_arm2/d.merge$rand_arm2 -> p.attr.arm2
  ifelse(ifelse(p.attr.arm1<.05, "Yes", "No") == "Yes" &
           ifelse(p.attr.arm2<.05, "Yes", "No") == "Yes",
         "Yes", "No") %>% {.[is.na(d.merge$attr_arm1)]="NI";.} %>% 
                          {.[is.na(d.merge$attr_arm2)]="NI";.} %>%
                          {.[is.na(d.merge$rand_arm1)]="NI";.} %>% 
                          {.[is.na(d.merge$attr_arm2)]="NI";.}-> d.merge$d3_10
  ifelse(ifelse(p.attr.arm1<=.3, "Yes", "No") == "Yes" &
           ifelse(p.attr.arm2<=.3, "Yes", "No") == "Yes",
         "Yes", "No") %>% {.[is.na(d.merge$attr_arm1)]="NI";.} %>% 
                          {.[is.na(d.merge$attr_arm2)]="NI";.} %>% 
                          {.[is.na(d.merge$rand_arm1)]="NI";.} %>% 
                          {.[is.na(d.merge$attr_arm2)]="NI";.}-> d.merge$d3_13
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
  
  # Create comparisons variable
  d.merge %>% split(.$study) %>% 
    lapply(function(x) {
     paste0(x$condition_arm1, 
           ifelse(!is.na(x$multi_arm1), paste0(" [", x$multi_arm1, "]"), "")) %>% 
        unique() %>% paste(collapse = "; ") -> tmp1
      paste0(x$condition_arm2, 
             ifelse(!is.na(x$multi_arm2), paste0(" [", x$multi_arm2, "]"), "")) %>% 
        unique() %>% paste(collapse = "; ") -> tmp2
      c(tmp1, tmp2) %>% paste(collapse = " | ") %>% 
        rep(nrow(x))
    }) %>% unlist() -> d.merge$.comparisons
  
  # Create instruments variable
  d.merge %>% split(.$study) %>% 
    lapply(function(x) {
        rep(unique(x$instrument) %>% paste(collapse = "; "), nrow(x))
    }) %>% unlist() -> d.merge$.instruments
  
  d.merge[colnames(data)] %>% 
    {.$rob = NULL;.} %>% 
    {cbind(., d.merge[c("d1", "d2", "d3", "d4", "d5", 
                        "rob", "rob_study_lvl")] %>% 
             rename("rob_d1" = "d1", "rob_d2" = "d2",
                    "rob_d3" = "d3", "rob_d3" = "d3",
                    "rob_d4" = "d4", "rob_d5" = "d5"))} -> database
  
  d.merge[
    d.merge$study %in% has.studies,
    c("study", ".comparisons", ".instruments", "d1_1", "d1_2", "d1_3", "d1_4", 
      "d1_notes", "d2_5", "d2_6", "d2_7", "d2_8", "d2_9", "d2_notes", "d3_10", 
      "d3_11", "d3_12", "d3_13", "d3_14", "d3_notes", "d4_15", "d4_16", "d4_17", 
      "d4_18", "d4_notes", "d5_19", "d5_20", "d5_21", "d5_22", "d5_23", "d5_24", 
      "d5_notes", "attr_arm1", "attr_arm2", "rand_arm1","rand_arm2", "d1", 
      "d2", "d3", "d4", "d5", "rob", "rob_study_lvl")] -> rob.data
  try({d.merge[d.merge$study %in% has.studies, ".id"] -> rob.data$.id}, 
      silent = TRUE)
  
  # Collapse rob.data per study
  rob.data[order(rob.data$study),] -> rob.data
  rob.data %>% split(.$study) %>% 
    {data.frame(
      d1_1 = lapply(.,collapseRatings, var="d1_1") %>% do.call(c,.),
      d1_2 = lapply(.,collapseRatings, var="d1_2") %>% do.call(c,.),
      d1_3 = lapply(.,collapseRatings, var="d1_3") %>% do.call(c,.),
      d1_4 = lapply(.,collapseRatings, var="d1_4") %>% do.call(c,.),
      d2_5 = lapply(.,collapseRatings, var="d2_5") %>% do.call(c,.),
      d2_6 = lapply(.,collapseRatings, var="d2_6") %>% do.call(c,.),
      d2_7 = lapply(.,collapseRatings, var="d2_7") %>% do.call(c,.),
      d2_8 = lapply(.,collapseRatings, var="d2_8") %>% do.call(c,.),
      d2_9 = lapply(.,collapseRatings, var="d2_9") %>% do.call(c,.),
      d3_10 = lapply(.,collapseRatings, var="d3_10") %>% do.call(c,.),
      d3_11 = lapply(.,collapseRatings, var="d3_11") %>% do.call(c,.),
      d3_12 = lapply(.,collapseRatings, var="d3_12") %>% do.call(c,.),
      d3_13 = lapply(.,collapseRatings, var="d3_13") %>% do.call(c,.),
      d3_14 = lapply(.,collapseRatings, var="d3_14") %>% do.call(c,.),
      d4_15 = lapply(.,collapseRatings, var="d4_15") %>% do.call(c,.),
      d4_16 = lapply(.,collapseRatings, var="d4_16") %>% do.call(c,.),
      d4_17 = lapply(.,collapseRatings, var="d4_17") %>% do.call(c,.),
      d4_18 = lapply(.,collapseRatings, var="d4_18") %>% do.call(c,.),
      d5_19 = lapply(.,collapseRatings, var="d5_19") %>% do.call(c,.),
      d5_20 = lapply(.,collapseRatings, var="d5_20") %>% do.call(c,.),
      d5_21 = lapply(.,collapseRatings, var="d5_21") %>% do.call(c,.),
      d5_22 = lapply(.,collapseRatings, var="d5_22") %>% do.call(c,.),
      d5_23 = lapply(.,collapseRatings, var="d5_23") %>% do.call(c,.),
      d5_24 = lapply(.,collapseRatings, var="d5_24") %>% do.call(c,.),
      d1 = lapply(.,collapseDomainRatings, var="d1") %>% do.call(c,.),
      d2 = lapply(.,collapseDomainRatings, var="d2") %>% do.call(c,.),
      d3 = lapply(.,collapseDomainRatings, var="d3") %>% do.call(c,.),
      d4 = lapply(.,collapseDomainRatings, var="d4") %>% do.call(c,.),
      d5 = lapply(.,collapseDomainRatings, var="d5") %>% do.call(c,.),
      attr_arm1 = lapply(., function(x) 
        suppressWarnings(max(x$attr_arm1/x$rand_arm1, na.rm = TRUE))) %>% 
        do.call(c,.) %>% {.[.==-Inf]=NA;round(.,3)*100},
      attr_arm2 = lapply(., function(x) 
        suppressWarnings(max(x$attr_arm2/x$rand_arm2, na.rm = TRUE))) %>% 
        do.call(c,.) %>% {.[.==-Inf]=NA;round(.,3)*100})} %>% 
    {cbind(rob.data[!colnames(rob.data) %in% colnames(.)] %>% 
             dplyr::distinct(study, .keep_all = TRUE),.)} %>% 
    {.[colnames(rob.data)]} %>% dplyr::select(-rand_arm1, -rand_arm2, -rob) %>% 
    rename(instruments = .instruments, comparisons = .comparisons) %>% 
    {rownames(.)=NULL;.} -> rob.data
  
  message("- ", crayon::green("[OK] "), "Done!")
  
  ret.obj = list(database = database,
                 rob.data = rob.data,
                 miss.studies = miss.studies,
                 miss.studies.rob = miss.studies.rob)
  
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
      ": database studies not found in extraction sheet (", 
      length(x$miss.studies), ")\n", sep = "")
  cat("- ", crayon::bold(crayon::green("miss.studies.rob")), 
      ": extraction sheet studies not found in database (", 
      length(x$miss.studies.rob), ")\n", sep = "")
}


#' Check for inconsistencies between two RoB extraction sheets
#'
#' Based on prepared extraction sheets by two independent raters, this function allows you 
#' to automatically check for differing studies and/or inconsistencies in the
#' ratings. Cohen's kappa is calculated automatically, differentiating between
#' "Yes" vs. "No"/"No Information" ratings to obtain the base rate of agreement.
#'
#' @usage checkRobDiscrepancies(data.1, data.2)
#'
#' @param data.1 An RoB extraction sheet. Columns of this file must have the same name as in the
#' [RoB extraction sheet template](https://www.metapsy.org/assets/files/rob-template.xlsx) 
#' provided by the Metapsy initative. If the Metapsy template has been
#' used, make sure to delete the top rows before importing, so that only the metapsyTools variables
#' remain as the column names. Required columns are: `study`, `d1_1`, `d1_2`, `d1_3`, `d1_4`, `d1_notes`, 
#' `d2_5`, `d2_6`, `d2_7`, `d2_8`, `d2_9`, `d2_notes`, `d3_10`, `d3_11`, `d3_12`, `d3_13`, 
#' `d3_14`, `d3_notes`, `d4_15`, `d4_16`, `d4_17`, `d4_18`, `d4_notes`, `d5_19`, 
#' `d5_20`, `d5_21`, `d5_22`, `d5_23`, `d5_24`, `d5_notes`. 
#' @param data.2 The same extraction sheet from another (i.e., second) rater.
#'
#' @return If discrepancies are found, the function will return a data frame with
#' the respective study/studies, along with the diverging ratings (`discrepancies` element). 
#' If different studies are included in both sheets, they will be saved under `diff.studies`.
#'
#' @examples
#' \dontrun{
#' library(readxl)
#' 
#' # Get example extraction sheet from metapsy.org/assets/files/rob_data.xlsx
#' rob_data <- read_excel("rob_data.xlsx")
#' 
#' # Create second sheet with partly different ratings
#' rob_data -> rob_data_2
#' rob_data_2[-1,] -> rob_data_2
#' rob_data_2[1,"d2_5"] = "NI"
#' rob_data_2[1,"d4_15"] = "No"
#' 
#' # Check for discrepancies
#' tmp <- metapsyTools:::checkRobDiscrepancies(rob_data, rob_data_2)
#' tmp
#' }
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Clara Miguel Sanz \email{clara.miguelsanz@@vu.nl}, 
#' Pim Cuijpers \email{p.cuijpers@@vu.nl}
#' 
#' @seealso \code{\link{createRobRatings}}
#'
#' @import dplyr
#' @importFrom crayon green yellow cyan bold
#' @keywords internal

checkRobDiscrepancies = function(data.1, data.2) {
  
  # Required variables
  req.variables = c("study", "d1_1", "d1_2", 
                    "d1_3", "d1_4", "d1_notes", 
                    "d2_5", "d2_6", "d2_7", 
                    "d2_8", "d2_9", "d2_notes", 
                    "d3_10", "d3_11", "d3_12", 
                    "d3_13", "d3_14", "d3_notes", 
                    "d4_15", "d4_16", "d4_17", 
                    "d4_18", "d4_notes", "d5_19", 
                    "d5_20", "d5_21", "d5_22", 
                    "d5_23", "d5_24", 
                    "d5_notes")
  
  # Check if required variables are included
  if (sum(!req.variables %in% colnames(data.1)) > 0) {
    stop("Required variables ", 
         paste(req.variables[!req.variables %in% colnames(data.1)], 
               collapse = ", "), 
         " not found in 'data.1'.",
         call. = FALSE)
  }
  if (sum(!req.variables %in% colnames(data.2)) > 0) {
    stop("Required variables ", 
         paste(req.variables[!req.variables %in% colnames(data.2)], 
               collapse = ", "), 
         " not found in 'data.2'.",
         call. = FALSE)
  }
  
  # Select required variables
  data.1[req.variables] -> data.1
  data.1[,-1] %>% as.list() %>% 
    lapply(function(x){x[!x %in% c("Yes", "No", "NI")] = NA; x}) %>% 
    as.data.frame() %>% cbind(data.1[,1],.) -> data.1
  data.2[req.variables] -> data.2
  data.2[,-1] %>% as.list() %>% 
    lapply(function(x){x[!x %in% c("Yes", "No", "NI")] = NA; x}) %>% 
    as.data.frame() %>% cbind(data.2[,1],.) -> data.2
  
  # Check if required variables are in data
  if (!"study" %in% colnames(data.1)) {
    stop("'study' variable not found in data.1.", call. = FALSE)
  }
  if (!"study" %in% colnames(data.2)) {
    stop("'study' variable not found in data.2.", call. = FALSE)
  }
  
  # Find inconsistent studies
  setdiff(union(data.1$study,data.2$study), 
          intersect(data.1$study,data.2$study)) -> diff.studies
  if (length(diff.studies)>0) {
    message("- ", crayon::yellow("[!] "), 
            "Studies differ between datasets (",
            length(diff.studies), ").") } else {
              message("- ", crayon::green("[OK] "), 
                      "Studies match between datasets.")}
  
  # Get common studies
  data.1[data.1$study %in% intersect(data.1$study,data.2$study),] -> data.1
  data.2[data.2$study %in% intersect(data.1$study,data.2$study),] -> data.2
  if (nrow(data.1)==0|nrow(data.2)==0) {
    stop("No overlapping studies found in both datasets. Check the 'study' variables.", 
         call. = FALSE)
  }
  data.1[order(data.1$study),] -> data.1
  data.2[order(data.2$study),] -> data.2
  list(data.1, data.2) %>% Reduce(`==`,.) %>% {!.} -> mat
  l = list()
  for (i in 1:nrow(mat)) {
    data.frame(
      data.1 = (data.1[i,mat[i,] & !is.na(mat[i,])]) %>% 
        paste(names(.), ., collapse = "; "),
      data.2 = (data.2[i,mat[i,] & !is.na(mat[i,])]) %>% 
        paste(names(.), ., collapse = "; ")) -> l[[i]]
  }
  cbind(study = data.1$study[rowSums(mat, na.rm=T)>0],
        do.call(rbind, l)[rowSums(mat, na.rm=T)>0,]) -> discrepancies
  
  # Calculate Cohen's kappa
  sum(!is.na(mat[,colnames(mat)!="study"]), na.rm = TRUE) -> n
  (n-sum(mat[,colnames(mat)!="study"], na.rm = TRUE))/n -> p0
  c((sum(data.1[,colnames(data.1)!="study"] == "Yes", na.rm = T)/n),
    (sum(data.2[,colnames(data.2)!="study"] == "Yes", na.rm = T)/n)) %>% 
    {prod(.)+prod(1-.)} -> pe
  kappa = (p0-pe)/(1-pe)
  
  if (nrow(discrepancies)==0) {
    message("- ", crayon::green("[OK] "), 
            "No discrepancies detected!")
    message("- ", crayon::green("[OK] "), 
            "Inter-rater agreement (Cohen's kappa): ", round(kappa,2))
    if (length(diff.studies)>0) {
      return(list(diff.studies = diff.studies,
                  kappa = kappa))
    }
  } else {
    message("- ", crayon::yellow("[!] "), 
            "Discrepancies detected (", nrow(discrepancies), ").")
    message("- ", crayon::green("[OK] "), 
            "Inter-rater agreement (Cohen's kappa): ", round(kappa,2))
    if (length(diff.studies)>0) {
      return(list(discrepancies = discrepancies,
                  diff.studies = diff.studies,
                  kappa = kappa))
    } else {
      return(list(discrepancies = discrepancies,
                  kappa = kappa)) 
    }
  }
}



