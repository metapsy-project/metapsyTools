#' Convert Metapsy database into network meta-analysis format
#'
#' This function converts a database following the 
#' [Metapsy data standard](https://docs.metapsy.org/data-preparation/format/)
#' into format that is suitable for network meta-analysis (e.g. using [netmeta::netmeta()]). 
#' 
#'
#' @param mean_arm1 Mean score in the first trial arm.
#' @param mean_arm2 Mean score in the second trial arm.
#' @param sd_arm1 Standard deviation in the first trial arm.
#' @param sd_arm2 Standard deviation in the second trial arm.
#' @param n_arm1 Sample size in the first trial arm.
#' @param n_arm2 Sample size in the second trial arm.
#' @param mean_change_arm1 Mean change scores in the first trial arm.
#' @param mean_change_arm2 Mean change scores in the second trial arm.
#' @param sd_change_arm1 Standard deviation of change scores in the first trial arm.
#' @param sd_change_arm2 Standard deviation of change scores in the second trial arm.
#' @param n_change_arm1 Sample size of change scores in the first trial arm. 
#' @param n_change_arm2 Sample size of change scores in the second trial arm.
#' @param event_arm1 Number of responders in the first trial arm.
#' @param event_arm2 Number of responders in the second trial arm.
#' @param totaln_arm1 Total number of participants in the first trial arm.
#' @param totaln_arm2 Total number of participants in the second trial arm.
#' @param condition_arm1 Treatment or format in the first trial arm.
#' @param condition_arm2 Treatment or format in the second trial arm.
#' @param study Study labels for each comparison.
#' @param ... Additional arguments. Can be used to specify additional columns to be included in the output (see Details).
#' @param data Dataset following the [Metapsy data standard](https://docs.metapsy.org/data-preparation/format/) (optional).
#'
#' @usage convertNetmeta(
#'                
#'                # Continuous outcomes (endpoint scores)
#'                mean_arm1, mean_arm2, sd_arm1, sd_arm2, 
#'                n_arm1, n_arm2,
#'                
#'                # Continuous outcomes (change scores)
#'                mean_change_arm1, mean_change_arm2, 
#'                sd_change_arm1, sd_change_arm2,
#'                n_change_arm1, n_change_arm2, 
#'                
#'                # Response (event counts)
#'                event_arm1, event_arm2, totaln_arm1, 
#'                totaln_arm2, 
#'                
#'                # Study characteristics
#'                condition_arm1, condition_arm2, study, 
#'                ..., data = NULL)
#'                
#'                
#' @details
#' This function converts a Metapsy database into a "wider" format dataset that can be used to
#' run network meta-analyses. Returned objects are optimized for [netmeta::netmeta()] and can 
#' be used "out-of-the box" in this package.
#' 
#' The function will perform an expansion of multi-arm trials, which is required
#' for most network meta-analysis implementations. Thus, the function will calculate all
#' three unique comparisons for three-arm trials, all six comparisons for four-arm trials, 
#' etc.  
#' 
#' Two additional formatting requirements must be met to conduct the conversion:
#' - Each comparison in a trial is only allowed to provide exactly one effect size/contrast.
#' This may be resolved by filtering the dataset beforehand using [filterPoolingData] or
#' [filterPriorityRule]. The function will return an informative error message if
#' non-unique comparisons are found.
#' - The function can only use raw continuous outcome and binary response data to calculate SMDs for 
#' each comparison. Rows with other effect size information (e.g. pre-calculated effects based on
#' _t_ or _F_-tests) will not be included. If comparisons had to be removed, affected studies will
#' be printed into the console and should be added manually. 
#' 
#' It is also possible to add additional columns (e.g., columns included in the dataset provided in `data`) to
#' the final data frame. These columns have to be specified as additional arguments
#' in the function call, where the argument name will be used as the column name (see Examples).
#' 
#'
#' @return 
#' Returns a `data.frame` in wide-format, containing the calculated effect sizes
#' (standardized mean differences, SMDs) for each required comparison. The following
#' columns will be included in all outputs:
#' 
#' - `studlab`: Study label for each comparison. If response counts were used as outcome,
#' "`(response)`" is appended to the study name.
#' - `treat1`: Condition or format used in the first trial arm.
#' - `treat2`: Condition or format used in the second trial arm.
#' - `TE`: The calculated effect size (SMD).
#' - `seTE`: Standard error of the calculated effect size.
#' 
#' Depending on whether continuous or binary outcomes (or both) were used, 
#' the dataset will also include further columns containing the raw data used
#' to obtain the effect size:
#' 
#' - `n1`: Sample size in the first trial arm.
#' - `n2`: Sample size in the second trial arm.
#' - `mean1`: Mean (change) scores in the first trial arm.
#' - `mean2`: Mean (change) scores in the second trial arm.
#' - `sd1`: Standard deviation in the first trial arm.
#' - `sd2`: Standard deviation in the second trial arm.
#' - `event1`: Responders in the first trial arm.
#' - `event2`: Responders in the second trial arm.
#' 
#' Studies for which effect sizes could not be calculated will be saved as
#' a character vector in the `removed.studies` attribute. They can be extracted
#' using `attr(res, "removed.studies")`, where `res` is the returned data frame.
#' 
#' @examples
#' \dontrun{
#' 
#' # Filter database so that only unique comparisons remain for each study.
#' data <- depressionPsyCtr %>% 
#'   filterPriorityRule(
#'     instrument = c("ces-d", "phq-9", "scl", 
#'                    "hdrs", "bdi-2", "scid")) %>% 
#'   filterPoolingData(year >= 1985, study != "Barrett, 2001") 
#' 
#' # Convert endpoint, change score, and response data.
#' dat.netmeta <- convertNetmeta(
#'   
#'   # Continuous outcome data
#'   mean_arm1, mean_arm2, sd_arm1, sd_arm2, n_arm1, n_arm2,
#'   mean_change_arm1, mean_change_arm2, sd_change_arm1, 
#'   sd_change_arm2, n_change_arm1, n_change_arm2, 
#'   
#'   # Response data
#'   event_arm1 = event_arm1, event_arm2 = event_arm2, 
#'   totaln_arm1 = totaln_arm1, totaln_arm2 = totaln_arm2, 
#'   
#'   # Treatments to be used in NMA
#'   condition_arm1 = condition_arm1, 
#'   condition_arm2 = condition_arm2, 
#'   
#'   # Additional column to be added
#'   scale = instrument,
#'   
#'   # Study label and data
#'   study = study, data = data)
#' 
#' # Load netmeta and perform NMA
#' library(netmeta)
#' netmeta(TE, seTE, treat1, treat2, studlab, 
#'         data = dat.netmeta, reference.group = "wl")
#'         
#'         
#' # Example using metapsyData database
#' library(metapsyData)
#' d <- getData("depression-psyctr", version = "22.0.2")
#' 
#' d$data %>% 
#'   filterPriorityRule(instrument = c("phq-9", "ces-d", "hdrs")) %>% 
#'   filterPoolingData(!study %in%
#'                       c('Baumgartner, 2021', 'Brown, 1984', 'Fann, 2015', 
#'                         'Fledderus, 2012', 'Floyd, 2004', 'Kleiboer, 2015', 
#'                         'Lemma, 2013', 'Mohr, 2013', 'Nezu, 1989', 
#'                         'NystrÃ¶m, 2017', 'Pecheur, 1984', 'Propst, 1992', 
#'                         'Rehm, 1981', 'Rohan, 2007', 'Scogin, 1989', 'Selmi, 1990', 
#'                         'Smith, 2017a', 'Titov, 2010', 'Tomasino, 2017', 
#'                         'Watt, 2000', 'Westerhof, 2019', 'Araya, 2021', 
#'                         'Choi, 2014')) %>% 
#'   convertNetmeta(mean_arm1, mean_arm2, sd_arm1, sd_arm2, n_arm1, n_arm2,
#'                  event_arm1 = event_arm1, event_arm2 = event_arm2, 
#'                  totaln_arm1 = totaln_arm1, totaln_arm2 = totaln_arm2,
#'                  condition_arm1 = condition_arm1, condition_arm2 = condition_arm2,
#'                  study = study, format = format, data = .) -> dat.netmeta
#' 
#' # Extract studies for which no effect sizes could be calculated
#' attr(dat.netmeta, "removed.studies")
#' 
#' # Run network meta-analysis
#' netmeta(TE, seTE, treat1, treat2, studlab, data = dat.netmeta)
#' 
#' 
#' # Multi-arm expansion with single trials:
#' # - Using continuous outcome
#' convertNetmeta(mean_arm1 = c(4.12, 5.74),
#'                mean_arm2 = c(5.74, 6.41),
#'                sd_arm1 = c(4.22, 5.15),
#'                sd_arm2 = c(5.15, 2.79),
#'                n_arm1 = c(50, 50),
#'                n_arm2 = c(50, 50),
#'                condition_arm1 = c("cbt", "dyn"),
#'                condition_arm2 = c("dyn", "wl"),
#'                study = c("Doe, 1999", "Doe, 1999"))
#' 
#' # - using response outcome
#' convertNetmeta(event_arm1 = c(22, 12),
#'                event_arm2 = c(12, 5),
#'                totaln_arm1 = c(87, 89),
#'                totaln_arm2 = c(89, 92),
#'                condition_arm1 = c("cbt", "dyn"),
#'                condition_arm2 = c("dyn", "wl"),
#'                study = c("Doe, 1999", "Doe, 1999"))
#'                
#' # - using study format instead of conditions
#' format.data <- data.frame(study = c("Doe, 1999", "Doe, 1999", "Miller, 2000",
#'                                     "Willms, 2017", "Willms, 2017"),
#'                           format_arm1 = c("gsh", "ush", "ush", "gsh", "ush"),
#'                           format_arm2 = c("ush", "wl",  "cau", "ush", "cau"),
#'                           mean_arm1 = c(4.12, 5.74, 3.21, 4.99, 6.23),
#'                           mean_arm2 = c(5.74, 6.41, 6.29, 6.23, 6.41),
#'                           sd_arm1 = c(4.22, 5.15, 4.21, 4.00, 5.92),
#'                           sd_arm2 = c(5.15, 2.79, 4.52, 5.92, 3.12),
#'                           n_arm1 = c(50, 50, 76, 30, 30),
#'                           n_arm2 = c(50, 50, 75, 30, 30))
#' 
#' convertNetmeta(mean_arm1, mean_arm2, sd_arm1, sd_arm2, n_arm1, n_arm2,
#'                condition_arm1 = format_arm1, condition_arm2 = format_arm2,
#'                study = study, data = format.data)
#' }
#' 
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com}, Paula Kuper \email{paula.r.kuper@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso [netmeta::netmeta()], [checkDataFormat()]
#' @importFrom dplyr select rename bind_rows any_of arrange distinct
#' @importFrom metafor to.long
#' @importFrom stringr str_remove_all
#' @importFrom crayon bold cyan
#' @export convertNetmeta
#' @keywords internal


convertNetmeta = function(mean_arm1, mean_arm2,
                          sd_arm1, sd_arm2,
                          n_arm1, n_arm2,
                          mean_change_arm1, mean_change_arm2,
                          sd_change_arm1, sd_change_arm2,
                          n_change_arm1, n_change_arm2,
                          event_arm1, event_arm2,
                          totaln_arm1, totaln_arm2,
                          condition_arm1, condition_arm2,
                          study, ..., data = NULL) {
  
  
  # Check if netmeta is available
  if (!requireNamespace("netmeta", quietly = TRUE))
    stop("Package \"netmeta\" must be installed to use this function.", call. = FALSE)

  # get data; extract from parent frame if necessary
  isNullData = is.null(data)
  parentFrame = sys.frame(sys.parent())
  mC = match.call()
  mCd = match.call(expand.dots = FALSE)$...
  if (isNullData) { data = parentFrame}
  
  # get data & check missings
  mean_arm1 = try({catch("mean_arm1", mC, data, parentFrame)}, silent = TRUE)
  mean_arm2 = try({catch("mean_arm2", mC, data, parentFrame)}, silent = TRUE)
  sd_arm1 = try({catch("sd_arm1", mC, data, parentFrame)}, silent = TRUE)
  sd_arm2 = try({catch("sd_arm2", mC, data, parentFrame)}, silent = TRUE)
  n_arm1 = try({catch("n_arm1", mC, data, parentFrame)}, silent = TRUE)
  n_arm2 = try({catch("n_arm2", mC, data, parentFrame)}, silent = TRUE)
  mean_change_arm1 = try({catch("mean_change_arm1", mC, data, parentFrame)}, silent = TRUE)
  mean_change_arm2 = try({catch("mean_change_arm2", mC, data, parentFrame)}, silent = TRUE)
  sd_change_arm1 = try({catch("sd_change_arm1", mC, data, parentFrame)}, silent = TRUE)
  sd_change_arm2 = try({catch("sd_change_arm2", mC, data, parentFrame)}, silent = TRUE)
  n_change_arm1 = try({catch("n_change_arm1", mC, data, parentFrame)}, silent = TRUE)
  n_change_arm2 = try({catch("n_change_arm2", mC, data, parentFrame)}, silent = TRUE)
  condition_arm1 = try({catch("condition_arm1", mC, data, parentFrame)}, silent = TRUE)
  condition_arm2 = try({catch("condition_arm2", mC, data, parentFrame)}, silent = TRUE)
  event_arm1 = try({catch("event_arm1", mC, data, parentFrame)}, silent = TRUE)
  event_arm2 = try({catch("event_arm2", mC, data, parentFrame)}, silent = TRUE)
  totaln_arm1 = try({catch("totaln_arm1", mC, data, parentFrame)}, silent = TRUE)
  totaln_arm2 = try({catch("totaln_arm2", mC, data, parentFrame)}, silent = TRUE)
  study = try({catch("study", mC, data, parentFrame)}, silent = TRUE)
  
  # Throw error if essential elements are missing
  nulls = c(is.null(condition_arm1), is.null(condition_arm2), is.null(study))
  if (sum(nulls)>0) {
    c("condition_arm1", "condition_arm2", "study")[nulls] -> missings
    stop("Argument(s) ", paste(missings, sep = ", "), 
         " missing, with no default.", call. = FALSE)
  }
  
  # Get additional variables
  l = list()
  extra.vars = sapply(mCd[sapply(mCd, is.name)], deparse)
  for (i in extra.vars) {
    l[[i]] = try({eval(as.symbol(i), data, enclos = parentFrame)}, 
                 silent = TRUE)
  }; do.call(cbind, l[lapply(l, class) != "try-error"]) -> extra.dat
  colnames(extra.dat)[names(mCd) != ""] = names(mCd)[names(mCd) != ""]
  extra.vars = colnames(extra.dat)
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Continuous outcome module                         #
  # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  if (!any(is.null(mean_arm1), is.null(mean_arm2), 
           is.null(sd_arm1), is.null(sd_arm2),
           is.null(n_arm1), is.null(n_arm2))) {
    
    # Where post-test is missing, replace with change
    if (!any(is.null(mean_change_arm1), is.null(mean_change_arm2), 
             is.null(sd_change_arm1), is.null(sd_change_arm2),
             is.null(n_change_arm1), is.null(n_change_arm2))) {
      mask = is.na(mean_arm1)
      mean_arm1[mask] = mean_change_arm1[mask]
      mean_arm2[mask] = mean_change_arm2[mask]
      sd_arm1[mask] = sd_change_arm1[mask]
      sd_arm2[mask] = sd_change_arm2[mask]
      n_arm1[mask] = n_change_arm1[mask]
      n_arm2[mask] = n_change_arm2[mask]
    }
    
    # Combine all variables in data frame
    dat = cbind(mean_arm1, mean_arm2, sd_arm1, sd_arm2,
                n_arm1, n_arm2, condition_arm1, condition_arm2, 
                study, extra.dat)
    
    # Pivot data frame to long format
    suppressWarnings({metafor::to.long("SMD", 
                                       m1i = as.numeric(mean_arm1), sd1i = as.numeric(sd_arm1), 
                                       n1i = as.numeric(n_arm1),
                                       m2i = as.numeric(mean_arm2), sd2i = as.numeric(sd_arm2), 
                                       n2i = as.numeric(n_arm2), slab = study,
                                       data = dat)}) -> dat.long.smd
    dat.long.smd <- dat.long.smd[, !duplicated(colnames(dat.long.smd))]
    
    rbind(
      dat.long.smd[dat.long.smd$group==1,] %>% 
        dplyr::rename(condition = "condition_arm1") %>% 
        dplyr::select(-condition_arm2),
      dat.long.smd[dat.long.smd$group==2,] %>% 
        dplyr::rename(condition = "condition_arm2") %>% 
        dplyr::select(-condition_arm1)
    ) %>% {.[order(.$study),]} -> dat.long.smd
    
    # Re-generate clean slabs
    within(dat.long.smd, {
      study.clean = stringr::str_remove_all(study, "\\.[0-9]")
    }) %>% dplyr::distinct(study.clean, condition, 
                           mean, .keep_all = TRUE) -> dat.long.smd
    
    try({netmeta::pairwise(condition, mean = mean, sd = sd, n = n,
                           studlab = study.clean, data = dat.long.smd,
                           sm = "SMD")}, silent = TRUE) -> pw
    
    # Send error if pairwise fails
    if (class(pw)[1] == "try-error") {
      message = paste(
        "To transform the data, each trial arm is only allowed",
        "to provide", 
        crayon::bold(crayon::cyan("exactly one outcome.")), 
        crayon::bold(crayon::cyan("\n\n[!] Please make sure that:")),
        "\n- variables",
        catchName("condition_arm1", mC, data, parentFrame), 
        "and", 
        catchName("condition_arm2", mC, data, parentFrame), 
        "have unique labels in multiarm trials.",
        "\n- only one outcome (i.e., one instrument at one",
        "specific assessment point) is used.\n \n",
        attr(pw, "condition")$message)
      stop(message, call. = FALSE)
    }
    
  } else {
    pw = NULL
  }
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Dichotomous outcome module                        #
  # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  if (!any(is.null(event_arm1), is.null(event_arm2), 
           is.null(totaln_arm1), is.null(totaln_arm2))) {

    
    # Combine all variables in data frame
    dat = cbind(totaln_arm1, totaln_arm2, event_arm1, event_arm2,
                condition_arm1, condition_arm2, 
                study, extra.dat)
    
    # Pivot data frame to long format
    suppressWarnings({metafor::to.long("OR", 
                                       ai = as.numeric(event_arm1), 
                                       n1i = as.numeric(totaln_arm1),
                                       ci = as.numeric(event_arm2),
                                       n2i = as.numeric(totaln_arm2),
                                       slab = study,
                                       data = dat)}) -> dat.long.smd
    dat.long.smd <- dat.long.smd[, !duplicated(colnames(dat.long.smd))]
    
    rbind(
      dat.long.smd[dat.long.smd$group==1,] %>% 
        dplyr::rename(condition = "condition_arm1") %>% 
        dplyr::select(-condition_arm2),
      dat.long.smd[dat.long.smd$group==2,] %>% 
        dplyr::rename(condition = "condition_arm2") %>% 
        dplyr::select(-condition_arm1)
    ) %>% {.[order(.$study),]} -> dat.long.smd
    
    # Re-generate clean slabs
    within(dat.long.smd, {
      study.clean = stringr::str_remove_all(study, "\\.[0-9]")
    }) %>% dplyr::distinct(study.clean, condition, 
                           out1, .keep_all = TRUE) -> dat.long.smd
    
    try({netmeta::pairwise(condition, event = out1, n = out1 + out2,
                           studlab = study.clean, data = dat.long.smd,
                           sm = "OR")}, silent = TRUE) -> pw.dich
    
    # Send error if pairwise fails
    if (class(pw.dich)[1] == "try-error") {
      message = paste(
        "To transform the data, each trial arm is only allowed",
        "to provide", 
        crayon::bold(crayon::cyan("exactly one outcome.")), 
        crayon::bold(crayon::cyan("\n\n[!] Please make sure that:")),
        "\n- variables",
        catchName("condition_arm1", mC, data, parentFrame), 
        "and", 
        catchName("condition_arm2", mC, data, parentFrame), 
        "have unique labels in multiarm trials.",
        "\n- only one outcome (i.e., one instrument at one",
        "specific assessment point) is used.\n \n",
        attr(pw.dich, "condition")$message)
      stop(message, call. = FALSE)
    }
    
    # Create SMD using Chinn transformation
    within(pw.dich, {
      studlab = paste(studlab, "(reponse)")
      TE = (TE/(pi/sqrt(3)))*-1
      seTE = sqrt(seTE^2/((pi^2)/3))
    }) -> pw.dich
    
  } else {
    pw.dich = NULL
  }
  
  # Bind both results together
  pw = dplyr::bind_rows(pw, pw.dich)
  
  # Return
  if (nrow(pw) > 0) {
    miss.studies = unique(study)[!unique(study) %in% unique(pw$study)]
    if (length(miss.studies) > 0) {
      message("- ", crayon::yellow("[!] "), "Effect sizes could not be ",
              "calculated for the following studies:\n", 
              paste(miss.studies, collapse = " - ")) 
    }
    selvars = c("studlab", "treat1", "treat2", "TE", "seTE", "n1", "mean1", "sd1",
                "n2", "mean2", "sd2", "event1", "event2", "n1", "n2")
    if (length(extra.vars)>0) {
      res = dplyr::select(pw, dplyr::any_of(selvars), starts_with(extra.vars))
    } else { res = dplyr::select(pw, dplyr::any_of(selvars)) }
    attributes(res)$removed.studies = miss.studies
    
    return(dplyr::arrange(res, studlab))
  } else {
    stop("- ", crayon::yellow("[!] "), 
            "Effect sizes could not be ",
            "calculated for any study. ",
            "Did you provide the correct effect size data?") 
  }
  
}








