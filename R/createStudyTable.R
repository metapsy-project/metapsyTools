#' Create study table
#'
#' This function creates an overview table containing selected study information.
#'
#' @param .data Meta-analysis study data; typically data created by the
#' \code{\link{expandMultiarmTrials}} or \code{\link{calculateEffectSizes}} function.
#' Trial arm-specific information (e.g. sample size in each group) can be added via
#' \code{\link{addTrialArmInfo}}. See 'Details'.
#' @param ... <\link[dplyr]{dplyr_data_masking}>. The name of several columns
#' (included in \code{.data}) that should be added to the study table. Also allows to
#' alter individual values/factor labels within a variable. See 'Details'.
#' @param .round.by.digits named \code{list}. Should contain the number of digits by which to
#' round a numeric column in \code{.data}. The name of the column must be specified in the
#' list element's name. Set to \code{NULL} if no rounding should be performed (default).
#' @param .column.names named \code{list}. If variable names should be renamed when producing the
#' study table, the new name should be included in this list. The original column name must
#' be specified as the name of the list element. Set to \code{NULL} if no renaming should
#' be performed (default).
#' @param .na.replace \code{character} to replace \code{NA} values with; \code{"nr"} by
#' default.
#' @param .html \code{logical}. Should an HTML table be produced? \code{TRUE} by default. See
#' 'Details'.
#'
#' @usage createStudyTable(.data, ...,
#'                  .round.by.digits = NULL,
#'                  .column.names = NULL,
#'                  .na.replace = "nr",
#'                  .html = TRUE)
#'
#' @return Returns a \code{data.frame} with all the selected variables. If \code{.html} is
#' \code{TRUE}, an HTML table will also be produced.
#'
#' @details \strong{General Purpose}: This function allows to select variables to be included in a study table.
#' Such study tables are typically part of a meta-analysis report/article. Variables are
#' included by adding their names to the function call and separating
#' them with commas. The columns will appear in the exact same order as specified in
#' the function.
#'
#' \strong{Trial-Arm Variables}: Before producing the final table, \code{createStudyTable} will filter out all redundant
#' rows based on the selected variables. If you want to include information that differs
#' between the (two or more) trial arms (e.g. the sample size of each group) as separate columns, you have to
#' use \code{\link{addTrialArmInfo}} first. This ensures that individual columns are created
#' for both the intervention and control group (e.g. \code{N_ig} and \code{N_cg}),
#' which can then be included in the call to \code{createStudyTable}.
#'
#'
#' \strong{Changing Values}: The function also allows to change specified values within a variable. Factor levels
#' encoded as numbers, for example (e.g. \code{country = 1} for European studies, and so
#' forth) can be changed by adding a concatenated (\code{\link[base]{c}}) vector to the
#' name of the variable. This vector should contain the \emph{new} value as a \code{character}
#' on the \emph{left} side, and the \emph{old} value on the \emph{right} side,
#' separated by '\code{=}' (e.g. \code{country = c("Europe" = "1")}).
#' The values will then be recoded before producing the table.
#'
#' \strong{HTML Table}: By default, \code{createStudyTable} produces an HTML table using
#' \code{\link[knitr]{kable}}. The HTML table makes copy & paste easier, particularly
#' when working with MS Word, since the table formatting is kept.
#'
#' For more details see the help vignette: \code{vignette("metapsyTools")}.
#'
#' @examples
#' \dontrun{
#' # Filter out all primary outcomes, check data,
#' # expand multiarm trials, calculate effect sizes,
#' # expand information that differs between trial arms;
#' # then produce study table using selected information.
#'
#' inpatients %>%
#'   filter(primary == 1) %>%
#'   checkDataFormat() %>%
#'   checkConflicts() %>%
#'   expandMultiarmTrials() %>%
#'   calculateEffectSizes() %>%
#'   createStudyTable(
#'     study,
#'     diag = c("Cutoff" = "3", "Mood" = "2",
#'              "MDD" = "1"),
#'     agegrp2,
#'     ADD_setting = c("Nursing" = "nursing",
#'                     "Other" = "other",
#'                     "Psychiatric" = "psych"),
#'     meanage, propwomen,
#'     type_trt1 = c("CBT" = "cbt", "Other" = "other", "PST" = "pst",
#'                   "BA" = "bat", "LR" = "lrt",
#'                   "PDT" = "dyn", "IPT" = "ipt"),
#'     nsess_trt1, Post_N_trt1, Post_N_trt2, Cond_spec_trt1, Cond_spec_trt2,
#'     country = c("Canada" = "4", "Europe" = "3", "USA" = "1",
#'                 "Asia" = "6", "Middle East" = "7", "Australia" = "5"),
#'     sg, ac, ba, itt,
#'     .round.by.digits = list(meanage = 0, Post_N_trt1 = 0,
#'                             Post_N_trt2 = 0),
#'     .column.names = list(agegrp2 = "age group",
#'                          ADD_setting = "setting",
#'                          Post_N_trt1 = "N_ig", Post_N_trt2 = "N_cg",
#'                          Cond_spec_trt1 = "format_ig",
#'                          Cond_spec_trt2 = "format_cg")) -> table
#' }
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @seealso \code{\link{addTrialArmInfo}}, \code{\link[kableExtra]{kable_styling}}
#' @importFrom dplyr distinct enquos as_label all_of select filter
#' @importFrom tidyr replace_na
#' @importFrom purrr map map_chr map_dfr map2_dfr
#' @importFrom rlang eval_tidy
#' @importFrom forcats fct_recode
#' @importFrom kableExtra kable_styling column_spec
#' @importFrom knitr kable
#' @importFrom stats dffits model.matrix rnorm rstudent
#' @importFrom utils combn
#' @importFrom crayon green yellow cyan bold
#' @export createStudyTable

createStudyTable = function(.data, ...,
                            .round.by.digits = NULL,
                            .column.names = NULL,
                            .na.replace = "nr",
                            .html = TRUE){

  if (class(.data)[1] == "runMetaAnalysis")
    .data = .data$data

  # Get tidy-select columns
  convert.vars = dplyr::enquos(...) %>%
    purrr::map(~ dplyr::as_label(.)) %>% names() %>% {!(. == "")}

  # Select non-convert variables
  non.convert.vars = dplyr::enquos(...) %>%
    {.[!convert.vars]} %>% purrr::map_chr(~ dplyr::as_label(.))

  vars = dplyr::enquos(...) %>% names()

  args = dplyr::enquos(...) %>% {.[convert.vars]} %>%
    purrr::map(~ rlang::eval_tidy(.))

  # Converted data
  .data %>%
    dplyr::select(dplyr::all_of(names(args))) %>%
    purrr::map_dfr(~ as.character(.) %>% as.factor()) %>%
    as.list() %>%
    purrr::map2_dfr(., args,
       function(x,y) forcats::fct_recode(x, !!!y)) -> data.convert

  # Merge with nonconverted data
  vars[vars == ""] = non.convert.vars
  if (ncol(.data[non.convert.vars]) > 0 & ncol(data.convert) > 0){
    cbind(.data[non.convert.vars], data.convert) %>%
      dplyr::select(dplyr::all_of(vars)) -> data
  } else if (ncol(.data[non.convert.vars]) > 0){
    cbind(.data[non.convert.vars]) %>%
      dplyr::select(dplyr::all_of(vars)) -> data
  } else if (ncol(data.convert) > 0){
    cbind(data.convert) %>%
      dplyr::select(dplyr::all_of(vars)) -> data}

  if (!is.null(.round.by.digits)){

    # Round digits
    data %>% dplyr::select(dplyr::all_of(names(.round.by.digits))) %>%
      as.list() %>% purrr::map2_dfr(.round.by.digits, function(x,y){
        suppressWarnings({as.numeric(x)}) -> x
        round(x, y)
      }) -> rounded

    # Replace
    for (i in 1:length(rounded)){
      data[[names(rounded)[i]]] == "NA" |
        is.na(data[[names(rounded)[i]]]) -> mask
      data[[names(rounded)[i]]][!mask] = rounded[[i]][!mask]
    }
  }

  # Remove redundant rows
  data.d = dplyr::distinct(data)


  # Replace NA values
  apply(data.d, 2, function(x){
    x[x == "NA"] = NA
    tidyr::replace_na(x, .na.replace) -> x
    return(x)}) %>%
    data.frame() -> data.d

  # Rename columns
  if (!is.null(.column.names)){
    colnames(data.d)[which(colnames(data.d) %in% names(.column.names))] -> col.mask
    as.data.frame(.column.names) %>%
      dplyr::select(dplyr::all_of(col.mask)) %>%
      {.[1,]} %>% as.character() -> new.cnames
    colnames(data.d)[which(colnames(data.d) %in% names(.column.names))] = new.cnames
  }

  rownames(data.d) = NULL

  if (.html == TRUE){
    message(crayon::cyan(crayon::bold("- Creating HTML table...")))
    data.d %>% knitr::kable() %>%
      kableExtra::kable_styling(font_size = 8, full_width = FALSE) %>%
      kableExtra::column_spec(1, bold = TRUE, width_min = "13em") %>% print()
  }

  message("- ", crayon::green("[OK] "), "Study table created successfully.")
  return(data.d)

}


