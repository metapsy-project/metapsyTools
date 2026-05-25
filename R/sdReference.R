# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# SD reference standardization for metapsyTools.                          #
#                                                                         #
# Exports : listSDReferences()                                            #
# Internal: applySDReference()                                            #
#                                                                         #
# The bundled data object `sd.reference` (see data-raw/build-sd-reference #
# .R) provides one row per psychometric instrument with the headline      #
# random-effects pooled SD from the Metapsy SD-reference catalogue.       #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

utils::globalVariables(c("sdReference"))


#' Display available SD-reference instruments
#'
#' Prints the catalogue of patient or clinician-rated outcome 
#' measurement instruments for which a
#' reference standard deviation is available within `metapsyTools`. The
#' shorthand in the `instr` column is the lookup key used by
#' \code{\link{calculateEffectSizes}} when \code{sd.reference = "fill"}
#' or \code{"override"}: by default, this code is matched (case-insensitively)
#' against the `instrument` column of the data.
#'
#' @param pattern \code{character}. Optional regular-expression filter
#'   applied to both the \code{instr} shorthand and the \code{full}
#'   instrument name. Matching is case-insensitive. \code{NULL}
#'   (default) returns the entire catalogue.
#' @param print \code{logical}. If \code{TRUE} (default), the catalogue
#'   is pretty-printed to the console. The (possibly filtered) table is
#'   always returned invisibly so it can be captured into a variable.
#'
#' @return Invisibly, a \code{data.frame} with one row per instrument and
#' the columns \code{instr}, \code{full}, \code{k}, \code{k_es},
#' \code{sd}, \code{lo}, \code{hi}, \code{tau2}, \code{i2}. The
#' \code{"version"}, \code{"version_date"}, \code{"source_url"} and
#' \code{"citation"} attributes of the bundled catalogue are preserved.
#'
#' @details
#' The bundled catalogue is derived from the Metapsy SD-reference
#' database and keeps only the headline random-effects pooled SD per
#' instrument (no subgroup, no sensitivity analysis). For subgroup- or
#' sensitivity-specific estimates, browse the full reference database at
#' [metapsy.org/database/sd-reference.html](https://metapsy.org/database/sd-reference.html).
#'
#' @examples
#' \dontrun{
#' # Full catalogue
#' listSDReferences()
#'
#' # Filter
#' listSDReferences("depression")
#' listSDReferences("^BDI")
#'
#' # Capture without printing
#' tbl <- listSDReferences(print = FALSE)
#' }
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com}
#' @seealso \code{\link{calculateEffectSizes}}
#' @export listSDReferences
listSDReferences = function(pattern = NULL, print = TRUE){
  
  tbl.full = sdReference
  tbl = tbl.full
  
  if (!is.null(pattern)){
    keep = grepl(pattern, tbl$instr, ignore.case = TRUE) |
           grepl(pattern, tbl$full, ignore.case = TRUE)
    tbl = tbl[keep, , drop = FALSE]
  }
  # Preserve attributes carried by the bundled object
  for (a in c("version", "version_date", "source_url",
              "citation", "built_at")) {
    attr(tbl, a) = attr(tbl.full, a)
  }
  
  if (!print) return(invisible(tbl))
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Pretty-print to the console                                         #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  ver = attr(tbl.full, "version")
  vdate = attr(tbl.full, "version_date")
  n.all = nrow(tbl.full)
  n.shown = nrow(tbl)
  
  # Header
  cat("\n")
  cat(crayon::bold(crayon::cyan("  Metapsy SD-Reference Catalogue")), "\n", sep = "")
  cat(crayon::silver(paste0("  Version ", ver, " (", vdate, ")  |  ",
                            n.all, " instruments")),
      "\n", sep = "")
  cat(crayon::silver("  https://www.metapsy.org/tools/standardizers"),
      "\n\n", sep = "")
  
  if (n.shown == 0){
    cat(crayon::yellow(paste0("  No instruments match pattern '", pattern,
                              "'.\n\n")))
    return(invisible(tbl))
  }
  if (!is.null(pattern)){
    cat(crayon::silver(paste0("  Showing ", n.shown, " of ", n.all,
                              " (pattern: '", pattern, "')")),
        "\n\n", sep = "")
  }
  
  # Column widths (sd-fixed 8, k-fixed 4, k_es-fixed 5)
  w.instr = max(nchar(tbl$instr), nchar("instr"))
  full.max = 50L  # cap full-name width for readability
  w.full = min(max(nchar(tbl$full), nchar("instrument")), full.max)
  
  # Column header row
  cat("  ",
      crayon::bold(formatC("instr", width = w.instr, flag = "-")), "  ",
      crayon::bold(formatC("instrument", width = w.full, flag = "-")), "  ",
      crayon::bold(formatC("sd", width = 8)), "  ",
      crayon::bold(formatC("k", width = 4)), "  ",
      crayon::bold(formatC("estimates", width = 9)),
      "\n", sep = "")
  cat("  ",
      strrep("-", w.instr + 2 + w.full + 2 + 8 + 2 + 4 + 2 + 9),
      "\n", sep = "")
  
  # Data rows
  for (i in seq_len(n.shown)){
    full.disp = tbl$full[i]
    if (nchar(full.disp) > w.full){
      full.disp = paste0(substr(full.disp, 1, w.full - 1), "\u2026")
    }
    sd.disp = formatC(tbl$sd[i], width = 8, format = "f", digits = 3)
    cat("  ",
        crayon::cyan(formatC(tbl$instr[i], width = w.instr, flag = "-")), "  ",
        formatC(full.disp, width = w.full, flag = "-"), "  ",
        crayon::green(sd.disp), "  ",
        formatC(tbl$k[i], width = 4), "  ",
        formatC(tbl$k_es[i], width = 5),
        "\n", sep = "")
  }
  
  # Footer
  cat("\n")
  cat(crayon::silver(
        "  Use the 'instr' shorthand in your data's `instrument` column."),
      "\n", sep = "")
  cat(crayon::silver(
        "  Pass sd.reference = \"fill\" or \"override\" to calculateEffectSizes()."),
      "\n\n", sep = "")
  
  invisible(tbl)
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Internal helper: apply reference SDs to a Metapsy-format data.frame.    #
#                                                                         #
# mode = "none"     -> no-op (returns data unchanged)                     #
# mode = "fill"     -> only fill sd_arm1/sd_arm2 cells that are NA, where #
#                      mean_arm1/2 and n_arm1/2 are non-NA                #
# mode = "override" -> set sd_arm1/sd_arm2 to the reference for every     #
#                      row whose instrument matches the catalogue and    #
#                      where mean_arm1/2 and n_arm1/2 are non-NA          #
#                                                                         #
# Only touches the cross-sectional path (sd_arm1, sd_arm2). The change-   #
# score SDs (sd_change_arm1/2) are NEVER modified, since the catalogue    #
# refers to cross-sectional SDs of the scale.                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

applySDReference = function(data, mode = c("none", "fill", "override"),
                            instr.col = "instrument", sd.ref = NULL){
  mode = match.arg(mode)
  if (mode == "none") return(data)
  
  if (!instr.col %in% colnames(data)){
    message("- ", crayon::yellow("[!] "),
            "Column '", instr.col, "' not found in data. ",
            "Reference SD standardization skipped.")
    return(data)
  }
  required.cols = c("mean_arm1", "mean_arm2", "n_arm1", "n_arm2",
                    "sd_arm1", "sd_arm2")
  missing.cols = setdiff(required.cols, colnames(data))
  if (length(missing.cols) > 0){
    message("- ", crayon::yellow("[!] "),
            "Required column(s) missing: ", paste(missing.cols, collapse = ", "),
            ". Reference SD standardization skipped.")
    return(data)
  }
  
  if (is.null(sd.ref)) sd.ref = sdReference
  
  # Case-insensitive lookup against the instrument shorthand
  lookup = setNames(sd.ref$sd, tolower(sd.ref$instr))
  ref.sd = unname(lookup[tolower(as.character(data[[instr.col]]))])
  
  has.match = !is.na(ref.sd)
  has.mean = !is.na(data$mean_arm1) & !is.na(data$mean_arm2) &
             !is.na(data$n_arm1) & !is.na(data$n_arm2)
  
  changed.rows = logical(nrow(data))
  
  if (mode == "fill"){
    fill.arm1 = has.match & has.mean & is.na(data$sd_arm1)
    fill.arm2 = has.match & has.mean & is.na(data$sd_arm2)
    data$sd_arm1[fill.arm1] = ref.sd[fill.arm1]
    data$sd_arm2[fill.arm2] = ref.sd[fill.arm2]
    changed.rows = fill.arm1 | fill.arm2
    # Rows that could have benefited but did not match (had mean+n,
    # at least one SD missing, no catalogue match)
    need.sub = has.mean & (is.na(data$sd_arm1) | is.na(data$sd_arm2))
  } else if (mode == "override"){
    target = has.match & has.mean
    data$sd_arm1[target] = ref.sd[target]
    data$sd_arm2[target] = ref.sd[target]
    changed.rows = target
    need.sub = has.mean
  }
  
  # Identify unique instrument values in rows that needed substitution but
  # did not match the catalogue. Used below for the diagnostic.
  unmatched.vals = unique(as.character(data[[instr.col]][need.sub & !has.match]))
  unmatched.vals = unmatched.vals[!is.na(unmatched.vals) & nzchar(unmatched.vals)]
  unmatched.vals = sort(unmatched.vals)
  n.unmatched = length(unmatched.vals)
  
  format.unmatched = function(vals, n.max = 10){
    if (length(vals) > n.max) {
      paste0(paste(vals[seq_len(n.max)], collapse = ", "),
             ", ... (+", length(vals) - n.max, " more)")
    } else {
      paste(vals, collapse = ", ")
    }
  }
  
  if (any(changed.rows)){
    tab = sort(table(as.character(data[[instr.col]][changed.rows])),
               decreasing = TRUE)
    ver = attr(sd.ref, "version")
    if (is.null(ver)) ver = "unknown"
    message("- ", crayon::green("[OK] "),
            "Reference SDs applied to ", sum(changed.rows),
            " row(s) (mode = '", mode, "', catalogue v", ver, "). ",
            "Instruments: ",
            paste0(names(tab), " (", as.integer(tab), ")", collapse = ", "))
    # Inform about rows that needed substitution but did not match
    if (n.unmatched > 0){
      message("- ", crayon::yellow("[!] "),
              sum(need.sub & !has.match), " row(s) needed substitution but the '",
              instr.col, "' value was not in the catalogue. ",
              "Unique unmatched value(s): ", format.unmatched(unmatched.vals), ".")
      message("  ", crayon::silver(
              "    See listSDReferences() for the available shorthand codes."))
    }
  } else {
    if (n.unmatched > 0){
      message("- ", crayon::yellow("[!] "),
              "No rows matched the SD reference catalogue (mode = '", mode, "'). ",
              "Unique '", instr.col, "' value(s) in rows needing substitution: ",
              format.unmatched(unmatched.vals), ".")
      message("  ", crayon::silver(
              "    See listSDReferences() for the available shorthand codes."))
    } else {
      message("- ", crayon::yellow("[!] "),
              "No rows in the data needed SD substitution (mode = '", mode, "').")
    }
  }
  data
}