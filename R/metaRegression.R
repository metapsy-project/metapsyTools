#' Meta-Regression method for objects of class 'runMetaAnalysis'
#'
#' Serves as a wrapper for \code{metareg} or \code{update.rma}, depending
#' on the class of the fitted model.
#'
#' @param x A model extracted from an object of class \code{runMetaAnalysis}.
#' @param ... Additional arguments.
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#' @examples
#' \dontrun{
#' metaRegression(res$model.combined, ~ rob + scale(year))
#' }
#' @export metaRegression

metaRegression = function(x, ...){
  UseMethod("metaRegression", x)
}


#' Meta-Regression method for objects of class 'runMetaAnalysis'
#'
#' Serves as a wrapper for \code{metareg}.
#'
#' @param x A model extracted from an object of class \code{runMetaAnalysis}.
#' @param formula A \code{formula} object describing the predictor(s) to
#' be added to the model. Default is \code{NULL}.
#' @param ... Additional arguments.
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @importFrom meta metareg
#' @export
#' @method metaRegression meta

metaRegression.meta = function(x, formula = NULL, ...){
  if (identical(x$.type.es, "RR")){
    message("- ", crayon::green("[OK] "), 
     "Coefficent estimates are based on log-transformed risk ratios.")
    message("- ", crayon::green("[OK] "), 
     "Use the exp() function to transform estimates into regular RR values.")
  }
  m.rerun = meta::metareg(x, as.formula(formula), ...)
  return(m.rerun)
}



#' Meta-Regression method for objects of class 'runMetaAnalysis'
#'
#' Prints information that models have to be extracted from \code{runMetaAnalysis}
#' results objects to run meta-regression.s
#'
#' @param x An object of class \code{runMetaAnalysis}.
#' @param formula not used.
#' @param ... Additional arguments.
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @export
#' @method metaRegression runMetaAnalysis


metaRegression.runMetaAnalysis = function(x, formula = NULL, ...){
  stop("It seems like you have supplied a 'runMetaAnalysis' result. ",
       "To run a meta-regression, please extract a specific model ",
       "(e.g. x$model.overall).")
}



#' Meta-Regression method for objects of class 'runMetaAnalysis'
#'
#' Serves as a wrapper for \code{update.rma}.
#'
#' @param x A model extracted from an object of class \code{runMetaAnalysis}.
#' @param formula A \code{formula} object describing the predictor(s) to
#' be added to the model. Default is \code{NULL}.
#' @param ... Additional arguments.
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Paula Kuper \email{paula.r.kuper@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#'
#' @importFrom metafor rma.mv
#' @export
#' @method metaRegression rma

metaRegression.rma = function(x, formula = NULL, ...){
  
  if (identical(x$.type.es, "RR")){
    message("- ", crayon::green("[OK] "), 
     "Coefficent estimates are based on log-transformed risk ratios.")
    message("- ", crayon::green("[OK] "), 
     "Use the exp() function to transform estimates into regular RR values.")
  }

  if (class(x$V)[1] == "dgCMatrix" |
      class(x$V)[1] == "dsCMatrix"){

    if (is.null(formula)){
      formula.new = x$formula.yi
    } else {
      formula.new = paste(as.character(x$formula.yi)[2],
                          "~", as.character(x$formula.yi)[3], "+",
                          as.character(formula)[2])}

    with(x, {metafor::rma.mv(yi = as.formula(formula.new), V = V,
                             test = test, data = legacy$data,
                             method = method,
                             random = legacy$formula.rnd,
                             sparse = TRUE, ...)}) -> m.rerun
    return(m.rerun)

  } else {

    if (is.null(formula)){
      formula.new = paste(x$legacy$es.var, "~ 1")
    } else {
      formula.new = paste(x$legacy$es.var,
                          "~ 1 +", as.character(formula)[2])
    }
    with(x, {metafor::rma.mv(as.formula(formula.new), V = V,
                             test = test, data = legacy$data,
                             method = method, random = legacy$formula.rnd,
                             sparse = TRUE, ...)}) -> m.rerun
    return(m.rerun)
  }
}

