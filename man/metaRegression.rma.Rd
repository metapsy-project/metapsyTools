% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaRegression.R
\name{metaRegression.rma}
\alias{metaRegression.rma}
\title{Meta-Regression method for objects of class 'runMetaAnalysis'}
\usage{
\method{metaRegression}{rma}(x, formula = NULL, ...)
}
\arguments{
\item{x}{A model extracted from an object of class \code{runMetaAnalysis}.}

\item{formula}{A \code{formula} object describing the predictor(s) to
be added to the model. Default is \code{NULL}.}

\item{...}{Additional arguments.}
}
\description{
Serves as a wrapper for \code{update.rma}.
}
\author{
Mathias Harrer \email{mathias.h.harrer@gmail.com},
Paula Kuper \email{paula.r.kuper@gmail.com}, Pim Cuijpers \email{p.cuijpers@vu.nl}
}
