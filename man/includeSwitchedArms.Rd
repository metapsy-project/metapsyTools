% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{includeSwitchedArms}
\alias{includeSwitchedArms}
\title{Include information of rows with switched reference group}
\usage{
includeSwitchedArms(dat, ...)
}
\arguments{
\item{dat}{Data set created by \code{calculateEffectSizes} in the wider format, which includes
calculated effect sizes and standard errors in columns \code{es} and \code{se}, respectively. Only
data sets created by \code{calculateEffectSizes} with \code{trt.indicator} set to \code{"trt"}
can be used.}

\item{...}{Further arguments (not used).}
}
\description{
Adds effect size data and study information of rows with switched reference arms.
}
\keyword{internal}
