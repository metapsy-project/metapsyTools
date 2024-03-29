% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createRobSummary.R
\name{createRobSummary}
\alias{createRobSummary}
\title{Create a summary risk of bias plot}
\usage{
createRobSummary(model, 
                 name.low, 
                 name.high, 
                 name.unclear, 
                 which.run = model$which.run[1])
}
\arguments{
\item{model}{An object of class \code{runMetaAnalysis}, created by the \code{\link[=runMetaAnalysis]{runMetaAnalysis()}} function.}

\item{name.low}{A \code{character} vector, specifying which code(s) have been used in the original data for
studies with a low risk of bias.}

\item{name.high}{A \code{character} vector, specifying which code(s) have been used in the original data for
studies with a high risk of bias.}

\item{name.unclear}{A \code{character} vector, specifying which code(s) have been used in the original data for
studies with unclear risk of bias.}

\item{which.run}{The model in \code{model} that should be used for the summary risk of bias plot.
Uses the default analysis in \code{model} if no value is specified by the user. Possible values are
\code{"overall"}, \code{"combined"}, \code{"lowest"}, \code{"highest"}, \code{"outliers"},
\code{"influence"} and \code{"rob"}.}
}
\value{
Creates a RevMan-type risk of bias summary plot.
}
\description{
If the \code{rob.data} argument has been specified, this function allows to create
a summary risk of bias plot for results of the \code{\link[=runMetaAnalysis]{runMetaAnalysis()}} function.
}
\examples{
\dontrun{


# Define ROB data to be added to the models
robData = list(
  # Names of ROB variables included in 'data'
  domains = c("sg", "ac", "ba", "itt"),
  # Long-format labels for each ROB domain
  domain.names = c("Sequence Generation", 
                   "Allocation Concealment", 
                   "Blinding of Assessors", 
                   "ITT Analyses"),
  # Codes used to rate the risk of bias (sr=self-report)
  categories = c("0", "1", "sr"),
  # Symbols that should be used for these codes in forest plots
  symbols = c("-", "+", "s"),
  # Colors to be used in forest plots for each of these codes
  colors = c("red", "green", "yellow"))

# Run meta-analyses with ROB data
res <- depressionPsyCtr \%>\% 
  filterPoolingData(condition_arm1 \%in\% c("cbt", "pst", "3rd")) \%>\% 
  runMetaAnalysis(rob.data = robData)

# Create a summary plot
createRobSummary(res, 
                 name.low = "1", 
                 name.high = "0", 
                 name.unclear = "sr")

# Create a summary plot for the "combined" model
# - Recode 'sr' (self-report) as low risk of bias
createRobSummary(res, 
                 name.low = c("1", "sr"), 
                 name.high = "0", 
                 name.unclear = NULL,
                 which.run = "combined")
                 
                 
}

}
\seealso{
\link[dmetar:rob.summary]{dmetar::rob.summary}, \link[robvis:rob_summary]{robvis::rob_summary}, \link[meta:rob]{meta::rob}
}
\author{
Mathias Harrer \email{mathias.h.harrer@gmail.com},
Paula Kuper \email{paula.r.kuper@gmail.com}, Pim Cuijpers \email{p.cuijpers@vu.nl}
}
