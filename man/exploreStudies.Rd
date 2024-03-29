% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exploreStudies.R
\name{exploreStudies}
\alias{exploreStudies}
\title{Explore included treatments and comparisons}
\usage{
exploreStudies(data,
               which = c("treatments", "comparisons"),
               
               # Metapsy standard variables
               .study.var = "study",
               .condition = "condition",
               .condition.specification = "multi",
               .groups.column.indicator = c("_arm1", "_arm2"),
               .trt.indicator = "arm",
               .n.vars = c("n", "n_change", "totaln", "N"),
               
               # Output
               html = TRUE)
}
\arguments{
\item{data}{\code{data.frame}. Effect size data in the wide format, as
created by \code{\link{calculateEffectSizes}}. For the other default settings to be
applicable, the data set should follow the \href{https://docs.metapsy.org/data-preparation/format/}{Metapsy data standard}.
Alternatively, one can also provide an \code{metapsyDatabase} object as returned by \code{\link[metapsyData:getData]{metapsyData::getData()}},
or a meta-analysis object returned by \code{\link[=runMetaAnalysis]{runMetaAnalysis()}}.}

\item{which}{Should the data set be summarized with respect to the included
treatments (\code{"treatments"}) or treatment comparisons (\code{"comparisons"})? Defaults to
\code{"treatments"}.}

\item{.study.var}{\code{character}. The name of the variable in the data set in which
the study labels are stored.}

\item{.condition}{\code{character}. The prefix of the two variables in \code{data} in
which the conditions (e.g. "guided iCBT", "waitlist") of the trial arm comparison are stored.}

\item{.condition.specification}{\code{character}. The prefix of the two variables in the dataset
which provide a "specification" of the trial arm condition in multiarm trials.}

\item{.groups.column.indicator}{\code{character}. A character vector with two elements,
representing the suffix used to differentiate between the first and second arm in a comparison.}

\item{.trt.indicator}{\code{character}. A character specifying the name used to indicate the treatment arm.}

\item{.n.vars}{\code{character}. A character vector which includes the names of all variables in the data set in
which sample size information is stored. Only the prefix is needed, where \code{.groups.column.indicator}
provides the suffixes.}

\item{html}{\code{logical}. Should an HTML table be created for the results? Default is \code{TRUE}.}
}
\value{
Returns an object of class \code{"exploreStudies"}. This object includes a list object called \code{summary}
in which the counts for distinct treatments (\code{conditions}) and comparisons (\code{comparisons}) are summarized, as well as
a \code{data.frame}  data. This data frame includes the initially provided data set collapsed by study
(so that each row represents one study). To this data set, variables are added that encode how many arms with a specific condition
are included in the trial (e.g. if \code{cbt=2}, this means that two CBT groups are included in the trial),
as well as the number of distinct comparisons, and the sample size of both (these columns all start with \code{n.}).
This can be helpful to perform further descriptive analyses.
}
\description{
This function allows to summarize included treatments and
treatment comparisons in a data set.
}
\details{
Using the variables provided in the \code{.n.vars} argument, \code{exploreStudies} calculates the arm- and study-specific
sample sizes. If no adequate information is provided, sample sizes cannot be calculated for a study. If this
is the case, a warning is printed, pointing to the studies with missing sample size information.
}
\examples{
\dontrun{
# Explore studies in built-in dataset
data("depressionPsyCtr")
exploreStudies(depressionPsyCtr, "treatments") 
exploreStudies(depressionPsyCtr, "comparisons") 

# - Extract metapsy database using metapsyData
# - Filter CBT and PST studies
# - Run a meta-analysis and explore synthesize studies
library(metapsyData)
getData("depression-psyctr", version="22.0.2") \%>\% 
  filterPoolingData(condition_arm1 \%in\% c("cbt", "pst")) \%>\% 
  runMetaAnalysis(which.run = c("combined")) -> res

exploreStudies(res)
exploreStudies(res, "comparisons")
}

}
\seealso{
\code{\link{createStudyTable}}, \code{\link{calculateEffectSizes}},
\code{\link{subgroupAnalysis}}, \code{\link{correctPublicationBias}},
\code{\link{metaRegression}}, \code{\link{runMetaAnalysis}}.
}
\author{
Mathias Harrer \email{mathias.h.harrer@gmail.com},
Paula Kuper \email{paula.r.kuper@gmail.com}, Pim Cuijpers \email{p.cuijpers@vu.nl}
}
