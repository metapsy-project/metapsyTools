% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateEffectSizes.R
\name{calculateEffectSizes}
\alias{calculateEffectSizes}
\title{Calculate effect sizes}
\usage{
calculateEffectSizes(data,
                     funcs.g = list(g.m.sd = g.m.sd,
                                    g.change.m.sd = g.change.m.sd,
                                    g.binary = g.binary,
                                    g.precalc = g.precalc),
                     funcs.rr = list(rr.binary = rr.binary,
                                     rr.precalc = rr.precalc),
                     include.switched.arms = FALSE,
                     change.sign = NULL,
                     vars.for.id = c("study", "outcome_type",
                                     "instrument", "time",
                                     "time_weeks",
                                     "rating"),
                     .condition = "condition",
                     .condition.specification = "multi",
                     .groups.column.indicator = c("_arm1", "_arm2"),
                     .trt.indicator = "arm")
}
\arguments{
\item{data}{Meta-analysis data set formatted using the
\href{https://tools.metapsy.org/articles/metapsyTools.html#required-data-structure}{Metapsy guidelines}.}

\item{funcs.g}{\code{list} of functions. These functions will be used to calculate the
effect sizes (Hedges' \emph{g}) based on the raw data (see Details).}

\item{funcs.rr}{\code{list} of functions. These functions will be used to calculate risk ratios
based on the raw event data (see Details).}

\item{include.switched.arms}{\code{logical}. Should all unique arm \emph{comparisons} (in lieu of unique arm \emph{combinations}) be
calculated? Default is \code{FALSE}.}

\item{change.sign}{\code{character}. Name of a \code{logical} column in \code{data}, encoding if the
sign of a calculated effect size should be reversed (\code{TRUE}) or not (\code{FALSE}). Set to \code{NULL} (default) if
no changes should be made.}

\item{vars.for.id}{\code{character} vector, containing column names of all variables
used to construct unique comparison IDs.}

\item{.condition}{\code{character}. The prefix of the two variables in \code{data} in
which the conditions (e.g. "guided iCBT", "waitlist") of the trial arm comparison are stored.}

\item{.condition.specification}{\code{character}. The prefix of the two variables in the dataset
which provide a "specification" of the trial arm condition in multiarm trials.}

\item{.groups.column.indicator}{\code{character}. A character vector with two elements,
representing the suffix used to differentiate between the first and second arm in a comparison.}

\item{.trt.indicator}{\code{character}. A character specifying the name used to indicate the treatment arm.}
}
\value{
\code{calculateEffectSizes} returns the meta-analysis data set as
class \code{data.frame} in wide format (if results are saved to a variable).
It also generates the following columns, wich are added to the data:
\itemize{
\item \code{.id}: Unique identifier for a trial arm comparison/row.
\item \code{.g}: Calculated effect size (Hedges' \emph{g}).
\item \code{.g_se}: Standard error of Hedges' \emph{g}.
\item \code{.log_rr}: Calculated effect size (\emph{logRR}).
\item \code{.log_rr_se}: Standard error of \emph{logRR}.
\item \code{.event_arm1}: Number of events (responders, remission, deterioration cases) in the first trial arm.
\item \code{.event_arm2}: Number of events (responders, remission, deterioration cases) in the second trial arm.
\item \code{.totaln_arm1}: Total sample size in the first trial arm.
\item \code{.totaln_arm2}: Total sample size in the second trial arm.
}
}
\description{
This is a function to calculate effect sizes of meta-analysis data
prepared in the \href{https://docs.metapsy.org/data-preparation/format/}{Metapsy data format}.
}
\details{
By default, \code{calculateEffectSizes} calculates the
small-sample bias corrected standardized mean difference  (Hedges' \emph{g})
and log-risk ratios, as well their respective standard errors, if
adequate raw effect size data is available for a comparison.

It is essential that the data set in \code{data} contains a few required columns
for this to work. An overview of the required data format is provided
on the \href{https://tools.metapsy.org/articles/metapsytools}{"Get Started"}
page of the \code{metapsyTools} documentation.

\strong{Standardized mean differences (Hedges' \emph{g})} can be calculated from the following
column types:
\itemize{
\item (1) Continuous Outcome Data
\itemize{
\item \strong{\code{mean_arm1}}: Mean of the outcome in the first arm at the measured time point.
\item \strong{\code{mean_arm2}}: Mean of the outcome in the second arm at the measured time point.
\item \strong{\code{sd_arm1}}: Standard deviation of the outcome in the first arm at the measured time point.
\item \strong{\code{sd_arm2}}: Standard deviation of the outcome in the second arm at the measured time point.
\item \strong{\code{n_arm1}}: Sample size in the first trial arm.
\item \strong{\code{n_arm2}}: Sample size in the second trial arm.
}
\item (2) Change Score Data
\itemize{
\item \strong{\code{mean_change_arm1}}: Mean score change between baseline and the measured time point in the first arm.
\item \strong{\code{mean_change_arm2}}: Mean score change between baseline and the measured time point in the second arm.
\item \strong{\code{sd_change_arm1}}: Standard deviation of the mean change in the first arm.
\item \strong{\code{sd_change_arm2}}: Standard deviation of the mean change in the second arm.
\item \strong{\code{n_change_arm1}}: Sample size in the first trial arm.
\item \strong{\code{n_change_arm2}}: Sample size in the second trial arm.
}
\item (3) Dichotomous Outcome Data
\itemize{
\item \strong{\code{event_arm1}}: Number of events (responders, remission, deterioration cases) in the first trial arm.
\item \strong{\code{event_arm2}}: Number of events (responders, remission, deterioration cases) in the second trial arm.
\item \strong{\code{totaln_arm1}}: Sample size in the first trial arm.
\item \strong{\code{totaln_arm2}}: Sample size in the second trial arm.
}
\item (4) Pre-calculated Hedges' \emph{g}
\itemize{
\item \strong{\code{precalc_g}}: The pre-calculated value of Hedges' \emph{g} (small-sample bias corrected standardized mean difference; \href{https://journals.sagepub.com/doi/10.3102/10769986006002107}{Hedges, 1981}).
\item \strong{\code{precalc_g_se}}: Standard error of \emph{g}.
}
}

The \strong{log-risk ratio} and its standard error can be calculated from the followin column types:
\itemize{
\item (1) Dichotomous Outcome Data
\itemize{
\item \strong{\code{event_arm1}}: Number of events (responders, remission, deterioration cases) in the first trial arm.
\item \strong{\code{event_arm2}}: Number of events (responders, remission, deterioration cases) in the second trial arm.
\item \strong{\code{totaln_arm1}}: Sample size in the first trial arm.
\item \strong{\code{totaln_arm2}}: Sample size in the second trial arm.
}
\item (2) Pre-calculated log-risk ratio
\itemize{
\item \strong{\code{precalc_log_rr}}: The pre-calculated value of the log-risk ratio logRR, comparing events in the first arm to events in the second arm.
\item \strong{\code{precalc_log_rr_se}}: The standard error of the log-risk ratio logRR, comparing events in the first arm to events in the second arm.
}
}

Other functions can be added to the list provided to \code{funcs.g} and \code{funcs.rr}.
However, results of the function must result in a \code{data.frame} that contains
the following columns:
\itemize{
\item \code{.id}: Unique identifier for a trial arm comparison/row.
\item \code{.g}: Calculated effect size (Hedges' \emph{g}).
\item \code{.g_se}: Standard error of Hedges' \emph{g}.
\item \code{.log_rr}: Calculated effect size (logRR).
\item \code{.log_rr_se}: Standard error of logRR.
\item \code{.event_arm1}: Number of events (responders, remission, deterioration cases) in the first trial arm.
\item \code{.event_arm2}: Number of events (responders, remission, deterioration cases) in the second trial arm.
\item \code{.totaln_arm1}: Total sample size in the first trial arm.
\item \code{.totaln_arm2}: Total sample size in the second trial arm.
}

It is possible to set one or several of these column entries to \code{NA}; but the columns
themselves must be included.

For more details see the \href{https://tools.metapsy.org/articles/metapsytools}{Get Started} vignette.
}
\examples{
\dontrun{
data("depressionPsyCtr")
depressionPsyCtr \%>\%
    checkDataFormat() \%>\%
    checkConflicts() \%>\% 
    calculateEffectSizes()
}

}
\seealso{
\code{\link{checkDataFormat}}
}
\author{
Mathias Harrer \email{mathias.h.harrer@gmail.com},
Paula Kuper \email{paula.r.kuper@gmail.com},
Pim Cuijpers \email{p.cuijpers@vu.nl}
}
