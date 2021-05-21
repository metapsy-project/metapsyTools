#'  The 'inpatients' dataset with 32 clinical trials
#'
#' An example dataset containing data of 32 clinical trials. Its format equals the format of data in the "metapsy" database that is similar to the long format in R. The dataset also contains columns with study characteristics that are important for effect size calculation and more.
#'
#' @format A \code{data.frame} with 179 rows and 45 variables:
#' \describe{
#'   \item{study}{\code{Character}, The study label containing the author(s) and year of the study.}
#'   \item{condition}{\code{Character}, The condition of the groups, either "intervention group" or "control group".}
#'   \item{Cond_spec}{\code{Character}, The specific intervention in conditions.}
#'   \item{is.multiarm}{\code{Numeric}, The dichotomized indication if a study has multiple arms or not.}
#'   \item{no.arms}{\code{Numeric}, The number of arms of a study.}
#'   \item{multiple.arms}{\code{Character}, The specification of arms in multiarm studies.}
#'   \item{Outc_type}{\code{Character}, The type of outcome.}
#'   \item{primary}{\code{Numeric}, The indication if a outcome is primary or not.}
#'   \item{Outc_measure}{\code{Character}, The outcome measure used.}
#'   \item{Time}{\code{Character}, The dichotomized time of assessment, either "post" or "FU".}
#'   \item{Time_weeks}{\code{Character}, The assessment time of FU in weeks.}
#'   \item{year}{\code{Numeric}, The year of the study.}
#'   \item{country}{\code{Character}, The country of the study.}
#'   \item{...}{}
#' }
#' @source
#' @author Mathias Harrer, Paula Kuper, Pim Cuijpers
#' @usage data(inpatients)
"inpatients"
