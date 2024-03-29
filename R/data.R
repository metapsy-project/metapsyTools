#'  The 'depressionPsyCtr' dataset
#'
#' An example data set containing a subset of the 2021 depression trials database.
#' The data set has been formatted following the 
#' \href{https://docs.metapsy.org/data-preparation/format/}{Metapsy data standard}.
#' 
#' The dataset can be used by the \code{runMetaAnalysis} function without any
#' additional preparation required.
#'
#'
#' @format A \code{data.frame} with 102 rows and 46 variables.
#' @author Mathias Harrer, Paula Kuper, Clara Miguel, Pim Cuijpers
#' @usage data("depressionPsyCtr")
"depressionPsyCtr"


#'  The 'inpatients' dataset with 32 clinical trials
#'
#' An example dataset containing data of 32 clinical trials. Its format equals
#' the format of data in the "metapsy" database that is similar to the long format
#' in \emph{R}. The dataset also contains columns with study characteristics that
#' are important for effect size calculation and more.
#'
#' @format A \code{data.frame} with 179 rows and 46 variables:
#' \describe{
#'   \item{study}{\code{character}, The study label containing the author(s) and year of the study.}
#'   \item{condition}{\code{character}, The condition of the groups, either "intervention group" or "control group".}
#'   \item{Cond_spec}{\code{character}, The specific intervention in conditions.}
#'   \item{is.multiarm}{\code{numeric}, The dichotomized indication if a study has multiple arms or not.}
#'   \item{no.arms}{\code{numeric}, The number of arms of a study.}
#'   \item{multiple.arms}{\code{character}, The specification of arms in multiarm studies.}
#'   \item{Outc_type}{\code{character}, The type of outcome.}
#'   \item{primary}{\code{numeric}, The indication if a outcome is primary or not.}
#'   \item{Outc_measure}{\code{character}, The outcome measure used.}
#'   \item{Time}{\code{character}, The dichotomized time of assessment, either "post" or "FU".}
#'   \item{Time_weeks}{\code{character}, The assessment time of FU in weeks.}
#'   \item{year}{\code{numeric}, The year of the study.}
#'   \item{country}{\code{character}, The country of the study.}
#'   \item{...}{}
#' }
#' @author Mathias Harrer, Paula Kuper, Pim Cuijpers
#' @usage data("inpatients")
#' @keywords internal
"inpatients"



#'  The 'psyCtrSubset' dataset
#'
#' An example dataset containing a subset of the 2021 depression trials database.
#' Its format equals the format of data in the "metapsy" database, which is similar to the
#' long format in \emph{R}. The dataset also contains columns with study characteristics that
#' are important for effect size calculation and more.
#'
#' The dataset has been included to showcase the correct formatting necessary to use \code{metapsyTools}.
#'
#' @format A \code{data.frame} with 503 rows and 60 variables.
#' @author Mathias Harrer, Paula Kuper, Pim Cuijpers
#' @usage data("psyCtrSubset")
#' @keywords internal
"psyCtrSubset"



#'  The 'psyCtrSubset' dataset (wide format)
#'
#' An example dataset containing a subset of the 2021 depression trials database.
#' In constrast to, e.g., the \code{psyCtrSubset} dataset, \code{psyCtrSubsetWide}
#' is structured in the "wide(r)" format. The dataset also contains columns with study characteristics that
#' are important for effect size calculation and more.
#'
#' The dataset has been included to showcase the correct formatting necessary to use \code{metapsyTools}.
#'
#' @format A \code{data.frame} with 175 rows and 64 variables.
#' @author Mathias Harrer, Paula Kuper, Clara Miguel, Pim Cuijpers
#' @usage data("psyCtrSubsetWide")
#' @keywords internal
"psyCtrSubsetWide"


#'  The 'database2021Subset' dataset
#'
#' An example dataset containing a subset of the 2021 depression trials database.
#' Its format equals the format of data in the "Metapsy" database.
#'
#' The dataset has been included to showcase the correct formatting necessary to use \code{metapsyTools}.
#'
#' @format A \code{data.frame} with 381 rows and 63 variables.
#' @author Mathias Harrer, Paula Kuper, Pim Cuijpers
#' @usage data("database2021Subset")
#' @keywords internal
"database2021Subset"




