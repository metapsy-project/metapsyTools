% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imputeResponse.R
\name{imputeResponse}
\alias{imputeResponse}
\title{Impute response rates based on continuous outcome data}
\usage{
imputeResponse(m.trt.pre, m.trt.post, sd.trt.post, n.trt, 
               m.ctr.pre, m.ctr.post, sd.ctr.post, n.ctr,
               sd.trt.pre, sd.ctr.pre, 
               rho = 0.8, cutoff = NULL, lower.is.better = TRUE,
               method = c("cutpoint", "rci"))
}
\arguments{
\item{m.trt.pre}{Pre-test mean in the (treatment) arm.}

\item{m.trt.post}{Post-test mean in the (treatment) arm.}

\item{sd.trt.post}{Post-test standard deviation in the (treatment) arm.}

\item{n.trt}{Sample size in the (treatment) arm.}

\item{m.ctr.pre}{\emph{Optional}. Pre-test mean in the control arm.}

\item{m.ctr.post}{\emph{Optional}. Post-test mean in the control arm.}

\item{sd.ctr.post}{\emph{Optional}. Post-test standard deviation in the control arm.}

\item{n.ctr}{\emph{Optional}. Sample size in the control arm.}

\item{sd.trt.pre}{Pre-test standard deviation in the (treatment)
group. Required if \code{method = "rci"}.}

\item{sd.ctr.pre}{\emph{Optional}. Pre-test standard deviation in the control
group. Required if \code{method = "rci"}.}

\item{rho}{Reliability coefficient, used to calculate the reliable change when
\code{method = "rci"}. Set to 0.8 by default.}

\item{cutoff}{User defined cut-off score for response. \code{NULL} by default,
which means that 50\% symptom decrease is used as the response criterion.}

\item{lower.is.better}{Do lower values indicate better outcomes (e.g. less
depression)? \code{TRUE} by default.}

\item{method}{Method to define response. Either \code{"cutpoint"} (default; uses
50\% symptom decrease as criterion, or a user-specified score in \code{cutoff}) or
\code{"rci"}, which means that the reliable change index (Jacobson & Truax, 1991) is used.}
}
\value{
If only values for the treatment arm are specified, \code{imputeResponse}
returns a \code{data.frame} with two columns:
\itemize{
\item \code{trtResponder}: Number of responders
\item \code{nTrt}: Total sample size of the treatment arm
}

If values are also specified for the control group, six additional columns are
returned:
\itemize{
\item \code{ctrResponder}: Number of responders in the control group
\item \code{nCtr}: Total sample size of the treatment arm
\item \code{logRR}: Log-risk ratio of the difference in response rates in both arms
\item \code{seLogRR}: Standard error of the log-risk ratio
\item \code{logOR}: Log-odds ratio of the difference in response rates in both arms
\item \code{seLogOR}: Standard error of the log-odds ratio
}
}
\description{
This function allows to impute response rates based on the post-test
means, standard deviations and sample size using the normal-approximation
method by Furukawa et al. (2005). Response rates can be imputed using the 50\%
symptom decrease threshold, user-defined cut-off values, or the reliable
change index (Jacobson & Truax, 1991).
}
\details{
For more details see the \href{https://tools.metapsy.org/articles/metapsytools}{Get Started} vignette.
}
\examples{
\dontrun{

# Calculate response using 50\% method for one group.
imputeResponse(m.trt.pre = 20, m.trt.post = 11,
               sd.trt.post = 7, n.trt = 100)

# Calculate response for two groups and calculate effect sizes
imputeResponse(m.trt.pre = 20, m.trt.post = 11,
               sd.trt.post = 7, n.trt = 100,
               m.ctr.pre = 20, m.ctr.post = 18,
               sd.ctr.post = 6, n.ctr = 120)

# Calculate using user-defined response threshold
imputeResponse(m.trt.pre = 20, m.trt.post = 11,
               sd.trt.post = 7, n.trt = 100,
               m.ctr.pre = 20, m.ctr.post = 18,
               sd.ctr.post = 6, n.ctr = 120,
               cutoff = 15)

# Assuming higher outcomes are better...
imputeResponse(m.trt.pre = 20, m.trt.post = 11,
               sd.trt.post = 7, n.trt = 100,
               m.ctr.pre = 20, m.ctr.post = 18,
               sd.ctr.post = 6, n.ctr = 120,
               cutoff = 15, lower.is.better = FALSE)

# Using RCI, with an assumed reliability of 0.88
imputeResponse(m.trt.pre = 20, m.trt.post = 11,
               sd.trt.post = 7, n.trt = 100,
               m.ctr.pre = 20, m.ctr.post = 18,
               sd.ctr.post = 6, n.ctr = 120,
               method = "rci", rho = 0.88,
               sd.trt.pre = 8)

# Using a different pre-test SD in the control
imputeResponse(m.trt.pre = 20, m.trt.post = 11,
               sd.trt.post = 7, n.trt = 100,
               m.ctr.pre = 20, m.ctr.post = 18,
               sd.ctr.post = 6, n.ctr = 120,
               method = "rci", rho = 0.88,
               sd.trt.pre = 8, sd.ctr.pre = 9)

# Calculating many results at the same time
set.seed(123)
imputeResponse(m.trt.pre = runif(10, 10, 30), m.trt.post = runif(10, 8, 14),
               sd.trt.post = runif(10, 5, 8), n.trt = rpois(10, 150),
               m.ctr.pre = runif(10, 10, 30), m.ctr.post = runif(10, 13, 20),
               sd.ctr.post = runif(10, 5, 8), n.ctr = rpois(10, 150),
               method = "rci", rho = round(runif(10, 70, 90))/100,
               sd.trt.pre = runif(10, 7, 9))

set.seed(123)
imputeResponse(m.trt.pre = runif(10, 10, 30), m.trt.post = runif(10, 8, 14),
               sd.trt.post = runif(10, 5, 8), n.trt = rpois(10, 150),
               m.ctr.pre = runif(10, 10, 30), m.ctr.post = runif(10, 13, 20),
               sd.ctr.post = runif(10, 5, 8), n.ctr = rpois(10, 150))

set.seed(123)
imputeResponse(m.trt.pre = runif(10, 10, 30), m.trt.post = runif(10, 8, 14),
               sd.trt.post = runif(10, 5, 8), n.trt = rpois(10, 150),
               cutoff = 11:20)
}

}
\references{
Furukawa, T. A., Cipriani, A., Barbui, C., Brambilla, P., & Watanabe, N. (2005).
Imputing response rates from means and standard deviations in meta-analyses.
\emph{International Clinical Psychopharmacology, 20}(1), 49-52.

Jacobson, N. S., & Truax, P. (1991). Clinical significance: a statistical
approach to defining meaningful change in psychotherapy research. \emph{Journal of
Consulting and Clinical Psychology, 59}(1), 12–19.
https://doi.org/10.1037//0022-006x.59.1.12
}
\seealso{
\code{\link{calculateEffectSizes}}
}
\author{
Mathias Harrer \email{mathias.h.harrer@gmail.com}
}
