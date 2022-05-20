#' Impute response rates based on continuous outcome data
#'
#' This function allows to impute response rates based on the post-test 
#' means, standard deviations and sample size using the normal-approximation 
#' method by Furukawa et al. (2005). Response rates can be imputed using the 50\% 
#' symptom decrease threshold, user-defined cut-off values, or the reliable 
#' change index (Jacobson & Truax, 1991).
#'
#' @usage imputeResponse(m.trt.pre, m.trt.post, sd.trt.post, n.trt, 
#'                m.ctr.pre, m.ctr.post, sd.ctr.post, n.ctr,
#'                sd.trt.pre, sd.ctr.pre, 
#'                rho = 0.8, cutoff = NULL, lower.is.better = TRUE,
#'                method = c("cutpoint", "rci"))
#'
#' @param m.trt.pre Pre-test mean in the (treatment) arm.
#' @param m.trt.post Post-test mean in the (treatment) arm.
#' @param sd.trt.post Post-test standard deviation in the (treatment) arm.
#' @param n.trt Sample size in the (treatment) arm.
#' @param m.ctr.pre \emph{Optional}. Pre-test mean in the control arm.
#' @param m.ctr.post \emph{Optional}. Post-test mean in the control arm.
#' @param sd.ctr.post \emph{Optional}. Post-test standard deviation in the control arm.
#' @param n.ctr \emph{Optional}. Sample size in the control arm.
#' @param sd.trt.pre Pre-test standard deviation in the (treatment)
#' group. Required if \code{method = "rci"}.
#' @param sd.ctr.pre \emph{Optional}. Pre-test standard deviation in the control
#' group. Required if \code{method = "rci"}.
#' @param rho Reliability coefficient, used to calculate the reliable change when
#' \code{method = "rci"}. Set to 0.8 by default.
#' @param cutoff User defined cut-off score for response. \code{NULL} by default,
#' which means that 50\% symptom decrease is used as the response criterion.
#' @param lower.is.better Do lower values indicate better outcomes (e.g. less
#' depression)? \code{TRUE} by default.
#' @param method Method to define response. Either \code{"cutpoint"} (default; uses
#' 50\% symptom decrease as criterion, or a user-specified score in \code{cutoff}) or
#' \code{"rci"}, which means that the reliable change index (Jacobson & Truax, 1991) is used.
#'
#' @return If only values for the treatment arm are specified, \code{imputeResponse}
#' returns a \code{data.frame} with two columns:
#' \itemize{
#'  \item{\code{trtResponder}: }{Number of responders}
#'  \item{\code{nTrt}: }{Total sample size of the treatment arm}
#' }
#' If values are also specified for the control group, six additional columns are
#' returned:
#' \itemize{
#'  \item{\code{ctrResponder}: }{Number of responders in the control group}
#'  \item{\code{nCtr}: }{Total sample size of the treatment arm}
#'  \item{\code{logRR}: }{Log-risk ratio of the difference in response rates in both arms}
#'  \item{\code{seLogRR}: }{Standard error of the log-risk ratio}
#'  \item{\code{logOR}: }{Log-odds ratio of the difference in response rates in both arms}
#'  \item{\code{seLogOR}: }{Standard error of the log-odds ratio}
#' }
#' 
#'
#' @examples
#' \dontrun{
#'
#' # Calculate response using 50% method for one group.
#' imputeResponse(m.trt.pre = 20, m.trt.post = 11,
#'                sd.trt.post = 7, n.trt = 100)
#' 
#' # Calculate response for two groups and calculate effect sizes
#' imputeResponse(m.trt.pre = 20, m.trt.post = 11,
#'                sd.trt.post = 7, n.trt = 100,
#'                m.ctr.pre = 20, m.ctr.post = 18,
#'                sd.ctr.post = 6, n.ctr = 120)
#' 
#' # Calculate using user-defined response threshold
#' imputeResponse(m.trt.pre = 20, m.trt.post = 11,
#'                sd.trt.post = 7, n.trt = 100,
#'                m.ctr.pre = 20, m.ctr.post = 18,
#'                sd.ctr.post = 6, n.ctr = 120,
#'                cutoff = 15)
#' 
#' # Assuming higher outcomes are better...
#' imputeResponse(m.trt.pre = 20, m.trt.post = 11,
#'                sd.trt.post = 7, n.trt = 100,
#'                m.ctr.pre = 20, m.ctr.post = 18,
#'                sd.ctr.post = 6, n.ctr = 120,
#'                cutoff = 15, lower.is.better = FALSE)
#' 
#' # Using RCI, with an assumed reliability of 0.88
#' imputeResponse(m.trt.pre = 20, m.trt.post = 11,
#'                sd.trt.post = 7, n.trt = 100,
#'                m.ctr.pre = 20, m.ctr.post = 18,
#'                sd.ctr.post = 6, n.ctr = 120,
#'                method = "rci", rho = 0.88,
#'                sd.trt.pre = 8)
#' 
#' # Using a different pre-test SD in the control
#' imputeResponse(m.trt.pre = 20, m.trt.post = 11,
#'                sd.trt.post = 7, n.trt = 100,
#'                m.ctr.pre = 20, m.ctr.post = 18,
#'                sd.ctr.post = 6, n.ctr = 120,
#'                method = "rci", rho = 0.88,
#'                sd.trt.pre = 8, sd.ctr.pre = 9)
#' 
#' # Calculating many results at the same time
#' set.seed(123)
#' imputeResponse(m.trt.pre = runif(10, 10, 30), m.trt.post = runif(10, 8, 14),
#'                sd.trt.post = runif(10, 5, 8), n.trt = rpois(10, 150),
#'                m.ctr.pre = runif(10, 10, 30), m.ctr.post = runif(10, 13, 20),
#'                sd.ctr.post = runif(10, 5, 8), n.ctr = rpois(10, 150),
#'                method = "rci", rho = round(runif(10, 70, 90))/100,
#'                sd.trt.pre = runif(10, 7, 9))
#' 
#' set.seed(123)
#' imputeResponse(m.trt.pre = runif(10, 10, 30), m.trt.post = runif(10, 8, 14),
#'                sd.trt.post = runif(10, 5, 8), n.trt = rpois(10, 150),
#'                m.ctr.pre = runif(10, 10, 30), m.ctr.post = runif(10, 13, 20),
#'                sd.ctr.post = runif(10, 5, 8), n.ctr = rpois(10, 150))
#' 
#' set.seed(123)
#' imputeResponse(m.trt.pre = runif(10, 10, 30), m.trt.post = runif(10, 8, 14),
#'                sd.trt.post = runif(10, 5, 8), n.trt = rpois(10, 150),
#'                cutoff = 11:20)
#' }
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com}
#'
#' @seealso \code{\link{calculateEffectSizes}}
#' 
#' @references 
#' Furukawa, T. A., Cipriani, A., Barbui, C., Brambilla, P., & Watanabe, N. (2005). 
#' Imputing response rates from means and standard deviations in meta-analyses. 
#' \emph{International Clinical Psychopharmacology, 20}(1), 49-52.
#' 
#' Jacobson, N. S., & Truax, P. (1991). Clinical significance: a statistical 
#' approach to defining meaningful change in psychotherapy research. \emph{Journal of 
#' Consulting and Clinical Psychology, 59}(1), 12–19. 
#' https://doi.org/10.1037//0022-006x.59.1.12
#' 
#'
#' @details For more details see the help vignette: \code{vignette("metapsyTools")}.
#'
#' @importFrom stats dffits model.matrix rnorm rstudent pnorm qnorm
#' @export imputeResponse
#' 


imputeResponse = function(m.trt.pre, m.trt.post, sd.trt.post,
                          n.trt, m.ctr.pre, m.ctr.post, sd.ctr.post, n.ctr,
                          sd.trt.pre, sd.ctr.pre, rho = 0.8,
                          cutoff = NULL, lower.is.better = TRUE,
                          method = c("cutpoint", "rci")){
  

  if (method[1] == "cutpoint"){
    if (is.null(cutoff)){
      if (missing(m.trt.pre) | missing(m.trt.post) | missing(sd.trt.post) |
          missing(n.trt)){
        warning(paste("One of 'm.trt.pre', 'm.trt.post', 'sd.trt.post' or 'n.trt'",
                      "is missing. Returning NA."))
        m.trt.pre = NA
        m.trt.post = NA
        sd.trt.post = NA
        n.trt = NA
        n.resp = NA
      }

      if (lower.is.better[1] == TRUE){
        round(pnorm(((m.trt.pre/2)-m.trt.post)/
                      sd.trt.post)*n.trt) -> n.resp
        message("50% symptom decrease assumed as clinically relevant threshold.")
      } else {
        round(pnorm(((m.trt.pre*2)-m.trt.post)/
                      sd.trt.post)*n.trt) -> n.resp
        message("50% increase in values assumed as clinically relevant threshold.")
      }
    } else {
      if (missing(m.trt.post) | missing(sd.trt.post) | missing(n.trt)){
        warning(paste("One of 'm.trt.post', 'sd.trt.post' or 'n.trt'",
                      "is missing. Returning NA."))
        m.trt.post = NA
        sd.trt.post = NA
        n.trt = NA
        n.resp = NA
      }
      if (lower.is.better[1] == TRUE){
        round(pnorm((cutoff-m.trt.post)/
                      sd.trt.post)*n.trt) -> n.resp
        if (length(cutoff) == 1){
          message(paste("A value of", cutoff,
                        "was assumed as the clinically relevant threshold."))}
      } else {
        round(pnorm((m.trt.post-cutoff)/
                      sd.trt.post)*n.trt) -> n.resp
        if (length(cutoff) == 1){
          message(paste("A value of", cutoff,
                        "was assumed as the clinically relevant threshold",
                        "(higher values indicate better outcomes)."))}
      }
    }
  }

  if (method[1] == "rci"){
    if (missing(m.trt.post) | missing(sd.trt.post) | missing(n.trt)){
      warning(paste("One of 'm.trt.post', 'sd.trt.post' or 'n.trt'",
                    "is missing. Returning NA."))
      m.trt.post = NA
      sd.trt.post = NA
      n.trt = NA
      n.resp = NA
    }
    if (missing(sd.trt.pre)){
      warning("Variable 'sd.trt.pre' is missing.",
              " Calculation using RCI as threshold will return NA.")
      sd.trt.pre = NA
    }
    if (sum(rho < 0) > 1 | sum(rho > 1) > 1){
      error("Value of 'reliability' must be between 0 and 1.")
    }
    sdiff = sqrt(2*(sd.trt.pre*sqrt(1-rho))^2)
    rcrit = qnorm(0.975)*sdiff

    if (lower.is.better[1] == TRUE){
      rcrit.m = m.trt.pre - rcrit
      round(pnorm((rcrit.m-m.trt.post)/
                    sd.trt.post)*n.trt) -> n.resp
      if (length(rcrit) == 1){
        message(paste0("A reliable change value of ", round(rcrit, 3),
                       " was assumed as the clinically relevant threshold,",
                       " using a reliability of ryy=", rho, "."))}
    } else {
      rcrit.m = m.trt.pre + rcrit
      round(pnorm((m.trt.post-rcrit.m)/
                    sd.trt.post)*n.trt) -> n.resp
      if (length(rcrit) == 1){
        message(paste0("A reliable change value of ", round(rcrit, 3),
                       " was assumed as the clinically relevant threshold,",
                       " using a reliability of ryy=", rho,
                       " (higher values indicate better outcomes)."))}
    }}



  # Effect Size Calculation
  if (method[1] == "cutpoint"){
    if (missing(m.ctr.pre) | missing(m.ctr.post) |
        missing(sd.ctr.post) | missing(n.ctr)){
      m.ctr.pre = NA
      m.ctr.post = NA
      sd.ctr.post = NA
      n.ctr = NA
      esCalc = FALSE
    } else {
      if (method[1] == "cutpoint"){
        if (is.null(cutoff)){

          if (lower.is.better[1] == TRUE){
            round(pnorm(((m.ctr.pre/2)-m.ctr.post)/
                          sd.ctr.post)*n.ctr) -> n.resp.ctr
          } else {
            round(pnorm(((m.ctr.pre*2)-m.ctr.post)/
                          sd.ctr.post)*n.ctr) -> n.resp.ctr
          }
        } else {

          if (lower.is.better[1] == TRUE){
            round(pnorm((cutoff-m.ctr.post)/
                          sd.ctr.post)*n.ctr) -> n.resp.ctr
          } else {
            round(pnorm((m.ctr.post-cutoff)/
                          sd.ctr.post)*n.ctr) -> n.resp.ctr
          }
        }
      }
      esCalc = TRUE
    }}


  if (method[1] == "rci"){

    if (missing(m.ctr.pre) | missing(m.ctr.post) |
        missing(sd.ctr.post) | missing(n.ctr)){
      m.ctr.pre = NA
      m.ctr.post = NA
      sd.ctr.post = NA
      n.ctr = NA
      esCalc = FALSE
    }

    if (missing(sd.ctr.pre)){
      message("Assuming same pre-test SD in control and treatment group.")
      sd.ctr.pre = sd.trt.pre
    }

    sdiff = sqrt(2*(sd.ctr.pre*sqrt(1-rho))^2)
    rcrit = qnorm(0.975)*sdiff

    if (lower.is.better[1] == TRUE){
      rcrit.m = m.ctr.pre - rcrit
      round(pnorm((rcrit.m-m.ctr.post)/
                    sd.ctr.post)*n.ctr) -> n.resp.ctr
    } else {
      rcrit.m = m.ctr.pre + rcrit
      round(pnorm((m.ctr.post-rcrit.m)/
                    sd.ctr.post)*n.ctr) -> n.resp.ctr
    }
    esCalc = TRUE}

  if (esCalc){
    log((n.resp/n.trt)/(n.resp.ctr/n.ctr)) -> logRR
    sqrt((1/n.resp) + (1/n.resp.ctr) - (1/n.trt) - (1/n.ctr)) -> seLogRR

    log((n.resp/(n.trt-n.resp))/(n.resp.ctr/(n.ctr-n.resp.ctr))) -> logOR
    sqrt(n.resp^-1 + (n.trt-n.resp)^-1 +
           n.resp.ctr^-1 + (n.ctr-n.resp.ctr)^-1) -> seLogOR

    return(data.frame(trtResponder = n.resp, nTrt = n.trt,
                      ctrResponder = n.resp.ctr, nCtr = n.ctr,
                      logRR = logRR, seLogRR = seLogRR,
                      logOR = logOR, seLogOR = seLogOR))

  } else {
    return(data.frame(trtResponder = n.resp, nTrt = n.trt))
  }
}
