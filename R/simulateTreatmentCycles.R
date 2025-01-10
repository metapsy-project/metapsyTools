#' Simulate the number of treatment cycles and "excess treatments"
#'
#' @description This function allows you to simulate the total number of treatment cycles
#' needed to achieve response in 99% of all patients in a population, and the
#' number of "excess" treatments due to non-response.
#' \loadmathjax
#'
#' @usage simulateTreatmentCycles(response.rate, decay, max.cycles = 1000)
#'
#' @param response.rate The response rate of the treatment(s). Can be either a single numeric value
#' between 0.01 and 0.99, or a vector of single numeric values, indicating the response rate after the
#' first, second, third, etc. treatment cycle.
#' @param decay The proportion (percentage) by which the response rate decreases after each unsuccessful treatment
#' cycle. Must be a numeric value between -0.99 and 0.99. Negative values indicate an _increase_ in response
#' rates after each unsuccessful treatment cycle. A value of 0 indicates stable response rates in each cycle.
#' If `response.rate` includes multiple values, the decay is applied after the last user-provided response rate.
#' @param max.cycles By default, the simulation is terminated if 99% cumulative response in the population has
#' not been reached after 1000 cycles. This can be set to a lower number to model more realistic scenarios 
#' (e.g., no more treatment attempts after 10 cycles, so that `max.cycles=10`).
#'
#' @return Returns an object of class \code{"simulateTreatmentCycles"}. This object includes, among other things,
#' a \code{data.frame} with the name \code{data}, in which all simulation results are stored.
#' 
#' Other objects are the total number of treatments provided (`total.treatments`), the total number of
#' "excess" treatments (`excess.treatments`), and the average number of treatment cycles per person 
#' (`avg.no.treatments`), assuming a population of 100 patients. A `plot` S3 method is also defined for the
#' `simulateTreatmentCycles` object.
#' 
#' @details
#' In the “simple” scenario, the repeated treatment cycles can be modeled as a Markov chain with absorbing state 
#' \mjeqn{S_0}{S_0} (patient responds), as well as transition probabilities 
#' \mjeqn{p_k=P(S_0│S_k)}{p_k=P(S_0│S_k)} and \mjeqn{1-p_k=P(S_{k+1}|S_k)}{1-p_k=P(S_{k+1}|S_k)}. 
#' Here, \mjeqn{k}{k} is the current treatment cycle, and \mjeqn{p_k}{p_k} is the probability of responding to the \mjeqn{k}{k}th treatment, 
#' which is held constant across all treatment cycles. In this scenario, the formula to calculate the cumulative 
#' response \mjeqn{C_n}{C_n} after \mjeqn{n}{n} treatment cycles reduces to:
#' 
#' \mjtdeqn{C_n=1-(1-p)^n}{C_n=1-(1-p)^n}{C_n=1-(1-p)^n}
#' 
#' The number of “excess” treatments \mjeqn{N^E}{N^E} can be obtained using this formula (assuming that 100 patients are treated initially):
#' 
#' \mjtdeqn{N^E\approx\frac{100×(1-(1-p)^n)}{p}-100}{N^E\approx\frac{100×(1-(1-p)^n)}{p}-100}{N^E\approx\frac{100×(1-(1-p)^n)}{p}-100}
#' 
#' If we additionally consider “decay” of the treatment response (e.g., response rates decrease by 10% with each additional treatment attempt), we obtain the following formula for the cumulative response:
#' 
#' \mjtdeqn{C_n = \sum_{i=1}^{n} \left( \prod_{j=1}^{i-1} \left( 1 - p(1-d)^{j-1} \right) \right) p(1-d)^{i-1}}{C_n = \sum_{i=1}^{n} \left( \prod_{j=1}^{i-1} \left( 1 - p(1-d)^{j-1} \right) \right) p(1-d)^{i-1}}{C_n = \sum_{i=1}^{n} \left( \prod_{j=1}^{i-1} \left( 1 - p(1-d)^{j-1} \right) \right) p(1-d)^{i-1}}
#' 
#' where \mjeqn{C_n}{C_n} is the cumulative treatment response, \mjeqn{p}{p} is the response rate of treatments, \mjeqn{d}{d} encodes the proportional “decay” with each treatment cycle, and \mjeqn{n}{n} is the total number of treatment cycles. 
#' 
#' @import mathjaxr
#'
#' @examples
#' \dontrun{
#' # "Simple" scenario: 50% response rate at each cycle
#' res <- simulateTreatmentCycles(response.rate = 0.5, decay = 0)
#' res; plot(res)
#' 
#' # Define manual response rates. Since decay=0, rates stay at 0.4 afterwards
#' res <- simulateTreatmentCycles(response.rate = c(0.57, 0.4), decay = 0)
#' res; plot(res)
#' 
#' # Define manual response rates up to the 4th cycles. Afterwards, response decays by 10%
#' res <- simulateTreatmentCycles(response.rate = c(0.56, 0.3, 0.28, 0.4), decay = .1)
#' res; plot(res)
#' 
#' # Run the simulation for a fixed number of cycles. Here: 10
#' res <- simulateTreatmentCycles(response.rate = c(0.56, 0.3, 0.28, 0.1), decay = .1, max.cycles = 10)
#' res; plot(res)
#' 
#' # Go crazy
#' res <- simulateTreatmentCycles(response.rate = c(.2, .8, .2, .8, .2, .8, .2), decay = -.2)
#' res; plot(res)
#' }
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com},
#' Pim Cuijpers \email{p.cuijpers@@vu.nl}, Toshi Furukawa
#'
#' @details For more details on the metapsyTools package, see the [Get Started](https://tools.metapsy.org/articles/metapsytools) vignette.
#'
#' @export simulateTreatmentCycles
#' @keywords internal


simulateTreatmentCycles = function(response.rate, decay, max.cycles = 1000) {
  
  if (missing(response.rate)) {
    stop("'response.rate' must be provided.", call. = FALSE)
  }
  if (missing(decay)) {
    stop("'decay' must be provided.")
  }
  if (!is.numeric(response.rate[1]) | response.rate[1] < 0.01 | response.rate[1] > 0.99) {
    stop("'response.rate' must be numeric and between 0.01 and 0.99.", call. = FALSE)
  }
  if (!is.numeric(decay[1]) | decay[1] < -0.99 | decay[1] > 0.99) {
    stop("'decay' must be numeric and between -0.99 and 0.99.", call. = FALSE)
  }
  if (!is.numeric(max.cycles[1]) | max.cycles[1] < 1) {
    stop("'max.cycles' must be numeric and >= 1.", call. = FALSE)
  }
  
  if (length(response.rate)>1) {
    response.rate -> R
    manual = TRUE
  } else {
    response.rate[1] -> R
    manual = FALSE
  }
  decay[1] -> p
  max.cycles[1] -> max.cycles
  
  cum.resp <- function(R, p, n) {
    cum.resp <- 0  
    for (i in 1:n) {
      resp.i <- R * (1 - p)^(i - 1)
      if (resp.i > 1) {resp.i <- 1}
      prob.non.resp <- 1  
      if (i > 1) {
        for (j in 1:(i - 1)) {
          prob.non.resp <- prob.non.resp * (1 - R * (1 - p)^(j - 1))
        }
      }
      cum.resp <- cum.resp + (prob.non.resp * resp.i)
    }
    return(c(cum.resp = cum.resp, 
             prob.non.resp = prob.non.resp, 
             prob.resp = resp.i))
  }
  
  cum.resp.manual <- function(R) {
    cum.resp <- 0
    prob.non.resp <- 1; l = list()
    for (i in 1:length(R)) {
      resp.i <- prob.non.resp * R[i]
      cum.resp = cum.resp + resp.i
      c(step = i,
        prop.treated = prob.non.resp, prob.response = R[i],  
        cum.response = cum.resp) -> l[[i]]
      prob.non.resp <- prob.non.resp * (1-R[i])
    }
    return(as.data.frame(do.call(rbind,l)))
  }
  
  if (manual) {
    cum.resp.manual(R) -> res.manual
    within(res.manual, {treated.per.100 = round(prop.treated*100)}) -> res.manual
  }
  
  response = 0; n = 1; l = list()
  while(response < .99 & n < (max.cycles+1)) {
    tmp = cum.resp(R[length(R)], p, n)
    l[[n]] = cbind(step = n, 
                   prop.treated = tmp[["prob.non.resp"]],
                   prob.response = tmp[["prob.resp"]],
                   cum.response = tmp[["cum.resp"]])
    n = n+1; response <- tmp[["cum.resp"]]
    if (n==100) {
      message("Large number of simulations needed. Computations may take several seconds...")
    }
  }
  as.data.frame(do.call(rbind,l)) -> res
  within(res, {treated.per.100 = round(prop.treated*100)}) -> res
  warn = FALSE; capped = FALSE
  if (nrow(res)==1000) {warn = TRUE}
  if (max.cycles != 1000) {capped = TRUE}
  
  if (manual) {
    within(res, {
      prop.treated = prop.treated*res.manual$prop.treated[nrow(res.manual)]
      cum.response = cum.response*res.manual$prop.treated[nrow(res.manual)] + 
        res.manual$cum.response[nrow(res.manual)-1]
      treated.per.100 = round(treated.per.100*res.manual$prop.treated[nrow(res.manual)])
    }) -> res
    rbind(res.manual, res[-1,]) -> res; res$step = 1:nrow(res)
    which.max(res$cum.response > 0.99) -> prob.filter
    ifelse(prob.filter == 1, nrow(res), prob.filter) -> prob.filter
    res[1:prob.filter, ] -> res
    rownames(res) = NULL
    if (nrow(res) > max.cycles) {res[1:max.cycles,] -> res}
    if (nrow(res)==1000) {warn = TRUE} else {warn = FALSE}
  }
  
  total.steps <- max(res$step)
  total.treatments <- sum(res$treated.per.100)
  excess.treatments <- sum(res$treated.per.100)-res$treated.per.100[1]
  avg.no.treatments <- with(res, sum((step * treated.per.100) / sum(treated.per.100)))
  names(res) = c("cycle", "prop.treated", "prob.response", "cum.response", "treated.per.100")
  
  reslist <- list(
    total.cycles = total.steps,
    total.treatments = total.treatments,
    excess.treatments = excess.treatments,
    avg.no.treatments = avg.no.treatments,
    data = res,
    warn = warn,
    capped = capped,
    max.cycles = max.cycles
  )
  
  class(reslist) = c("simulateTreatmentCycles", "list")
  
  return(reslist)
}



#' Print method for objects of class 'simulateTreatmentCycles'
#'
#' Print S3 method for objects of class \code{simulateTreatmentCycles}.
#'
#' @param x An object of class \code{simulateTreatmentCycles}.
#' @param ... Additional arguments (not used).
#'
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com}, Pim Cuijpers 
#' \email{p.cuijpers@@vu.nl}
#'
#'
#' @export
#' @method print simulateTreatmentCycles

print.simulateTreatmentCycles = function(x, ...) {
  if (x$capped[1] & (x$total.cycles == x$max.cycles)) {
    cap = "(capped)"} else { cap = "" }
  if (x$warn) {cap = "(stopped)"}
  cat("Sequential treatment until 100/100 people respond:", fill = TRUE)
  cat("-----------------------------------------------------", fill = TRUE)
  cat("Total number of treatment cycles needed:", x$total.cycles, cap, fill = TRUE)
  cat("Total number of treatments to be provided:", x$total.treatments, fill = TRUE)
  cat("Excess number of treatments to be provided:", x$excess.treatments, fill = TRUE)
  cat("Average number of treatments until response is reached:", round(x$avg.no.treatments,3), fill = TRUE)
  cat("\n")
  if (x$total.cycles <= 20) {
    print(x$data, row.names=FALSE)
  } else {
    print(head(x$data, 5), row.names=FALSE)
    cat("...", fill = TRUE)
    print(tail(x$data, 5), row.names=FALSE)
  }
  if (x$warn) warning("Simulations were stopped after 1000 treatment cycles. Results should be interpreted with caution.", call. = FALSE)
  if (x$capped[1] & (x$total.cycles == x$max.cycles)) warning("Simulations were capped after ", x$max.cycles," treatment cycles. Results should be interpreted with caution.", call. = FALSE)
}



#' Plot method for objects of class 'simulateTreatmentCycles'
#'
#' Plot S3 method for objects of class \code{simulateTreatmentCycles}.
#'
#' @param x An object of class \code{simulateTreatmentCycles}.
#' @param ... Additional arguments (not used).
#' @author Mathias Harrer \email{mathias.h.harrer@@gmail.com}, Pim Cuijpers \email{p.cuijpers@@vu.nl}
#' @importFrom graphics axTicks legend matplot
#' @export
#' @method plot simulateTreatmentCycles

plot.simulateTreatmentCycles = function(x, ...) {
  x$data -> res
  Y = res[c("prop.treated", "cum.response", "prob.response")]
  cols = c("black", "dodgerblue", "darkgreen")
  matplot(res$cycle, Y, type = "l", lty = 1, col = cols, 
          xlab = "Treatment Cycle", ylab = NA, yaxt = "n",
          lwd = 2)
  y.ticks <- axTicks(2)
  axis(2, at = y.ticks, labels = paste0(round(y.ticks * 100), "%"))
  legend("right", legend = c("Needing Treatment", "Responded", "Response Probability"), 
         col = cols, lty = 1, cex = 0.6, lwd = 2)
}

