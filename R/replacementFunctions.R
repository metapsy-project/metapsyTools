#' Replacement functions for "runMetaAnalysis" results objects
#'
#' Once a model has been fitted using `runMetaAnalysis`,
#' replacement functions are defined for each function argument. This 
#' allows to quickly tweak one or more analysis settings, which are implemented
#' once the `rerun` function is called.
#' 
#' @name Replacement functions
#'
#' @param x Object of class `runMetaAnalysis`.
#' @param value Value of one of the arguments in `runMetaAnalysis` or 
#' `correctPublicationBias`
#' @param m (Adapted) object of class `runMetaAnalysis`.
#' @usage 
#' data(x) <- value
#' which.run(x) <- value
#' es.measure(x) <- value
#' es.type(x) <- value
#' es.var(x) <- value
#' se.var(x) <- value
#' es.binary.raw.vars(x) <- value
#' method.tau(x) <- value
#' i2.ci.threelevel(x) <- value
#' nsim.boot(x) <- value
#' hakn(x) <- value
#' study.var(x) <- value
#' arm.var.1(x) <- value
#' arm.var.2(x) <- value
#' measure.var(x) <- value
#' low.rob.filter(x) <- value
#' method.tau.ci(x) <- value
#' which.combine(x) <- value
#' which.combine.var(x) <- value
#' which.outliers(x) <- value
#' which.influence(x) <- value
#' which.rob(x) <- value
#' nnt.cer(x) <- value
#' rho.within.study(x) <- value
#' phi.within.study(x) <- value
#' w1.var(x) <- value
#' w2.var(x) <- value
#' vcov(x) <- value
#' near.pd(x) <- value
#' use.rve(x) <- value
#' html(x) <- value
#' lower.is.better(x) <- value
#' selmodel.steps(x) <- value
#' rerun(m)
#' 
#' @export data<- which.run<- which.run near.pd near.pd<- es.measure es.measure<- es.type es.type<- es.var es.var<- se.var se.var<- es.binary.raw.vars es.binary.raw.vars<- method.tau method.tau<- hakn hakn<- study.var study.var<- arm.var.1 arm.var.1<- arm.var.2 arm.var.2<- measure.var measure.var<- low.rob.filter low.rob.filter<- method.tau.ci method.tau.ci<- which.combine which.combine<- which.combine.var which.combine.var<- which.outliers which.outliers<- which.influence which.influence<- which.rob which.rob<- nnt.cer nnt.cer<- rho.within.study rho.within.study<- phi.within.study phi.within.study<- w1.var w1.var<- w2.var w2.var<- vcov<- use.rve use.rve<- html html<- lower.is.better lower.is.better<- i2.ci.threelevel<- i2.ci.threelevel nsim.boot<- nsim.boot selmodel.steps selmodel.steps<- rerun
#' @aliases data which.run near.pd es.measure es.type es.var se.var es.binary.raw.vars 
#' method.tau hakn study.var arm.var.1 arm.var.2 measure.var low.rob.filter 
#' method.tau.ci round.digits which.combine which.combine.var which.outliers 
#' which.influence which.rob nnt.cer rho.within.study phi.within.study 
#' w1.var w2.var time.var vcov use.rve html lower.is.better selmodel.steps
#' rerun data<- which.run<-
#' es.measure<- es.type<- es.var<- se.var<- es.binary.raw.vars<- 
#' method.tau<- hakn<- study.var<- arm.var.1<- arm.var.2<- 
#' measure.var<- low.rob.filter<- method.tau.ci<- round.digits<- 
#' which.combine<- which.combine.var<- which.outliers<- 
#' which.influence<- which.rob<- nnt.cer<- rho.within.study<- 
#' phi.within.study<- w1.var<- w2.var<- time.var<- vcov<- 
#' use.rve<- html<- lower.is.better<- selmodel.steps<- near.pd near.pd<-
#' i2.ci.threelevel<- i2.ci.threelevel nsim.boot<- nsim.boot
#' @seealso \code{\link{runMetaAnalysis}}, \code{\link{correctPublicationBias}}
#' @examples 
#' \dontrun{
#' data("depressionPsyCtr")
#' depressionPsyCtr %>%
#'   checkDataFormat() %>%
#'   checkConflicts() %>%
#'   calculateEffectSizes() %>% 
#'   filterPoolingData(condition_arm2 %in% 
#'                       c("wl", "other ctr")) -> data
#' 
#' m <- runMetaAnalysis(data, "combined")
#' 
#' # Compare results when other tau^2 estimator is used
#' method.tau(m) <- "DL"
#' rerun(m)
#' }
`data<-` = function(x, value){
  x$data = value
  return(x)
}

`which.run<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$which.run = value
  x$call = as.call(as.call(origArgs))
  return(x)
}
which.run = `which.run<-`

`es.measure<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$es.measure = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
es.measure = `es.measure<-`

`es.type<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$es.type = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
es.type = `es.type<-`

`es.var<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$es.var = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
es.var = `es.var<-`

`se.var<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$se.var = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
se.var = `se.var<-`

`es.binary.raw.vars<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$es.binary.raw.vars = value
  x$call = as.call(as.call(origArgs))
  return(x)
}
es.binary.raw.vars = `es.binary.raw.vars<-`


`method.tau<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$method.tau = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
method.tau = `method.tau<-`

`i2.ci.threelevel<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$i2.ci.threelevel = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
i2.ci.threelevel = `i2.ci.threelevel<-`

`nsim.boot<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$nsim.boot = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
nsim.boot = `nsim.boot<-`

`hakn<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$hakn = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
hakn = `hakn<-`

`study.var<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$study.var = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
study.var = `study.var<-`

`arm.var.1<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$arm.var.1 = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
arm.var.1 = `arm.var.1<-`

`arm.var.2<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$arm.var.2 = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
arm.var.2 = `arm.var.2<-`

`measure.var<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$measure.var = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
measure.var = `measure.var<-`

`low.rob.filter<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$low.rob.filter = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
low.rob.filter = `low.rob.filter<-`

`method.tau.ci<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$method.tau.ci = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
method.tau.ci = `method.tau.ci<-`

`round.digits<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$round.digits = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
round.digits = `round.digits<-`

`which.combine<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$which.combine = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
which.combine = `which.combine<-`

`which.combine.var<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$which.combine.var = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
which.combine.var = `which.combine.var<-`

`which.outliers<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$which.outliers = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
which.outliers = `which.outliers<-`

`which.influence<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$which.influence = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
which.influence = `which.influence<-`

`which.rob<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$which.rob = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
which.rob = `which.rob<-`

`nnt.cer<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$nnt.cer = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
nnt.cer = `nnt.cer<-`

`rho.within.study<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$rho.within.study = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
rho.within.study = `rho.within.study<-`

`phi.within.study<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$phi.within.study = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
phi.within.study = `phi.within.study<-`

`w1.var<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$w1.var = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
w1.var = `w1.var<-`

`w2.var<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$w2.var = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
w2.var = `w2.var<-`

`time.var<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$time.var = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
time.var = `time.var<-`

`vcov<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$vcov = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
vcov = `vcov<-`

`near.pd<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$near.pd = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
near.pd = `near.pd<-`

`use.rve<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$use.rve = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
use.rve = `use.rve<-`

`html<-` = function(x, value){
  origArgs = as.list(x$call)
  origArgs$html = value[1]
  x$call = as.call(as.call(origArgs))
  return(x)
}
html = `html<-`


`lower.is.better<-` = function(x, value){
  if (is.null(x$correctPublicationBias)){
    stop("'correctPublicationBias' must be run first.")
  }
  origArgs = as.list(x$correctPublicationBias$call)
  origArgs$lower.is.better = value[1]
  x$correctPublicationBias$call = as.call(origArgs)
  return(x)
}
lower.is.better = `lower.is.better<-`

`selmodel.steps<-` = function(x, value){
  if (is.null(x$correctPublicationBias)){
    stop("'correctPublicationBias' must be run first.")
  }
  origArgs = as.list(x$correctPublicationBias$call)
  origArgs$selmodel.steps = value
  x$correctPublicationBias$call = as.call(origArgs)
  return(x)
}
selmodel.steps = `selmodel.steps<-`


rerun = function(m){
  d = m$data
  m$call$data = d
  rMA = eval(m$call)
  if (!is.null(m$correctPublicationBias)){
    m$correctPublicationBias$call$model = rMA
    rMA = eval(m$correctPublicationBias$call)
  }
  return(rMA)
}


