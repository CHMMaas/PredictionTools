#' @title fun.event
#' @description This function translates the linear predictor and baseline hazard to predictions.
#'
#' @param lp vector of linear predictor
#' @param h0 baseline hazard value at the time point the prediction needs to be made
#'
#' @return The predicted probability at the time point of the provided baseline hazard
#'
#' p
#'
#' Predicted probability at the time point of the provided baseline hazard
#'
#' @export
#'
#' @examples
#' library(PredictionTools)
#' library(survival)
#' library(rms)
#' set.seed(1)
#'
#' # load data
#' data(pbc, package="survival")
#'
#' # set survival time
#' S <- survival::Surv(time=pbc$time/365, event=as.numeric(pbc$status!=0))
#'
#' # limit survival to time horizon
#' horizon <- 5
#' S.5 <- S
#' S.5[S[,1]>horizon, 1] <- horizon # if time is more than horizon, set to horizon
#' S.5[S[,1]>horizon, 2] <- 0       # if time is more than horizon, censor
#'
#' # fit model
#' model <- cph(S.5 ~ sex + age + bili, data=pbc, x=TRUE, y=TRUE, se.fit=TRUE)
#'
#' # obtain linear predictor
#' lp <- predict(model, type="lp")
#'
#' # obtain baseline hazard
#' f.basehaz <- survival::basehaz(model)
#'
#' # obtain baseline hazard at horizon
#' h0 <- f.basehaz$hazard[f.basehaz$time==max(f.basehaz$time[f.basehaz$time<=horizon])]
#'
#' # make predictions at time horizon
#' p <- PredictionTools::fun.event(lp=lp, h0=h0)
fun.event<-function(lp, h0){
  stopifnot("est must be numeric" = is.numeric(lp))
  stopifnot("se must be numeric" = is.numeric(h0))
  stopifnot("est must be a vector" = is.vector(lp))

  h <- h0*exp(lp)
  p <- 1-exp(-h)

  return(p)
}
