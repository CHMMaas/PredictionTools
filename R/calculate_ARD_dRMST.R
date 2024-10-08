#' @title Calculate ARD and dRMST
#' @description This function calculates the absolute risk difference (ARD) and difference in restricted mean survival time (dRMST)
#'
#' @param S the outcome, a Surv() object
#' @param W treatment indicator, 1 for treatment and 0 for control
#' @param horizon the prediction horizon
#'
#' @return The output of the calculate.ARD.dRMST function is a "list" with the following components.
#'
#' LRT.p
#'
#' p-value of the log-rank test
#'
#'
#' ARD
#'
#' absolute risk difference
#'
#'
#' ARD_se
#'
#' standard error of absolute risk difference
#'
#'
#' dRMST
#'
#' difference in restricted mean survival time
#'
#'
#' dRMST_se
#'
#' standard error of differnece in restricted mean survival time
#'
#'
#' @export
#'
#' @examples
#' library(PredictionTools)
#' library(survival)
#'
#' # load data
#' data(pbc, package="survival")
#'
#' # set survival time
#' S <- survival::Surv(time=pbc$time/365, event=as.numeric(pbc$status!=0))
#'
#' # treatment indicator
#' W <- pbc$trt
#'
#' # calculate ARD and dRMST
#' calculate.ARD.dRMST(S=S, W=W, horizon=5)
calculate.ARD.dRMST <- function(S=NULL, W=NULL, horizon=0){
  # KM survival curve
  KM <- survival::survfit(S ~ W, data=data.frame(S, W))

  # log-rank test
  LRT.p <- survival::survdiff(S ~ W, data=data.frame(S, W))$pvalue

  # ARD
  KM.summary <- summary(KM, times=horizon, extend=TRUE)
  ARD <- diff(KM.summary$surv)
  ARD_se <- sqrt(sum((KM.summary$std.err)^2))

  # dRMST
  RMST.summary <- summary(KM, rmean=horizon, extend=TRUE)$table
  dRMST <- as.numeric(diff(RMST.summary[, "rmean"]))
  dRMST_se <- sqrt(sum((RMST.summary[, "se(rmean)"])^2))

  return(list(LRT.p=LRT.p,
              ARD=ARD, ARD_se=ARD_se,
              dRMST=dRMST, dRMST_se=dRMST_se))
}
