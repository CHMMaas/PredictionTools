#' @title val.surv.mi
#' @description This function calculates intercept, slope, and C-index for survival risk predictions of multiple imputed data set(s).
#'
#' @importFrom rms Predict
#' @importFrom rms cph
#' @importFrom rms rcs
#' @importFrom survival coxph
#' @importFrom survival survfit
#' @importFrom survival concordance
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
#' @importFrom graphics abline
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics points
#' @importFrom graphics polygon
#' @importFrom graphics segments
#' @importFrom mice mice
#' @importFrom mice complete
#' @importFrom stats qnorm
#' @importFrom stats vcov
#' @importFrom utils tail
#' @importFrom timeROC timeROC
#'
#' @param p Matrix with predicted probabilities for imputation i in columns (complete case analysis: one column)
#' @param y Time to event outcome as Surv object (time,status), either unrestricted or restricted follow-up
#' @param g Number of risk groups; default=5
#' @param time Time point at which to evaluate the predicted probabilities, default=NULL (not entered), the maximum time point will be taken. Please note that AUC doesn't compute at maximum follow-up time, you can use show.metrics to omit these results from the plot.
#' @param main Plot label, default=""
#' @param lim limit, default=NULL
#' @param dist distribution, default=TRUE
#' @param smoothed.curve plot smoothed calibration curve with 95 percent confidence interval, default=FALSE
#' @param CI.metrics plot confidence intervals of calibration intercept, calibration slope, and Harrell's C-index, Uno's C-index, and the area under the time-dependent ROC curve (AUC), default=FALSE.
#' @param show.metrics TRUE/FALSE vector of length 6 indicating if plot should show (1) sample size, (2) calibration intercept, (3) calibration slope, (4) Harrell's C-index possibly corrected with optimism specified in optimism.C, (5) Uno's C-index, (6) the area under time-dependent ROC curve (AUC) as defined by Blanche et al., default=rep(TRUE, 6)
#' @param optimism.C optimism-correction for Harrel's C-index in plot, default=0
#' @param n.sim number of simulations for computing simultaneous confidence bands for the time-dependent AUC, which heavily influence computation time, default=2000
#'
#' @return The output of the val_surv_mi function is a "list" with the following components.
#'
#' main
#'
#' Main title of plots.
#'
#'
#' n
#'
#' number of observations
#'
#'
#' quants
#'
#' quantiles
#'
#'
#' p.mi
#'
#' predicted survival for each imputed data set.
#'
#'
#' obs.mi
#'
#' observed survival for each imputed data set.
#'
#'
#' obs.mi.lower
#'
#' lower bound of 95\% confidence interval of observed survival.
#'
#'
#' obs.mi.upper
#'
#' upper bound of 95\% confidence interval of observed survival.
#'
#'
#' int
#'
#' intercept
#'
#'
#' int.lower
#'
#' lower bound of 95\% confidence interval of intercept.
#'
#'
#' int.upper
#'
#' upper bound of 95\% confidence interval of intercept.
#'
#'
#' slope
#'
#' slope estimate
#'
#'
#' slope.lower
#'
#' lower bound of 95\% confidence interval of slope.
#'
#'
#' slope.upper
#'
#' upper bound of 95\% confidence interval of slope.
#'
#'
#' HarrellC
#'
#' Harrell's C-index, uncorrected for optimism.
#'
#'
#' HarrellC.se
#'
#' standard error of Harrell's C-index
#'
#'
#' HarrellC.lower
#'
#' lower bound of 95\% confidence interval of Harrell's C-index
#'
#' HarrellC.upper
#'
#' upper bound of 95\% confidence interval of Harrell's C-index.
#'
#'
#' UnoC
#'
#' Uno's C-index for time-to-event outcomes, uncorrected for optimism
#'
#'
#' UnoC.se
#'
#' standard error of Uno's C-index, based one a single imputation
#'
#'
#' UnoC.lower
#'
#' lower bound of 95\% confidence interval of Uno's C-index
#'
#'
#' UnoC.lower
#'
#' lower bound of 95\% confidence interval of Uno's C-index
#'
#'
#' AUC
#'
#' time-dependent AUC, calculated using inverse propensity score weighting by the timeROC package, uncorrected for optimism
#'
#'
#' AUC.se
#'
#' standard error of time-dependent AUC, based one a single imputation
#'
#'
#' AUC.lower
#'
#' lower bound of 95\% confidence interval of time-dependent AUC
#'
#'
#' AUC.upper
#'
#' upper bound of 95\% confidence interval of time-dependent AUC
#'
#'
#' @export
#'
#' @examples
#' library(PredictionTools)
#' library(survival)
#' library(rms)
#' library(mice)
#' set.seed(1)
#'
#' # load data
#' data(pbc, package="survival")
#'
#' # Multiple imputation
#' m <- 5
#' candidate.predictors <- c("sex", "age", "bili", "ascites", "hepato")
#' imp <- mice::mice(pbc[, candidate.predictors], m=m, print=FALSE)
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
#' p <- c()
#' for (i in 1:m){
#'  # fit model
#'   model.i <- rms::cph(S.5 ~ sex + age + bili + ascites + hepato,
#'                  data=mice::complete(imp, i),
#'                  x=TRUE, y=TRUE, se.fit=TRUE)
#'
#'   # obtain linear predictor
#'   lp.i <- predict(model.i, type="lp")
#'
#'   # obtain baseline hazard
#'   f.basehaz.i <- survival::basehaz(model.i)
#'
#'   # obtain baseline hazard at horizon
#'   h0.i <- f.basehaz.i$hazard[f.basehaz.i$time==max(f.basehaz.i$time[f.basehaz.i$time<=horizon])]
#'
#'   # make predictions at time horizon
#'   p.i <- PredictionTools::fun.event(lp=lp.i, h0=h0.i)
#'
#'   # add to p
#'   p <- cbind(p, p.i)
#' }
#' # internally validate predictions
#' g <- 4
#' main <- paste("Calibration plot for predictions at time",  horizon)
#' show.metrics <- rep(TRUE, 6)
#' CI.metrics <- TRUE
#' PredictionTools::val.surv.mi(p=as.matrix(p), y=S,
#'                              g=g, time=horizon,
#'                              main=main, lim=c(0, 1),
#'                              dist=TRUE, smoothed.curve=TRUE,
#'                              CI.metrics=CI.metrics,
#'                              show.metrics=show.metrics,
#'                              n.sim=100)
val.surv.mi<-function(p, y, g=5, time=NULL,
                      main="", lim=c(0,1), dist=TRUE, smoothed.curve=FALSE,
                      CI.metrics=FALSE, show.metrics=rep(TRUE, 6),
                      optimism.C=0, n.sim=2000){
  stopifnot("p must be numeric" = is.numeric(p))
  stopifnot("y must be numeric" = is.numeric(y))
  stopifnot("g must be numeric" = is.numeric(g))
  stopifnot("optimism.C must be numeric" = is.numeric(optimism.C))

  stopifnot("dist must be a boolean (TRUE or FALSE)" = isTRUE(dist)|isFALSE(dist))
  stopifnot("smoothed.curve must be a boolean (TRUE or FALSE)" = isTRUE(smoothed.curve)|isFALSE(smoothed.curve))
  stopifnot("CI.metrics must be a boolean (TRUE or FALSE)" = isTRUE(CI.metrics)|isFALSE(CI.metrics))

  stopifnot("p must be a vector or matrix" = is.vector(p) | is.matrix(p))
  stopifnot("y must be a matrix created by survival::Surv()" = is.matrix(y))

  stopifnot("p and y must be the same length" = nrow(p)==nrow(y) | length(p)==nrow(y))

  y.orig <- y
  if (!is.null(time)){
    y[y[,1]>time,2]<-0
    y[y[,1]>time,1]<-time
  } else{
    time <- max(y[, 1])
  }

  lp<-log(-log(1-p))

  n<-length(y)
  m.imp.val<-ncol(lp)

  int<-rep(0,m.imp.val)
  int.se<-rep(0,m.imp.val)
  slope<-rep(0,m.imp.val)
  slope.se<-rep(0,m.imp.val)
  HarrellC<-rep(0,m.imp.val)
  HarrellC.se<-rep(0,m.imp.val)
  UnoC<-rep(0, m.imp.val)
  UnoC.se<-rep(0, m.imp.val)
  AUC<-rep(0, m.imp.val)
  AUC.se<-rep(0, m.imp.val)

  p.groups<-array(rep(0,g*m.imp.val),dim=c(m.imp.val,g),dimnames=list(1:m.imp.val,1:g))
  y.groups<-array(rep(0,2*g*m.imp.val),dim=c(m.imp.val,g,2),dimnames=list(1:m.imp.val,1:g,c("obs","se")))

  for (i in 1:m.imp.val){
    lp.val<-lp[,i]
    f.val<-survival::coxph(y~lp.val)
    f.val.offset<-survival::coxph(y~offset(lp.val))

    if (show.metrics[2]){
      # calibration intercept
      sf<-survival::survfit(f.val.offset,conf.type="log-log")
      log.H<-log(-log(utils::tail(sf$surv,1)))
      log.H.upper<-log(-log(utils::tail(sf$upper,1)))
      int[i]<-log.H-mean(f.val.offset$linear.predictors)
      int.se[i]<-(log.H-log.H.upper)/stats::qnorm(.975)
    }

    if (show.metrics[3]){
      # calibration slope
      slope[i]<-f.val$coefficients[[1]]
      slope.se[i]<-sqrt(stats::vcov(f.val)[[1,1]])
    }

    if (show.metrics[4]){
      # Harrell's C-index
      rc.H<-survival::concordance(f.val, timewt="n")
      HarrellC[i]<-rc.H$concordance
      HarrellC.se[i]<-sqrt(rc.H$var)
    }

    if (show.metrics[5]){
      # Uno's C-index
      rc.U<-survival::concordance(f.val, timewt="n/G2")
      UnoC[i]<-rc.U$concordance
      UnoC.se[i]<-sqrt(rc.U$var)
    }

    if (show.metrics[6]){
      # if horizon is same as maximum, add small time
      max.time <- y.orig[,1]==time
      if (sum(max.time)>0){
        y.orig[max.time, 1] <- y.orig[max.time, 1] + 0.01
      }

      # time-dependent ROC with CI
      if (CI.metrics){
        AUC.i <- timeROC::timeROC(T=y.orig[, 1],
                                        delta=y.orig[, 2],
                                        marker=lp.val,
                                        times=time,
                                        cause=1,
                                        iid=TRUE)
        AUC[i] <- AUC.i$AUC[2]
        AUC.CI <- stats::confint(AUC.i, n.sim=n.sim)
        AUC.se[i] <- as.numeric((as.numeric(AUC.CI$CI_AUC[2])/100-as.numeric(AUC.i$AUC[2]))/AUC.CI$C.alpha)
      } else{
        # Time-dependent ROC without CI
        AUC[i] <- timeROC::timeROC(T=y.orig[, 1],
                                   delta=y.orig[, 2],
                                   marker=lp.val,
                                   times=time,
                                   cause=1,
                                   iid=FALSE)$AUC[2]
      }
    }

    p.val<-p[,i]
    quants<-stats::quantile(p.val,(1:(g-1))/g)
    cuts<-cut(p.val,breaks=c(0,quants,1))
    p.groups[i,]<-tapply(p.val,cuts,mean)
    for (j in 1:g){
      sub<-(cuts==levels(cuts)[j])
      if (sum(y[sub,2])>0){
        sf<-survival::survfit(y~1,subset=sub,conf.type="log-log")
        y.groups[i,j,1]<-log(-log(utils::tail(sf$surv,1)))
        y.groups[i,j,2]<-(log(-log(utils::tail(sf$surv,1)))-log(-log(utils::tail(sf$upper,1))))/stats::qnorm(.975)} else {y.groups[i,j,]<- -Inf}
    }
  }

  p.mi<-colMeans(p.groups)
  obs.mi<-rep(0,g)
  obs.mi.lower<-rep(0,g)
  obs.mi.upper<-rep(0,g)
  for (j in 1:g)
  {
    RC<-Rubin.combine(y.groups[,j,1],y.groups[,j,2])
    obs.mi[j]<-1-exp(-exp(RC$est))
    obs.mi.lower[j]<-1-exp(-exp(RC$est+stats::qnorm(.025)*RC$se))
    obs.mi.upper[j]<-1-exp(-exp(RC$est+stats::qnorm(.975)*RC$se))
  }

  graphics::par(mar = c(5,5,2,1))
  graphics::par(mar = c(5,5,2,1))
  plot(lim,lim,type='l',xlab="Predicted probability",ylab="Observed frequency",
       main=main,lwd=1,bty='n')

  if (smoothed.curve){
    # lp of smoothed calibration curve of lp.range
    lp.range<-min(lp)+0:100*(max(lp)-min(lp))/100
    lp.sm<-array(rep(0,101*m.imp.val),dim=c(101,m.imp.val),dimnames=list(1:101,1:m.imp.val))
    lp.sm.se<-lp.sm
    for (i in 1:m.imp.val){
      lp.val<-lp[,i]
      f.val.rcs<-rms::cph(y~rms::rcs(lp.val,5),x=TRUE,y=TRUE)
      surv.sm<-rms::Predict(f.val.rcs,lp.val=lp.range,conf.int = 0.95, conf.type = "simultaneous",time=max(y[,1]))
      lp.sm[,i]<-log(-log(surv.sm$yhat))
      lp.sm.se[,i]<-(log(-log(surv.sm$upper))-log(-log(surv.sm$lower)))/(qnorm(.975)-qnorm(.025))
    }

    # predicted and observed of smoothed curve
    p.sm.mi<-1-exp(-exp(lp.range))
    obs.sm.mi<-rep(0,101)
    obs.sm.mi.lower<-rep(0,101)
    obs.sm.mi.upper<-rep(0,101)
    for (j in 1:101){
      RC<-Rubin.combine(lp.sm[j,],lp.sm.se[j,])
      obs.sm.mi[j]<-1-exp(-exp(RC$est))
      obs.sm.mi.lower[j]<-1-exp(-exp(RC$est+stats::qnorm(.025)*RC$se))
      obs.sm.mi.upper[j]<-1-exp(-exp(RC$est+stats::qnorm(.975)*RC$se))
    }

    # add to plot
    graphics::polygon(x=c(p.sm.mi, rev(p.sm.mi)),
                      y=c(obs.sm.mi.lower, rev(obs.sm.mi.upper)), border = NA,
                      col="Lightgray")
    graphics::lines(p.sm.mi, obs.sm.mi, lwd=2, col="Darkgray")
  }

  if (dist){
    line.bins <- 0.0
    length.seg <- 1
    dist.label <- 0.04
    dist.label2 <- 0.03
    d0lab <- 0
    d1lab <- 1
    cex.d01 <- 0.7
    x <- rowMeans(p)
    bins <- seq(0, min(1,max(lim[2])), length = 101)
    x <- x[x >= 0 & x <= 1]
    f0	<-table(cut(x,bins))
    j0	<-f0 > 0
    bins0 <-(bins[-101])[j0]
    f0	<-f0[j0]
    maxf <-max(f0)
    f0	<-(0.1*f0)/maxf
    graphics::segments(bins0, line.bins, bins0, length.seg*f0+line.bins, col="grey")
  }

  graphics::lines(lim, lim)
  graphics::abline(v=quants, col="darkgrey", lwd=1, lty=2)

  graphics::segments(p.mi, obs.mi.lower, p.mi, obs.mi.upper)
  graphics::points(p.mi, obs.mi, pch=20)

  int.mi<-Rubin.combine(int, int.se)
  slope.mi<-Rubin.combine(slope, slope.se)
  HarrellC.mi<-Rubin.combine(HarrellC, HarrellC.se)
  UnoC.mi<-Rubin.combine(UnoC, UnoC.se)
  if(CI.metrics){
    AUC.mi<-Rubin.combine(AUC, AUC.se)
  } else{
    AUC.mi <- list()
    AUC.mi$est <- mean(AUC)
    AUC.mi$se <- 0
  }

  legend.text <- c(paste("N =",format(n,big.mark=",")),
                   paste0("Intercept = ",format(round(int.mi$est,2),nsmall=2),
                          ifelse(CI.metrics,
                                 paste0(" [", format(round(int.mi$est+stats::qnorm(.025)*int.mi$se, 2), nsmall=2),
                                        "; ", format(round(int.mi$est+stats::qnorm(.975)*int.mi$se, 2), nsmall=2), "]"),
                                 "")),
                   paste0("Slope = ",format(round(slope.mi$est,2),nsmall=2),
                          ifelse(CI.metrics,
                                 paste0(" [", format(round(slope.mi$est+stats::qnorm(.025)*slope.mi$se, 2), nsmall=2),
                                        "; ", format(round(slope.mi$est+stats::qnorm(.975)*slope.mi$se, 2), nsmall=2), "]"),
                                 "")),
                   paste0("Harrell's C-index = ",format(round(HarrellC.mi$est-optimism.C,2),nsmall=2),
                          ifelse(CI.metrics,
                                 paste0(" [", format(round(HarrellC.mi$est+stats::qnorm(.025)*HarrellC.mi$se-optimism.C, 2), nsmall=2),
                                        "; ", format(round(HarrellC.mi$est+stats::qnorm(.975)*HarrellC.mi$se-optimism.C, 2), nsmall=2), "]"),
                                 "")),
                   paste0("Uno's C-index = ", format(round(UnoC.mi$est,2),nsmall=2),
                          ifelse(CI.metrics,
                                 paste0(" [", format(round(UnoC.mi$est+stats::qnorm(.025)*UnoC.mi$se-optimism.C, 2), nsmall=2),
                                        "; ", format(round(UnoC.mi$est+stats::qnorm(.975)*UnoC.mi$se-optimism.C, 2), nsmall=2), "]"),
                                 "")),
                   paste0("AUC = ", format(round(AUC.mi$est,2),nsmall=2),
                          ifelse(CI.metrics,
                                 paste0(" [", format(round(AUC.mi$est+stats::qnorm(.025)*AUC.mi$se, 2), nsmall=2),
                                        "; ", format(round(AUC.mi$est+stats::qnorm(.975)*AUC.mi$se, 2), nsmall=2), "]"),
                                 "")))
  if (sum(show.metrics)>0){
    graphics::legend(lim[1], lim[2], legend.text[show.metrics],
                   box.col="white",  bg = "white",cex=1)
  }

  return(list(main=main,
              n=n,quants=quants,
              p.mi=p.mi,
              obs.mi=obs.mi,
              obs.mi.lower=obs.mi.lower,
              obs.mi.upper=obs.mi.upper,
              int=int.mi$est,
              int.lower=int.mi$est+stats::qnorm(.025)*int.mi$se,
              int.upper=int.mi$est+stats::qnorm(.975)*int.mi$se,
              slope=slope.mi$est,
              slope.lower=slope.mi$est+stats::qnorm(.025)*slope.mi$se,
              slope.upper=slope.mi$est+stats::qnorm(.975)*slope.mi$se,
              HarrellC=HarrellC.mi$est,
              HarrellC.se=HarrellC.mi$se,
              HarrellC.lower=HarrellC.mi$est+stats::qnorm(.025)*HarrellC.mi$se,
              HarrellC.upper=HarrellC.mi$est+stats::qnorm(.975)*HarrellC.mi$se,
              UnoC=UnoC.mi$est,
              UnoC.se=UnoC.mi$se,
              UnoC.lower=UnoC.mi$est+stats::qnorm(.025)*UnoC.mi$se,
              UnoC.upper=UnoC.mi$est+stats::qnorm(.975)*UnoC.mi$se,
              AUC=AUC.mi$est,
              AUC.se=AUC.mi$se,
              AUC.lower=AUC.mi$est+stats::qnorm(.025)*AUC.mi$se,
              AUC.upper=AUC.mi$est+stats::qnorm(.975)*AUC.mi$se))
}
