% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/val_surv_mi.R
\name{val.surv.mi}
\alias{val.surv.mi}
\title{val.surv.mi}
\usage{
val.surv.mi(
  p,
  y,
  g = 5,
  time = NULL,
  main = "",
  lim = c(0, 1),
  dist = TRUE,
  smoothed.curve = FALSE,
  CI.metrics = FALSE,
  show.metrics = rep(TRUE, 7),
  optimism.C = 0,
  n.sim = 2000
)
}
\arguments{
\item{p}{Matrix with predicted probabilities for imputation i in columns (complete case analysis: one column)}

\item{y}{Time to event outcome as Surv object (time,status), either unrestricted or restricted follow-up}

\item{g}{Number of risk groups; default=5}

\item{time}{Time point at which to evaluate the predicted probabilities, default=NULL (not entered), the maximum time point will be taken. Please note that AUC doesn't compute at maximum follow-up time, you can use show.metrics to omit these results from the plot.}

\item{main}{Plot label, default=""}

\item{lim}{limit, default=NULL}

\item{dist}{distribution, default=TRUE}

\item{smoothed.curve}{plot smoothed calibration curve with 95 percent confidence interval, default=FALSE}

\item{CI.metrics}{plot confidence intervals of calibration intercept, calibration slope, and Harrell's C-index, Uno's C-index, and the area under the time-dependent ROC curve (AUC), default=FALSE.}

\item{show.metrics}{TRUE/FALSE vector of length 6 indicating if plot should show (1) sample size, (2) calibration intercept, (3) calibration slope, (4) Integrated Calibration Index, (5) Harrell's C-index possibly corrected with optimism specified in optimism.C, (6) Uno's C-index, (7) the area under time-dependent ROC curve (AUC) as defined by Blanche et al., default=rep(TRUE, 7)}

\item{optimism.C}{optimism-correction for Harrel's C-index in plot, default=0}

\item{n.sim}{number of simulations for computing simultaneous confidence bands for the time-dependent AUC, which heavily influence computation time, default=2000}
}
\value{
The output of the val_surv_mi function is a "list" with the following components.

main

Main title of plots.


n

number of observations


quants

quantiles


p.mi

predicted survival for each imputed data set.


obs.mi

observed survival for each imputed data set.


obs.mi.lower

lower bound of 95\% confidence interval of observed survival.


obs.mi.upper

upper bound of 95\% confidence interval of observed survival.


int

intercept


int.lower

lower bound of 95\% confidence interval of intercept.


int.upper

upper bound of 95\% confidence interval of intercept.


slope

slope estimate


slope.lower

lower bound of 95\% confidence interval of slope.


slope.upper

upper bound of 95\% confidence interval of slope.


ICI

Integrated Calibration Curve.


HarrellC

Harrell's C-index, uncorrected for optimism.


HarrellC.se

standard error of Harrell's C-index


HarrellC.lower

lower bound of 95\% confidence interval of Harrell's C-index

HarrellC.upper

upper bound of 95\% confidence interval of Harrell's C-index.


UnoC

Uno's C-index for time-to-event outcomes, uncorrected for optimism


UnoC.se

standard error of Uno's C-index, based one a single imputation


UnoC.lower

lower bound of 95\% confidence interval of Uno's C-index


UnoC.lower

lower bound of 95\% confidence interval of Uno's C-index


AUC

time-dependent AUC, calculated using inverse propensity score weighting by the timeROC package, uncorrected for optimism


AUC.se

standard error of time-dependent AUC, based one a single imputation


AUC.lower

lower bound of 95\% confidence interval of time-dependent AUC


AUC.upper

upper bound of 95\% confidence interval of time-dependent AUC
}
\description{
This function calculates intercept, slope, and C-index for survival risk predictions of multiple imputed data set(s).
}
\examples{
library(PredictionTools)
library(survival)
library(rms)
library(mice)
set.seed(1)

# load data
data(pbc, package="survival")

# Multiple imputation
m <- 5
candidate.predictors <- c("sex", "age", "bili", "ascites", "hepato")
imp <- mice::mice(pbc[, candidate.predictors], m=m, print=FALSE)

# set survival time
S <- survival::Surv(time=pbc$time/365, event=as.numeric(pbc$status!=0))

# limit survival to time horizon
horizon <- 5
S.5 <- S
S.5[S[,1]>horizon, 1] <- horizon # if time is more than horizon, set to horizon
S.5[S[,1]>horizon, 2] <- 0       # if time is more than horizon, censor

p <- c()
for (i in 1:m){
 # fit model
  model.i <- rms::cph(S.5 ~ sex + age + bili + ascites + hepato,
                 data=mice::complete(imp, i),
                 x=TRUE, y=TRUE, se.fit=TRUE)

  # obtain linear predictor
  lp.i <- predict(model.i, type="lp")

  # obtain baseline hazard
  f.basehaz.i <- survival::basehaz(model.i)

  # obtain baseline hazard at horizon
  h0.i <- f.basehaz.i$hazard[f.basehaz.i$time==max(f.basehaz.i$time[f.basehaz.i$time<=horizon])]

  # make predictions at time horizon
  p.i <- PredictionTools::fun.event(lp=lp.i, h0=h0.i)

  # add to p
  p <- cbind(p, p.i)
}
# internally validate predictions
g <- 4
main <- paste("Calibration plot for predictions at time",  horizon)
show.metrics <- rep(TRUE, 7)
CI.metrics <- TRUE
PredictionTools::val.surv.mi(p=as.matrix(p), y=S,
                             g=g, time=horizon,
                             main=main, lim=c(0, 1),
                             dist=TRUE, smoothed.curve=TRUE,
                             CI.metrics=CI.metrics,
                             show.metrics=show.metrics,
                             optimism.C=0,
                             n.sim=100)
}
