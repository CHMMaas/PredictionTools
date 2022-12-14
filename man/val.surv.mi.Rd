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
  main = "",
  time = NULL,
  lim = c(0, 1),
  dist = TRUE,
  CI = FALSE,
  df = 3
)
}
\arguments{
\item{p}{Matrix with predicted probabilities for imputation i in columns (complete case analysis: one column)}

\item{y}{Time to event outcome as Surv object (time,status)}

\item{g}{Number of risk groups; default=5}

\item{main}{Plot label, default=""}

\item{time}{Time point at which to evaluate the predicted probabilities, default=NULL (not entered), the maximum time point will be taken}

\item{lim}{limit, default=NULL}

\item{dist}{distribution, default=TRUE}

\item{CI}{plot confidence interval, default=FALSE}

\item{df}{degrees of freedom to compute confidence interval, default=3}
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

lower bound of 95% confidence interval of observed survival.


obs.mi.upper

upper bound of 95% confidence interval of observed survival.


int

intercept


int.lower

lower bound of 95% confidence interval of intercept.


int.upper

upper bound of 95% confidence interval of intercept.


slope

slope estimate


slope.lower

lower bound of 95% confidence interval of slope.


slope.upper

upper bound of 95% confidence interval of slope.


cindex

C-index


cindex.lower

lower bound of 95% confidence interval of C-index.

cindex.upper

upper bound of 95% confidence interval of C-index.
}
\description{
This function calculates intercept, slope, and C-index for survival risk predictions of multiple imputed data set(s).
}
\examples{
library(PredictionTools)
library(survival)
set.seed(1)
n <- 100
m <- 5 # number of imputations
p <- matrix(runif(n*m, 0, 1), n, m)
time_until_event <- runif(n, 0, 5)
status <- rbinom(n, 1, 0.5)
y <- Surv(time_until_event, status)
g <- 4
main <- "Plot label"
time <- 30
PredictionTools::val.surv.mi(p=p, y=y, g=g, main=main, time=time)
}
