# PredictionTools
This package is under development and for personal use only. It includes various metrics that can be used to validate predictions. Please contact me at carolienm@live.nl if you encounter any issues.

# Install
library(remotes)

remotes::install_github("CHMMaas/PredictionTools")

Note, if you want to update the package you can also use the above call.

# Updates
Release: 2023

Update 25 September 2023: added features
1. indicate which metrics to show in calibration plot with show.metrics (TRUE/FALSE vector of length 8);
2. pre-specify axis limits.

# Load
library(PredictionTools)

# val.surv.mi
This function outputs the calibration and discrimination of predictions of survival of multiple imputed data sets.

# val.prob.mi
This function outputs the calibration and discrimination of predictions of outcome probabilities of multiple imputed data sets.

# Rubin.combine
This function combines multiple estimates with multiples standard errors using the Rubin's rule.

# mb.c
This function calculates model-based concordance.
