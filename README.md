# PredictionTools
This package is under development and for personal use only. It includes various metrics that can be used to validate predictions.

# Install
library(remotes)

remotes::install_github("CHMMaas/PredictionTools")

Note, if you want to update the package you can also use the above call.

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
