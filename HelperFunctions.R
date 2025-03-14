###
### DEFINE FUNCTIONS
###
measures <- function(predictions=NULL, labels=NULL, N.measures=200, min.pred=0){
  # initialize cutoff values
  cutoffs <- seq(1, min.pred, length.out=N.measures)
  
  # calculate confusion matrix
  TP <- c()
  FP <- c()
  TN <- c()
  FN <- c()
  classified.positive <- c()
  classified.negative <- c()
  for (cutoff in cutoffs){
    TP <- c(TP, sum(predictions>=cutoff & labels==1)) # TODO: what to do with ties?
    FP <- c(FP, sum(predictions>=cutoff & labels==0))
    TN <- c(TN, sum(predictions<cutoff & labels==0))
    FN <- c(FN, sum(predictions<cutoff & labels==1))
    classified.positive <- c(classified.positive, sum(predictions>=cutoff)/length(predictions))
    classified.negative <- c(classified.negative, sum(predictions<cutoff)/length(predictions))
  }
  
  # calculate measures
  SENS <- TP/(TP+FN)
  SPEC <- TN/(TN+FP)
  PPV <- TP/(TP+FP)
  NPV <- TN/(TN+FN)
  
  return(list(cutoffs=cutoffs, SENS=SENS, SPEC=SPEC, PPV=PPV, NPV=NPV,
              classified.positive=classified.positive,
              classified.negative=classified.negative))
}

measures.with.CI <- function(predictions=NULL, labels=NULL, nboot=1,
                             N.measures=200, min.pred=0){
  # calculate measures on original data set
  original <- measures(predictions=predictions, labels=labels,
                       N.measures=N.measures, min.pred=min.pred)
  N.measures <- length(original$SENS)
  
  # initialize data
  df <- data.frame(predictions=predictions, labels=labels)
  N <- nrow(df)
  SENS.boot <- matrix(NA, N.measures, nboot)
  SPEC.boot <- matrix(NA, N.measures, nboot)
  PPV.boot <- matrix(NA, N.measures, nboot)
  NPV.boot <- matrix(NA, N.measures, nboot)
  for (n in 1:nboot){
    # take a bootstrap sample of predictions and labels
    bootstrapped.df <- df[sample(1:N, replace=TRUE),]
    
    # calculate measures on bootstrap samples
    boot.measures <- measures(predictions = bootstrapped.df$predictions,
                              labels = bootstrapped.df$labels,
                              N.measures=N.measures,
                              min.pred=min.pred)
    
    # save measures
    SENS.boot[, n] <- boot.measures$SENS
    SPEC.boot[, n] <- boot.measures$SPEC
    PPV.boot[, n] <- boot.measures$PPV
    NPV.boot[, n] <- boot.measures$NPV
  }
  
  # calculate the lower and upper limit of the confidence interval
  SENS.lower <- apply(SENS.boot, 1, function(x) quantile(x, 0.025, na.rm=TRUE))
  SENS.upper <- apply(SENS.boot, 1, function(x) quantile(x, 0.975, na.rm=TRUE))
  SPEC.lower <- apply(SPEC.boot, 1, function(x) quantile(x, 0.025, na.rm=TRUE))
  SPEC.upper <- apply(SPEC.boot, 1, function(x) quantile(x, 0.975, na.rm=TRUE))
  PPV.lower <- apply(PPV.boot, 1, function(x) quantile(x, 0.025, na.rm=TRUE))
  PPV.upper <- apply(PPV.boot, 1, function(x) quantile(x, 0.975, na.rm=TRUE))
  NPV.lower <- apply(NPV.boot, 1, function(x) quantile(x, 0.025, na.rm=TRUE))
  NPV.upper <- apply(NPV.boot, 1, function(x) quantile(x, 0.975, na.rm=TRUE))
  
  return(list(cutoffs=original$cutoffs,
              SENS=original$SENS,
              SENS.lower=SENS.lower, SENS.upper=SENS.upper,
              SPEC=original$SPEC,
              SPEC.lower=SPEC.lower, SPEC.upper=SPEC.upper,
              PPV=original$PPV,
              PPV.lower=PPV.lower, PPV.upper=PPV.upper,
              pred.pos=original$classified.positive,
              NPV=original$NPV,
              NPV.lower=NPV.lower, NPV.upper=NPV.upper,
              pred.neg=original$classified.negative))
}

create.figure <- function(measures=NULL, threshold=0.5, predictions=NULL,
                          min.pred=0, title=title,
                          file.path=file.path){
  # initialize figure
  grDevices::png(file=paste0(file.path, "Metrics ", 
                             title, '.png'),
                 width=700, height=1000, units="px")
  par(mfrow=c(5, 1), mar=c(4, 4, 1, 1), cex=1)
  
  # Sensitivity
  plot(measures$SENS~measures$cutoffs, type="l",
       xlab="", xlim=c(min.pred, 1),
       ylab="Sensitivity", ylim=c(0, 1), axes=FALSE)
  abline(h=seq(min.pred, 1, 0.125), lty=1, col="grey")
  abline(v=seq(min.pred, 1, 0.125), lty=1, col="grey")
  rect(threshold, 0, 1, 1, col=rgb(1.0, 0, 0, alpha=0.2), border=NA)
  polygon(x=c(measures$cutoffs, rev(measures$cutoffs)),
          y=c(measures$SENS.upper, rev(measures$SENS.lower)),
          col="grey", border=NA)
  lines(measures$cutoffs, measures$SENS, lwd=2)
  axis(side=1, at=seq(min.pred, 1, 0.125), labels=seq(min.pred*100, 100, 12.5)) # x-axis
  axis(side=2, at=seq(0, 1, 0.125), labels=seq(0, 100, 12.5)) # y-axis
  box(col="black")
  
  # Specificity
  plot(measures$SPEC~measures$cutoffs, type="l",
       xlab="", xlim=c(min.pred, 1),
       ylab="Specificity", ylim=c(0, 1), axes=FALSE)
  abline(h=seq(min.pred, 1, 0.125), lty=1, col="grey")
  abline(v=seq(min.pred, 1, 0.125), lty=1, col="grey")
  rect(threshold, 0, 1, 1, col=rgb(1.0, 0, 0, alpha=0.2), border=NA)
  polygon(x=c(measures$cutoffs, rev(measures$cutoffs)),
          y=c(measures$SPEC.upper, rev(measures$SPEC.lower)),
          col="grey", border=NA)
  lines(measures$cutoffs, measures$SPEC, lwd=2)
  axis(side=1, at=seq(min.pred, 1, 0.125), labels=seq(min.pred*100, 100, 12.5)) # x-axis
  axis(side=2, at=seq(0, 1, 0.125), labels=seq(0, 100, 12.5))
  box(col="black")
  
  # Positive predictive value
  plot(measures$PPV~measures$cutoffs, type="l",
       xlab="", xlim=c(min.pred, 1),
       ylab="PPV", ylim=c(0, 1), axes=FALSE)
  polygon(x=c(measures$cutoffs, rev(measures$cutoffs)),
          y=c(measures$pred.pos, rep(0, length(measures$pred.pos))),
          col="grey", border=NA)
  abline(h=seq(min.pred, 1, 0.125), lty=1, col="grey")
  abline(v=seq(min.pred, 1, 0.125), lty=1, col="grey")
  rect(threshold, 0, 1, 1, col=rgb(1.0, 0, 0, alpha=0.2), border=NA)
  polygon(x=c(measures$cutoffs, rev(measures$cutoffs)),
          y=c(measures$PPV.upper, rev(measures$PPV.lower)),
          col="grey", border=NA)
  lines(measures$cutoffs, measures$PPV, lwd=2)
  axis(side=1, at=seq(min.pred, 1, 0.125), labels=seq(min.pred*100, 100, 12.5)) # x-axis
  axis(side=2, at=seq(0, 1, 0.125), labels=seq(0, 100, 12.5))
  box(col="black")
  
  # Negative predictive value
  plot(measures$NPV~measures$cutoffs, type="l",
       xlab="Threshold", xlim=c(min.pred, 1),
       ylab="NPV", ylim=c(0, 1), axes=FALSE)
  polygon(x=c(measures$cutoffs, rev(measures$cutoffs)),
          y=c(measures$pred.neg, rep(0, length(measures$pred.neg))),
          col="grey", border=NA)
  abline(h=seq(min.pred, 1, 0.125), lty=1, col="grey")
  abline(v=seq(min.pred, 1, 0.125), lty=1, col="grey")
  rect(threshold, 0, 1, 1, col=rgb(1.0, 0, 0, alpha=0.2), border=NA)
  polygon(x=c(measures$cutoffs, rev(measures$cutoffs)),
          y=c(measures$NPV.upper, rev(measures$NPV.lower)),
          col="grey", border=NA)
  lines(measures$cutoffs, measures$NPV, lwd=2)
  axis(side=1, at=seq(min.pred, 1, 0.125), labels=seq(min.pred*100, 100, 12.5)) # x-axis
  axis(side=2, at=seq(0, 1, 0.125), labels=seq(0, 100, 12.5))
  box(col="black")
  
  # Histogram
  hist(predictions, breaks=100,
       main="", ylab="", col="black",
       xlab=paste0(title, " (n=", length(predictions), ")"), axes=FALSE, cex.lab=2)
  grDevices::dev.off()
}

# calculate measures at chosen threshold
print.measures <- function(threshold, predictions, labels){
  TP <- sum(predictions>=threshold & labels==1)
  FP <- sum(predictions>=threshold & labels==0)
  TN <- sum(predictions<threshold & labels==0)
  FN <- sum(predictions<threshold & labels==1)
  SENS <- TP/(TP+FN)
  SPEC <- TN/(TN+FP)
  PPV <- TP/(TP+FP)
  NPV <- TN/(TN+FN)
  n <- length(predictions)
  p <- sum(labels)/n
  return(c(n, paste0(sprintf("%.1f", p*100), "%"), 
           paste0(sprintf("%.1f", SENS*100), "%"), 
           paste0(sprintf("%.1f", SPEC*100), "%"), 
           paste0(sprintf("%.1f", PPV*100), "%"), 
           paste0(sprintf("%.1f", NPV*100), "%")))
}