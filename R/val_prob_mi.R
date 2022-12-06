#' @title val.prob.mi
#' @description This function calculates intercept, slope, and C-index for outcome risk predictions of multiple imputed data set(s).
#'
#' @importFrom rms Predict
#' @importFrom rms cph
#' @importFrom rms rcs
#' @importFrom Hmisc rcorr.cens
#' @importFrom grDevices dev.off
#' @importFrom grDevices dev.new
#' @importFrom grDevices savePlot
#' @importFrom grDevices png
#' @importFrom graphics abline
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics points
#' @importFrom graphics polygon
#' @importFrom graphics segments
#' @importFrom graphics text
#' @importFrom stats glm
#' @importFrom stats plogis
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' @importFrom stats loess
#' @importFrom stats predict
#' @importFrom stats vcov
#' @importFrom utils tail
#'
#' @param lp.mi Matrix with linear predictor for imputation i in columns (complete case analysis: one column)
#' @param y Outcome indicator, 1 if event, 0 otherwise
#' @param g Number of risk groups; default=5
#' @param main Plot label, default=""
#' @param dist distribution, default=TRUE
#'
#' @return The output of the val_prob_mi function is a "list" with the following components.
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
#' predicted probabilities for each imputed data set.
#'
#'
#' obs.mi
#'
#' observed probabilities for each imputed data set.
#'
#'
#' obs.mi.lower
#'
#' lower bound of 95% confidence interval of observed probabilities.
#'
#'
#' obs.mi.upper
#'
#' upper bound of 95% confidence interval of observed probabilities.
#'
#'
#' int
#'
#' intercept
#'
#'
#' int.lower
#'
#' lower bound of 95% confidence interval of intercept.
#'
#'
#' int.upper
#'
#' upper bound of 95% confidence interval of intercept.
#'
#'
#' slope
#'
#' slope estimate
#'
#'
#' slope.lower
#'
#' lower bound of 95% confidence interval of slope.
#'
#'
#' slope.upper
#'
#' upper bound of 95% confidence interval of slope.
#'
#'
#' cindex
#'
#' C-index
#'
#'
#' cindex.lower
#'
#' lower bound of 95% confidence interval of C-index.
#'
#' cindex.upper
#'
#' upper bound of 95% confidence interval of C-index.
#'
#'
#' mb.c
#'
#' Model-based concordance.
#'
#'
#' e.avg
#'
#' E-average.
#'
#'
#' e.90
#'
#' E-90
#'
#' @export
#'
#' @examples
#' library(PredictionTools)
#' set.seed(1)
#' n <- 100
#' m <- 5 # number of imputations
#' lp.val <- matrix(rnorm(n*m, 0, 1), n, m)
#' y.val <- rbinom(n, 1, 0.5)
#' g <- 4
#' main <- "Plot label"
#' dist <- TRUE
#' PredictionTools::val.prob.mi(lp.mi=lp.val, y=y.val, g=g, main=main, dist=dist)
val.prob.mi<-function(lp.mi, y, g=5, main="", dist=FALSE){
  stopifnot("lp.mi must be numeric" = is.numeric(lp.mi))
  stopifnot("y must be numeric" = is.numeric(y))
  stopifnot("g must be numeric" = is.numeric(g))
  stopifnot("dist must be a boolean (TRUE or FALSE)" = isTRUE(dist)|isFALSE(dist))

  stopifnot("lp.mi must be a vector or matrix" = is.vector(lp.mi) | is.matrix(lp.mi))
  stopifnot("y must be a vector" = is.vector(y))

  stopifnot("lp.mi and y must be the same length" = nrow(lp.mi)==length(y) | length(lp.mi)==length(y))

  # initialize number of observations
  n <- length(y)

  # initialize number of imputations
  m.imp.val<-ncol(lp.mi)
  if (is.null(m.imp.val)){
    m.imp.val <- 1
  }

  # initialize statistic vectors
  cindex<-rep(0,m.imp.val)
  cindex.se<-rep(0,m.imp.val)
  slope<-rep(0,m.imp.val)
  slope.se<-rep(0,m.imp.val)
  int<-rep(0,m.imp.val)
  int.se<-rep(0,m.imp.val)
  cindex<-rep(0,m.imp.val)
  cindex.se<-rep(0,m.imp.val)
  E.avg<-rep(0,m.imp.val)
  E.90<-rep(0,m.imp.val)
  mbc<-rep(0,m.imp.val)
  sm.y<-NULL
  p.groups<-array(rep(0,g*m.imp.val),dim=c(m.imp.val,g),dimnames=list(1:m.imp.val,1:g))
  y.groups<-array(rep(0,2*g*m.imp.val),dim=c(m.imp.val,g,2),dimnames=list(1:m.imp.val,1:g,c("obs","se")))

  # calculate metrics for each imputation
  for (i in 1:m.imp.val)
  {
    if (m.imp.val!=1){
      lp.val <- lp.mi[,i]
    }
    else{
      lp.val <- lp.mi
    }

    f.val<-stats::glm(y~lp.val,family='binomial')
    f.val.offset<-stats::glm(y~offset(lp.val),family='binomial')

    # C-index
    rc<-rcorr.cens(lp.val,y)
    cindex[i]<-rc["C Index"]
    cindex.se[i]<-rc["S.D."]/2
    #cindex[i]<-f.val$stats["C"]

    # slope
    slope[i]<-f.val$coefficients[[2]]
    slope.se[i]<-sqrt(vcov(f.val)[[2,2]])

    # intercept
    int[i]<-f.val.offset$coefficients[[1]]
    int.se[i]<-sqrt(vcov(f.val.offset)[[1,1]])

    # create g quantiles of predicted risk
    p.val<-stats::plogis(lp.val)
    quants<-stats::quantile(p.val,(1:(g-1))/g)
    cuts<-cut(p.val,breaks=c(0,quants,1))
    p.groups[i,]<-tapply(p.val,cuts,mean)
    for (j in 1:g) {
      sub<-(cuts==levels(cuts)[j])
      if (sum(y[sub])>0){
        f.0<-stats::glm(y~1,family = 'binomial',subset=sub)
        y.groups[i,j,1]<-f.0$coef
        y.groups[i,j,2]<-sqrt(vcov(f.0)[1,1])} else {y.groups[i,j,]<- -Inf}
    }

    # model-based concordance
    mbc[i]<-mb.c(as.vector(p.val))

    # calculate E-average and E-90
    Sm <- stats::loess(y~p.val,degree=2)
    sm.y<-cbind(sm.y,stats::predict(Sm,newdata=data.frame(p.val=(0:100)/100)))
    E.avg[i]<-mean(abs(Sm$x-Sm$fitted))
    E.90[i]<-stats::quantile(abs(Sm$x-Sm$fitted),.9)
  }

  # histogram of risk distribution
  line.bins <- 0.05
  length.seg <- 1
  dist.label <- 0.04
  dist.label2 <- 0.03
  d0lab <- 0
  d1lab <- 1
  cex.d01 <- 0.7

  # plot g quantiles by combining over imputations
  p.mi<-colMeans(p.groups)
  obs.mi<-rep(0,g)
  obs.mi.lower<-rep(0,g)
  obs.mi.upper<-rep(0,g)
  for (j in 1:g)
  {
    RC<-Rubin.combine(y.groups[,j,1],y.groups[,j,2])
    obs.mi[j]<-stats::plogis(RC$est)
    obs.mi.lower[j]<-stats::plogis(RC$est+qnorm(.025)*RC$se)
    obs.mi.upper[j]<-stats::plogis(RC$est+qnorm(.975)*RC$se)
  }

  lim<-c(0,1)
  graphics::par(mar = c(5,5,2,1))
  graphics::plot(lim,lim,type='l',xlab="Predicted probability",ylab="Observed frequency",main=main,lwd=1,bty='n')
  graphics::lines(lim,lim)
  graphics::abline(v=quants,col="darkgrey",lwd=1,lty=2)

  if (dist){
    if (m.imp.val > 1){
      x <- rowMeans(stats::plogis(lp.mi))
    } else{
      x <- stats::plogis(lp.mi)
    }
    bins <- seq(0, min(1,max(lim[2])), length = 101)
    x <- x[x >= 0 & x <= 1]
    f0	<-table(cut(x[y==0],bins))
    f1	<-table(cut(x[y==1],bins))
    j0	<-f0 > 0
    j1	<-f1 > 0
    bins0 <-(bins[-101])[j0]
    bins1 <-(bins[-101])[j1]
    f0	<-f0[j0]
    f1	<-f1[j1]
    maxf <-max(f0,f1)
    f0	<-(0.1*f0)/maxf
    f1	<-(0.1*f1)/maxf

    # verticle lines
    graphics::segments(bins1,line.bins,bins1,length.seg*f1+line.bins, col="grey")
    graphics::segments(bins0,line.bins,bins0,length.seg*-f0+line.bins, col="grey")
    # horizontal line
    graphics::lines(c(min(bins0,bins1)-0.01,max(bins0,bins1)+0.01),c(line.bins,line.bins), col="grey")
    # text indicating ones and zeros
    graphics::text(max(bins0,bins1)+dist.label,line.bins+dist.label2,d1lab,cex=cex.d01)
    graphics::text(max(bins0,bins1)+dist.label,line.bins-dist.label2,d0lab,cex=cex.d01)
  }

  sm.y.mi<-rowMeans(sm.y)
  lines((0:100)/100,sm.y.mi,col="dark gray",lwd=2,lty=2)

  graphics::segments(p.mi,obs.mi.lower,p.mi,obs.mi.upper)
  graphics::points(p.mi,obs.mi,pch=20)

  # combine statistics over m imputations
  int.mi<-Rubin.combine(int,int.se)
  slope.mi<-Rubin.combine(slope,slope.se)
  cindex.mi<-Rubin.combine(cindex,cindex.se)
  mbc.mi<-mean(mbc)       #mb.c(stats::plogis(as.vector(lp.mi)))
  E.avg.mi<-mean(E.avg)   ## Standard errors unclear
  E.90.mi<-mean(E.90)     ## Standard errors unclear

  # add statistics to plot
  graphics::legend(lim[1], lim[2], c(paste("n =",format(n,big.mark=",")),
                                     paste("p =",format(sum(y)/n,big.mark=",")),
                           paste("a =",format(round(int.mi$est,2),nsmall=2)),
                           paste("b =",format(round(slope.mi$est,2),nsmall=2)),
                           paste("c =",format(round(cindex.mi$est,2),nsmall=2)),
                           paste("mb.c =",format(round(mbc.mi,2),nsmall=2)),
                           paste("e.avg =",format(round(E.avg.mi,3),nsmall=3),
                                 "(", format(round(E.avg.mi/(sum(y)/n),3),nsmall=3), ")"),
                           paste("e.90 =",format(round(E.90.mi,3),nsmall=3),
                                 "(", format(round(E.90.mi/(sum(y)/n),3),nsmall=3), ")")),
         box.col="white",  bg = "white",cex=1)

  return(list(main=main,
              n=n,
              quants=quants,
              p.mi=p.mi,
              obs.mi=obs.mi,
              obs.mi.lower=obs.mi.lower,
              obs.mi.upper=obs.mi.upper,
              int=int.mi$est,
              int.lower=int.mi$est+qnorm(.025)*int.mi$se,
              int.upper=int.mi$est+qnorm(.975)*int.mi$se,
              slope=slope.mi$est,
              slope.lower=slope.mi$est+qnorm(.025)*slope.mi$se,
              slope.upper=slope.mi$est+qnorm(.975)*slope.mi$se,
              cindex=cindex.mi$est,
              cindex.lower=cindex.mi$est+qnorm(.025)*cindex.mi$se,
              cindex.upper=cindex.mi$est+qnorm(.975)*cindex.mi$se,
              mb.c=mbc.mi,
              e.avg=E.avg.mi,
              e.90=E.90.mi
              ))

}
